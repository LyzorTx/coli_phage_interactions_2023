#!/usr/bin/env python3
"""SX11: Potency loss-function ablation.

Four arms, each evaluated via the same 10-fold bacteria-stratified CV + bootstrap CI protocol
that SX01/SX04 use. All arms share identical feature engineering, RFE feature selection, and
fold assignments; the only variable is the training loss.

Arms:
  1. binary_baseline       — LGBMClassifier on any_lysis (reference, same as SX01 without per-phage blending)
  2. hurdle_two_stage      — LGBMClassifier on any_lysis * LGBMRegressor on MLC conditional on positive
  3. lambdarank            — LightGBM objective='lambdarank', query group=bacterium, relevance=MLC 0-3
  4. ordinal_all_threshold — three binary classifiers y>=1, y>=2, y>=3 combined via cumulative probabilities

Usage:
    python -m lyzortx.pipeline.autoresearch.sx11_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier, LGBMRanker, LGBMRegressor
from scipy.stats import spearmanr

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    build_st03_training_frame,
    load_module_from_path,
    load_st03_holdout_frame,
    safe_round,
)
from lyzortx.pipeline.autoresearch.derive_pairwise_depo_capsule_features import (
    compute_pairwise_depo_capsule_features,
)
from lyzortx.pipeline.autoresearch.derive_pairwise_receptor_omp_features import (
    compute_pairwise_receptor_omp_features,
)
from lyzortx.pipeline.autoresearch.gt03_eval import apply_rfe
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.spandex_metrics import evaluate_holdout_rows
from lyzortx.pipeline.autoresearch.sx01_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
    bootstrap_spandex_cis,
    load_mlc_scores,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx11_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")

LGBM_COMMON_PARAMS = {
    "n_estimators": 300,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "min_child_samples": 10,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
}

MAX_MLC = 3.0  # post-SX05
ARMS = ("binary_baseline", "hurdle_two_stage", "lambdarank", "ordinal_all_threshold")


def _build_fold_design(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """Shared feature-engineering + RFE step used by all arms.

    Returns (train_design, holdout_design, rfe_features, rfe_categorical).
    """
    host_slots = ["host_surface", "host_typing", "host_stats", "host_defense"]
    phage_slots = ["phage_projection", "phage_stats"]

    host_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=host_slots, entity_key="bacteria"
    )
    phage_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=phage_slots, entity_key="phage"
    )
    host_typed, _, host_categorical = candidate_module.type_entity_features(host_table, "bacteria")
    phage_typed, _, phage_categorical = candidate_module.type_entity_features(phage_table, "phage")

    train_design = candidate_module.build_raw_pair_design_matrix(
        training_frame, host_features=host_typed, phage_features=phage_typed
    )
    holdout_design = candidate_module.build_raw_pair_design_matrix(
        holdout_frame, host_features=host_typed, phage_features=phage_typed
    )

    compute_pairwise_depo_capsule_features(train_design)
    compute_pairwise_depo_capsule_features(holdout_design)
    compute_pairwise_receptor_omp_features(train_design)
    compute_pairwise_receptor_omp_features(holdout_design)

    prefixes = tuple(f"{s}__" for s in host_slots + phage_slots) + (
        "pair_depo_capsule__",
        "pair_receptor_omp__",
    )
    feature_columns = [c for c in train_design.columns if c.startswith(prefixes)]
    categorical_columns = [c for c in (host_categorical + phage_categorical) if c in feature_columns]

    y_train_binary = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train_binary, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    return train_design, holdout_design, rfe_features, rfe_categorical


def _train_binary(train_design: pd.DataFrame, features: list[str], cat: list[str], seed: int, device_type: str):
    model = LGBMClassifier(**LGBM_COMMON_PARAMS, random_state=seed, n_jobs=1, verbosity=-1, device_type=device_type)
    y = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    w = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)
    model.fit(train_design[features], y, sample_weight=w, categorical_feature=cat)
    return model


def _predict_binary(model, holdout_design: pd.DataFrame, features: list[str]) -> np.ndarray:
    return model.predict_proba(holdout_design[features])[:, 1]


def _arm_binary_baseline(train_design, holdout_design, features, cat, seed, device_type, mlc_lookup):
    model = _train_binary(train_design, features, cat, seed, device_type)
    preds = _predict_binary(model, holdout_design, features)
    return _result_rows("binary_baseline", holdout_design, preds, mlc_lookup, seed)


def _arm_hurdle_two_stage(train_design, holdout_design, features, cat, seed, device_type, mlc_lookup):
    """P(y>=1) via classifier, E[MLC | y>=1] via regressor on positives only; combine."""
    classifier = _train_binary(train_design, features, cat, seed, device_type)
    p_positive = _predict_binary(classifier, holdout_design, features)

    pos_mask = train_design["label_any_lysis"].astype(int) == 1
    if pos_mask.sum() < 10:
        # Degenerate fold — fall back to binary-only score.
        return _result_rows("hurdle_two_stage", holdout_design, p_positive, mlc_lookup, seed)

    train_pos = train_design.loc[pos_mask]
    mlc_train_pos = np.array(
        [float(mlc_lookup.get((str(r["bacteria"]), str(r["phage"])), 0.0)) for _, r in train_pos.iterrows()],
        dtype=float,
    )
    # Only train if there's variance in positive MLC grades.
    if mlc_train_pos.std() < 1e-6:
        combined = p_positive * float(mlc_train_pos.mean() or 1.0)
    else:
        regressor = LGBMRegressor(
            **LGBM_COMMON_PARAMS, random_state=seed, n_jobs=1, verbosity=-1, device_type=device_type
        )
        regressor.fit(
            train_pos[features],
            mlc_train_pos,
            sample_weight=train_pos["training_weight_v3"].astype(float).to_numpy(dtype=float),
            categorical_feature=cat,
        )
        conditional_mlc = regressor.predict(holdout_design[features])
        conditional_mlc = np.clip(conditional_mlc, 1.0, MAX_MLC)
        combined = p_positive * conditional_mlc

    # Normalize to [0, 1] for mAP/AUC/Brier by dividing by MAX_MLC.
    return _result_rows("hurdle_two_stage", holdout_design, combined / MAX_MLC, mlc_lookup, seed)


def _arm_lambdarank(train_design, holdout_design, features, cat, seed, device_type, mlc_lookup):
    """LightGBM LambdaRank; query group = bacterium, relevance = MLC 0-3."""
    mlc_train = np.array(
        [float(mlc_lookup.get((str(r["bacteria"]), str(r["phage"])), 0.0)) for _, r in train_design.iterrows()],
        dtype=float,
    )
    # LightGBM wants integer-relevance.
    relevance = mlc_train.astype(int)

    # Query groups must be contiguous; sort by bacteria and compute group sizes.
    train_sorted = train_design.assign(_rel=relevance).sort_values(["bacteria", "phage"]).reset_index(drop=True)
    group_sizes = train_sorted.groupby("bacteria", sort=False).size().tolist()

    model = LGBMRanker(
        **LGBM_COMMON_PARAMS,
        objective="lambdarank",
        random_state=seed,
        n_jobs=1,
        verbosity=-1,
        device_type=device_type,
    )
    model.fit(
        train_sorted[features],
        train_sorted["_rel"],
        group=group_sizes,
        categorical_feature=cat,
    )
    # Ranker predictions are unbounded scores; min-max normalize per bacterium to [0, 1] so AUC/mAP/Brier stay meaningful.
    raw = model.predict(holdout_design[features])
    normed = _minmax_per_bacterium(holdout_design["bacteria"].astype(str).to_numpy(), raw)
    return _result_rows("lambdarank", holdout_design, normed, mlc_lookup, seed)


def _arm_ordinal_all_threshold(train_design, holdout_design, features, cat, seed, device_type, mlc_lookup):
    """Fit 3 binary classifiers (y>=1, y>=2, y>=3); combine via cumulative probabilities."""
    mlc_train = np.array(
        [float(mlc_lookup.get((str(r["bacteria"]), str(r["phage"])), 0.0)) for _, r in train_design.iterrows()],
        dtype=float,
    )
    weights = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    threshold_probs = []
    for threshold in (1, 2, 3):
        y_thresh = (mlc_train >= threshold).astype(int)
        if y_thresh.sum() < 10 or (1 - y_thresh).sum() < 10:
            # Not enough positives/negatives for this threshold; fall back to zero probability for this level.
            threshold_probs.append(np.zeros(len(holdout_design)))
            continue
        model = LGBMClassifier(**LGBM_COMMON_PARAMS, random_state=seed, n_jobs=1, verbosity=-1, device_type=device_type)
        model.fit(train_design[features], y_thresh, sample_weight=weights, categorical_feature=cat)
        threshold_probs.append(model.predict_proba(holdout_design[features])[:, 1])

    # Combine: expected MLC = sum(P(y>=k) for k in 1..3); normalize to [0, 1].
    expected_mlc = np.sum(threshold_probs, axis=0)
    normed = np.clip(expected_mlc / MAX_MLC, 0.0, 1.0)
    return _result_rows("ordinal_all_threshold", holdout_design, normed, mlc_lookup, seed)


def _minmax_per_bacterium(bacteria: np.ndarray, scores: np.ndarray) -> np.ndarray:
    out = np.zeros_like(scores, dtype=float)
    for b in np.unique(bacteria):
        mask = bacteria == b
        vals = scores[mask]
        if vals.max() - vals.min() < 1e-12:
            out[mask] = 0.5
        else:
            out[mask] = (vals - vals.min()) / (vals.max() - vals.min())
    return out


def _result_rows(arm_id, holdout_design, preds, mlc_lookup, seed) -> list[dict[str, object]]:
    rows = []
    for row, pred in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        preds,
    ):
        key = (str(row["bacteria"]), str(row["phage"]))
        mlc = mlc_lookup.get(key)
        rows.append(
            {
                "arm_id": arm_id,
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(pred)),
                "mlc_score": float(mlc) if mlc is not None else None,
                "label_binary": int(mlc > 0) if mlc is not None else None,
            }
        )
    return rows


ARM_FUNCTIONS = {
    "binary_baseline": _arm_binary_baseline,
    "hurdle_two_stage": _arm_hurdle_two_stage,
    "lambdarank": _arm_lambdarank,
    "ordinal_all_threshold": _arm_ordinal_all_threshold,
}


def run_sx11(
    *,
    candidate_module: ModuleType,
    context: Any,
    device_type: str,
    output_dir: Path,
    arms: tuple[str, ...] = ARMS,
) -> None:
    start = datetime.now(timezone.utc)
    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    ambiguous = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean_frame), clean_frame["bacteria"].nunique())

    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    fold_assignments = assign_bacteria_folds(bacteria_to_cv_group_map(clean_frame))

    per_arm_predictions: dict[str, list[dict[str, object]]] = {arm: [] for arm in arms}

    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_fold = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_fold = clean_frame[clean_frame["bacteria"].isin(train_bacteria)].copy()

        LOGGER.info(
            "=== Fold %d: %d train bacteria (%d pairs), %d holdout bacteria (%d pairs) ===",
            fold_id,
            len(train_bacteria),
            len(training_fold),
            len(holdout_bacteria),
            len(holdout_fold),
        )

        train_design, holdout_design, rfe_features, rfe_cat = _build_fold_design(
            candidate_module=candidate_module,
            context=context,
            training_frame=training_fold,
            holdout_frame=holdout_fold,
        )

        for arm in arms:
            arm_fn = ARM_FUNCTIONS[arm]
            fold_rows: list[dict[str, object]] = []
            for seed in SEEDS:
                rows = arm_fn(train_design, holdout_design, rfe_features, rfe_cat, seed, device_type, mlc_lookup)
                for r in rows:
                    r["fold_id"] = fold_id
                fold_rows.extend(rows)
            # Aggregate seeds per (pair, arm).
            df = pd.DataFrame(fold_rows)
            agg = df.groupby(
                ["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis", "mlc_score", "label_binary"],
                as_index=False,
            ).agg({"predicted_probability": "mean"})
            agg["arm_id"] = arm
            per_arm_predictions[arm].extend(agg.to_dict(orient="records"))
            fold_metrics = evaluate_holdout_rows(agg.to_dict(orient="records"))
            LOGGER.info(
                "  Fold %d arm %s: nDCG=%.4f, mAP=%.4f, AUC=%.4f, Brier=%.4f",
                fold_id,
                arm,
                fold_metrics.get("holdout_ndcg") or 0,
                fold_metrics.get("holdout_map") or 0,
                fold_metrics.get("holdout_roc_auc") or 0,
                fold_metrics.get("holdout_brier_score") or 0,
            )

    output_dir.mkdir(parents=True, exist_ok=True)
    arm_summary_rows = []
    bootstrap_bundle: dict[str, dict[str, dict[str, float | None]]] = {}
    for arm, preds in per_arm_predictions.items():
        LOGGER.info("=== Bootstrap for arm %s (%d predictions) ===", arm, len(preds))
        bootstrap_results = bootstrap_spandex_cis(
            preds, bootstrap_samples=BOOTSTRAP_SAMPLES, bootstrap_random_state=BOOTSTRAP_RANDOM_STATE
        )
        bootstrap_bundle[arm] = {
            metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
            for metric, ci in bootstrap_results.items()
        }
        summary_row = {"arm_id": arm}
        for metric, ci in bootstrap_results.items():
            summary_row[f"{metric}"] = ci.point_estimate
            summary_row[f"{metric}_ci_low"] = ci.ci_low
            summary_row[f"{metric}_ci_high"] = ci.ci_high
        # Spearman for arms that produce a grade-proportional score.
        positives = [r for r in preds if r.get("mlc_score") and float(r["mlc_score"]) > 0]
        if positives:
            rho, pval = spearmanr(
                [float(r["mlc_score"]) for r in positives],
                [float(r["predicted_probability"]) for r in positives],
            )
            summary_row["spearman_rho_positives"] = float(rho)
            summary_row["spearman_pvalue"] = float(pval)
        arm_summary_rows.append(summary_row)

        pd.DataFrame(preds).to_csv(output_dir / f"{arm}_predictions.csv", index=False)

    with open(output_dir / "bootstrap_results.json", "w", encoding="utf-8") as f:
        json.dump(bootstrap_bundle, f, indent=2)
    pd.DataFrame(arm_summary_rows).to_csv(output_dir / "arm_comparison.csv", index=False)

    elapsed = (datetime.now(timezone.utc) - start).total_seconds()
    LOGGER.info("SX11 completed in %.0fs", elapsed)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--arms", type=str, default=",".join(ARMS))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    LOGGER.info("SX11 eval starting at %s", datetime.now(timezone.utc).isoformat())

    candidate_module = load_module_from_path("sx11_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)

    arms = tuple(a.strip() for a in args.arms.split(",") if a.strip())
    unknown = set(arms) - set(ARMS)
    if unknown:
        raise ValueError(f"Unknown arms requested: {sorted(unknown)}. Valid arms: {ARMS}")

    run_sx11(
        candidate_module=candidate_module,
        context=context,
        device_type=args.device_type,
        output_dir=args.output_dir,
        arms=arms,
    )


if __name__ == "__main__":
    main()
