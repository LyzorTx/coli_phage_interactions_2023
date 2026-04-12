#!/usr/bin/env python3
"""SX04: Ordinal lysis potency prediction.

LightGBM regression predicting MLC score (0-4) directly instead of binary
classification. Compares nDCG against the binary baseline from SX01.

The SX01 pre-flight confirmed Spearman rho=0.24 between predicted P(lysis)
and MLC grade among positives — the feature space separates potency levels.

Usage:
    python -m lyzortx.pipeline.autoresearch.sx04_eval --device-type cpu
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
from lightgbm import LGBMRegressor
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
from lyzortx.pipeline.autoresearch.spandex_metrics import (
    evaluate_holdout_rows,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bootstrap_spandex_cis,
    load_mlc_scores,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx04_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")

LGBM_REGRESSION_PARAMS = {
    "n_estimators": 300,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "min_child_samples": 10,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
}


def train_ordinal_fold(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    mlc_lookup: dict[tuple[str, str], float],
    seed: int,
    device_type: str,
) -> list[dict[str, object]]:
    """Train LightGBM regressor on MLC 0-4, predict on holdout."""
    # Build entity feature tables.
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

    # Add pairwise cross-terms.
    compute_pairwise_depo_capsule_features(train_design)
    compute_pairwise_depo_capsule_features(holdout_design)
    compute_pairwise_receptor_omp_features(train_design)
    compute_pairwise_receptor_omp_features(holdout_design)

    # Feature columns.
    prefixes = tuple(f"{s}__" for s in host_slots + phage_slots) + (
        "pair_depo_capsule__",
        "pair_receptor_omp__",
    )
    feature_columns = [col for col in train_design.columns if col.startswith(prefixes)]
    categorical_columns = [col for col in (host_categorical + phage_categorical) if col in feature_columns]

    # Build MLC target (ordinal 0-4) for training pairs.
    train_mlc = []
    for _, row in train_design.iterrows():
        key = (str(row["bacteria"]), str(row["phage"]))
        mlc = mlc_lookup.get(key, 0.0)
        train_mlc.append(mlc)
    y_train_mlc = np.array(train_mlc, dtype=float)

    # RFE on binary labels (same feature selection as SX01 for fair comparison).
    y_train_binary = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train_binary, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    LOGGER.info(
        "Ordinal fold: %d features after RFE, %d train, %d holdout",
        len(rfe_features),
        len(train_design),
        len(holdout_design),
    )

    # Train LightGBM regressor on MLC 0-4.
    estimator = LGBMRegressor(
        **LGBM_REGRESSION_PARAMS,
        random_state=seed,
        n_jobs=1,
        verbosity=-1,
        device_type=device_type,
    )
    estimator.fit(
        train_design[rfe_features],
        y_train_mlc,
        sample_weight=sample_weight,
        categorical_feature=rfe_categorical,
    )
    predictions = estimator.predict(holdout_design[rfe_features])
    # Clip to [0, 4] range.
    predictions = np.clip(predictions, 0, 4)

    # Build result rows.
    rows = []
    for row, pred in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        predictions,
    ):
        key = (str(row["bacteria"]), str(row["phage"]))
        mlc = mlc_lookup.get(key)
        rows.append(
            {
                "arm_id": "ordinal_regression",
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(pred) / 4.0),  # normalize to [0,1] for mAP/AUC
                "predicted_mlc": safe_round(float(pred)),
                "mlc_score": mlc,
                "label_binary": int(mlc > 0) if mlc is not None else None,
            }
        )
    return rows


def run_sx04_eval(
    *,
    candidate_module: ModuleType,
    context: Any,
    device_type: str,
    output_dir: Path,
) -> None:
    start_time = datetime.now(timezone.utc)

    # Load data.
    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    ambiguous_pairs = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous_pairs)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean_frame), clean_frame["bacteria"].nunique())

    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    # k-fold CV.
    all_bacteria = sorted(clean_frame["bacteria"].unique())
    fold_assignments = assign_bacteria_folds(all_bacteria)

    all_predictions: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_fold = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_fold = clean_frame[clean_frame["bacteria"].isin(train_bacteria)].copy()

        LOGGER.info("Fold %d: %d train, %d holdout", fold_id, len(training_fold), len(holdout_fold))

        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_ordinal_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=training_fold,
                holdout_frame=holdout_fold,
                mlc_lookup=mlc_lookup,
                seed=seed,
                device_type=device_type,
            )
            fold_rows.extend(rows)

        # Aggregate seeds.
        df = pd.DataFrame(fold_rows)
        agg = df.groupby(
            ["pair_id", "bacteria", "phage", "label_hard_any_lysis", "mlc_score", "label_binary"],
            as_index=False,
        ).agg({"predicted_probability": "mean", "predicted_mlc": "mean"})
        for _, row in agg.iterrows():
            all_predictions.append(dict(row))

        fold_metrics = evaluate_holdout_rows(list(agg.to_dict(orient="records")))
        LOGGER.info(
            "Fold %d: nDCG=%.4f, mAP=%.4f, AUC=%.4f",
            fold_id,
            fold_metrics.get("holdout_ndcg") or 0,
            fold_metrics.get("holdout_map") or 0,
            fold_metrics.get("holdout_roc_auc") or 0,
        )

    # Overall metrics.
    overall = evaluate_holdout_rows(all_predictions)
    LOGGER.info(
        "Overall ordinal regression: nDCG=%.4f, mAP=%.4f, AUC=%.4f, Brier=%.4f",
        overall.get("holdout_ndcg") or 0,
        overall.get("holdout_map") or 0,
        overall.get("holdout_roc_auc") or 0,
        overall.get("holdout_brier_score") or 0,
    )

    # Spearman between predicted MLC and true MLC among positives.
    pos = [r for r in all_predictions if r.get("mlc_score") is not None and float(r["mlc_score"]) > 0]
    if pos:
        true_mlc = [float(r["mlc_score"]) for r in pos]
        pred_mlc = [float(r["predicted_mlc"]) for r in pos]
        rho, pval = spearmanr(true_mlc, pred_mlc)
        LOGGER.info("Spearman (predicted vs true MLC among positives): rho=%.4f, p=%.2e", rho, pval)
    else:
        rho = None

    # Bootstrap CIs.
    LOGGER.info("Computing bootstrap CIs...")
    bootstrap_results = bootstrap_spandex_cis(
        all_predictions,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )
    LOGGER.info("Bootstrap CIs:")
    for metric, ci in bootstrap_results.items():
        LOGGER.info("  %s: %.4f [%.4f, %.4f]", metric, ci.point_estimate or 0, ci.ci_low or 0, ci.ci_high or 0)

    # Write outputs.
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_predictions).to_csv(output_dir / "ordinal_predictions.csv", index=False)
    bootstrap_json = {
        metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
        for metric, ci in bootstrap_results.items()
    }
    bootstrap_json["spearman_rho_positives"] = rho
    with open(output_dir / "bootstrap_results.json", "w", encoding="utf-8") as f:
        json.dump(bootstrap_json, f, indent=2)

    elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
    LOGGER.info("SX04 completed in %.0fs", elapsed)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    LOGGER.info("SX04 eval starting at %s", datetime.now(timezone.utc).isoformat())

    candidate_module = load_module_from_path("sx04_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)

    run_sx04_eval(
        candidate_module=candidate_module,
        context=context,
        device_type=args.device_type,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
