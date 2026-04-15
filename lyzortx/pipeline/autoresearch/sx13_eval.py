#!/usr/bin/env python3
"""SX13: Evaluate host OMP k-mer features (marginal, cross-term, fallback) against SX10 baseline.

Four arms, all evaluated via SX01 10-fold bacteria-stratified CV with bootstrap CIs:

  1. baseline             — SX10 configuration (reference)
  2. marginal             — + host_omp_kmer slot on host side
  3. cross_term           — + host_omp_kmer + within-fold top-100 phage_moriniere × host_omp_kmer pairs
  4. path1_cluster        — + host_omp_cluster slot (categorical) on host side, fallback per spec

Acceptance gate (per ticket): within-panel AUC ≥+2 pp OR Arm C AUC ≥+2 pp over SX10 baseline,
AND NILS53 / Dhillonvirus narrow-host rank improves measurably. Honest null is a valued outcome.

Cross-term arm leakage policy: the top-100 (phage_kmer, host_omp_kmer) pairs are selected by
correlation with `any_lysis` computed ONLY on the training fold's pairs. Holdout fold pairs are
never seen during selection.

Usage:
    python -m lyzortx.pipeline.autoresearch.build_host_omp_kmer_slot
    python -m lyzortx.pipeline.autoresearch.build_host_omp_cluster_slot
    python -m lyzortx.pipeline.autoresearch.sx13_eval --device-type cpu
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

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.build_host_omp_cluster_slot import (
    DEFAULT_SLOT_OUTPUT_PATH as DEFAULT_CLUSTER_SLOT_PATH,
    FEATURE_PREFIX as CLUSTER_PREFIX,
)
from lyzortx.pipeline.autoresearch.build_host_omp_kmer_slot import (
    DEFAULT_SLOT_OUTPUT_PATH as DEFAULT_KMER_SLOT_PATH,
    FEATURE_PREFIX as HOST_KMER_PREFIX,
)
from lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot import FEATURE_PREFIX as PHAGE_KMER_PREFIX
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
    bootstrap_spandex_cis,
    enrich_rows_with_mlc,
    load_mlc_scores,
    train_and_predict_fold,
)
from lyzortx.pipeline.autoresearch.sx12_eval import (
    DEFAULT_KMER_SLOT_PATH as DEFAULT_PHAGE_KMER_SLOT_PATH,
    MORINIERE_SLOT_NAME,
    attach_moriniere_slot,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx13_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")

HOST_OMP_KMER_SLOT = "host_omp_kmer"
HOST_OMP_CLUSTER_SLOT = "host_omp_cluster"
ARMS = ("baseline", "marginal", "cross_term", "path1_cluster")
TOP_N_CROSS_TERMS = 100
CROSS_TERM_PREFIX = "pair_motif_cross__"

LGBM_COMMON_PARAMS = {
    "n_estimators": 300,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "min_child_samples": 10,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_lambda": 1.0,
}


def attach_csv_slot(
    *,
    context: Any,
    features_path: Path,
    slot_name: str,
    feature_prefix: str,
    entity_key: str,
) -> int:
    """Attach a CSV-backed feature slot to context.slot_artifacts. Returns the feature column count."""
    from lyzortx.autoresearch.train import SlotArtifact

    if not features_path.exists():
        raise FileNotFoundError(f"Slot CSV not found at {features_path}; build it before running SX13")
    frame = pd.read_csv(features_path)
    feature_cols = tuple(c for c in frame.columns if c.startswith(feature_prefix))
    if entity_key not in frame.columns or not feature_cols:
        raise ValueError(
            f"Slot CSV {features_path} missing entity column '{entity_key}' or columns starting with '{feature_prefix}'"
        )
    artifact = SlotArtifact(
        slot_name=slot_name,
        entity_key=entity_key,
        feature_columns=feature_cols,
        frame=frame,
    )
    context.slot_artifacts[slot_name] = artifact
    LOGGER.info("Attached slot %s: %d %s × %d features", slot_name, len(frame), entity_key, len(feature_cols))
    return len(feature_cols)


def _select_top_cross_term_pairs(
    train_design: pd.DataFrame,
    phage_kmer_cols: list[str],
    host_kmer_cols: list[str],
    *,
    top_n: int = TOP_N_CROSS_TERMS,
) -> list[tuple[str, str]]:
    """Pick top-N (phage_kmer, host_kmer) pairs by absolute Pearson correlation with any_lysis.

    Computed within-fold (training rows only) to avoid leakage. Fully vectorized via two matrix
    multiplies, exploiting that both phage and host k-mer features are binary.

    For binary cross-term X = phage_p * host_h:
      sum(X)   = (P.T @ H)[p, h]
      sum(X*y) = (P.T @ (H * y))[p, h]
      mean(X)  = sum(X) / n
      var(X)   = mean(X) * (1 - mean(X))     # binary
      corr(X, y) = (sum(X*y)/n - mean(X)*mean(y)) / (std(X) * std(y))
    """
    y = train_design["label_any_lysis"].astype(np.float32).to_numpy()
    n = len(y)
    y_mean = float(y.mean())
    y_std = float(y.std())
    if y_std < 1e-9:
        return []

    P = train_design[phage_kmer_cols].astype(np.float32).to_numpy()  # n × n_p
    H = train_design[host_kmer_cols].astype(np.float32).to_numpy()  # n × n_h

    sum_x = P.T @ H  # n_p × n_h, count of co-occurrence
    sum_xy = P.T @ (H * y[:, None])  # n_p × n_h, count of triple-positive
    mean_x = sum_x / n
    # var_x = mean_x * (1 - mean_x) for binary X; clamp tiny negatives from float noise
    var_x = np.clip(mean_x * (1.0 - mean_x), 0.0, None)
    std_x = np.sqrt(var_x)
    with np.errstate(divide="ignore", invalid="ignore"):
        corr = (sum_xy / n - mean_x * y_mean) / (std_x * y_std)
    corr = np.where(np.isfinite(corr), corr, 0.0)
    abs_corr = np.abs(corr).ravel()

    if top_n >= abs_corr.size:
        order = np.argsort(-abs_corr)
    else:
        # argpartition is O(N), then sort the top-N for stable ranking
        idx = np.argpartition(-abs_corr, top_n)[:top_n]
        order = idx[np.argsort(-abs_corr[idx])]

    n_h = len(host_kmer_cols)
    return [(phage_kmer_cols[i // n_h], host_kmer_cols[i % n_h]) for i in order]


def _add_cross_term_columns(design: pd.DataFrame, pairs: list[tuple[str, str]]) -> list[str]:
    """Add cross-term product columns to design in place. Returns list of new column names."""
    new_cols: list[str] = []
    new_data: dict[str, np.ndarray] = {}
    for phage_col, host_col in pairs:
        col = f"{CROSS_TERM_PREFIX}{phage_col[len(PHAGE_KMER_PREFIX) :]}__x__{host_col[len(HOST_KMER_PREFIX) :]}"
        new_data[col] = design[phage_col].astype(float).to_numpy() * design[host_col].astype(float).to_numpy()
        new_cols.append(col)
    new_frame = pd.DataFrame(new_data, index=design.index)
    for col in new_cols:
        design[col] = new_frame[col]
    return new_cols


def train_and_predict_fold_cross_term(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    seed: int,
    device_type: str,
    host_slots: list[str],
    phage_slots: list[str],
) -> list[dict[str, object]]:
    """Custom train/predict for the cross-term arm. Mirrors train_and_predict_fold but injects
    within-fold-selected (phage_kmer × host_omp_kmer) cross-term columns before RFE.
    """
    from lyzortx.autoresearch.per_phage_model import fit_per_phage_models, predict_per_phage
    from lyzortx.pipeline.autoresearch.candidate_replay import temporary_module_attribute

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

    # Within-fold selection of top phage_kmer × host_omp_kmer pairs.
    phage_kmer_cols = [c for c in train_design.columns if c.startswith(PHAGE_KMER_PREFIX)]
    host_kmer_cols = [c for c in train_design.columns if c.startswith(HOST_KMER_PREFIX)]
    LOGGER.info(
        "Cross-term selection: %d phage k-mers × %d host k-mers = %d candidate pairs",
        len(phage_kmer_cols),
        len(host_kmer_cols),
        len(phage_kmer_cols) * len(host_kmer_cols),
    )
    if phage_kmer_cols and host_kmer_cols:
        top_pairs = _select_top_cross_term_pairs(train_design, phage_kmer_cols, host_kmer_cols, top_n=TOP_N_CROSS_TERMS)
        new_train_cols = _add_cross_term_columns(train_design, top_pairs)
        _add_cross_term_columns(holdout_design, top_pairs)
        LOGGER.info("Cross-term: added %d top-correlation product columns", len(new_train_cols))
    else:
        LOGGER.warning("Cross-term arm skipped: missing phage or host k-mer columns")

    prefixes = tuple(f"{s}__" for s in host_slots + phage_slots) + (
        "pair_depo_capsule__",
        "pair_receptor_omp__",
        CROSS_TERM_PREFIX,
    )
    feature_columns = [c for c in train_design.columns if c.startswith(prefixes)]
    categorical_columns = [c for c in (host_categorical + phage_categorical) if c in feature_columns]

    y_train = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    LOGGER.info(
        "Cross-term fold: %d features after RFE (from %d, including %d cross-terms), %d train / %d holdout",
        len(rfe_features),
        len(feature_columns),
        sum(1 for c in rfe_features if c.startswith(CROSS_TERM_PREFIX)),
        len(train_design),
        len(holdout_design),
    )

    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[rfe_features],
        y_train,
        sample_weight=sample_weight,
        categorical_feature=rfe_categorical,
    )
    all_pairs_predictions = estimator.predict_proba(holdout_design[rfe_features])[:, 1]

    host_feature_columns = [
        c
        for c in rfe_features
        if c.startswith(
            ("host_surface__", "host_typing__", "host_stats__", "host_defense__", HOST_KMER_PREFIX, CLUSTER_PREFIX)
        )
    ]
    per_phage_models = fit_per_phage_models(
        train_design, host_feature_columns, device_type=device_type, random_state=seed
    )
    predictions, _ = predict_per_phage(
        per_phage_models, holdout_design, host_feature_columns, all_pairs_predictions, blend_alpha=0.5
    )

    rows = []
    for row, probability in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        predictions,
    ):
        rows.append(
            {
                "arm_id": "cross_term",
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(probability)),
            }
        )
    return rows


def run_arm(
    *,
    arm_id: str,
    candidate_module: ModuleType,
    context: Any,
    clean_frame: pd.DataFrame,
    fold_assignments: dict[str, int],
    mlc_lookup: dict[tuple[str, str], float],
    device_type: str,
    host_slots: list[str],
    phage_slots: list[str],
) -> list[dict[str, object]]:
    """Run one arm across all folds × seeds; return enriched per-pair prediction rows."""
    all_predictions: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_fold = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        train_fold = clean_frame[clean_frame["bacteria"].isin(train_bacteria)].copy()
        LOGGER.info(
            "=== %s Fold %d: %d train bacteria (%d pairs), %d holdout bacteria (%d pairs) ===",
            arm_id,
            fold_id,
            len(train_bacteria),
            len(train_fold),
            len(holdout_bacteria),
            len(holdout_fold),
        )

        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            if arm_id == "cross_term":
                rows = train_and_predict_fold_cross_term(
                    candidate_module=candidate_module,
                    context=context,
                    training_frame=train_fold,
                    holdout_frame=holdout_fold,
                    seed=seed,
                    device_type=device_type,
                    host_slots=host_slots,
                    phage_slots=phage_slots,
                )
            else:
                rows = train_and_predict_fold(
                    candidate_module=candidate_module,
                    context=context,
                    training_frame=train_fold,
                    holdout_frame=holdout_fold,
                    seed=seed,
                    device_type=device_type,
                    host_slots=host_slots,
                    phage_slots=phage_slots,
                )
            for r in rows:
                r["fold_id"] = fold_id
                r["arm_id"] = arm_id
            fold_rows.extend(rows)

        df = pd.DataFrame(fold_rows)
        agg = (
            df.groupby(["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
                "predicted_probability"
            ]
            .mean()
            .sort_values(["bacteria", "phage"])
        )
        enriched = enrich_rows_with_mlc(agg.to_dict(orient="records"), mlc_lookup)
        for r in enriched:
            r["arm_id"] = arm_id
        fold_metrics = evaluate_holdout_rows(enriched)
        LOGGER.info(
            "  %s Fold %d: nDCG=%.4f, mAP=%.4f, AUC=%.4f, Brier=%.4f",
            arm_id,
            fold_id,
            fold_metrics.get("holdout_ndcg") or 0,
            fold_metrics.get("holdout_map") or 0,
            fold_metrics.get("holdout_roc_auc") or 0,
            fold_metrics.get("holdout_brier_score") or 0,
        )
        all_predictions.extend(enriched)
    return all_predictions


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--phage-kmer-slot-path", type=Path, default=DEFAULT_PHAGE_KMER_SLOT_PATH)
    parser.add_argument("--host-kmer-slot-path", type=Path, default=DEFAULT_KMER_SLOT_PATH)
    parser.add_argument("--host-cluster-slot-path", type=Path, default=DEFAULT_CLUSTER_SLOT_PATH)
    parser.add_argument("--arms", type=str, default=",".join(ARMS))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    start = datetime.now(timezone.utc)
    LOGGER.info("SX13 eval starting at %s", start.isoformat())

    arms = tuple(a.strip() for a in args.arms.split(",") if a.strip())
    unknown = set(arms) - set(ARMS)
    if unknown:
        raise ValueError(f"Unknown arms: {sorted(unknown)}. Valid: {ARMS}")

    candidate_module = load_module_from_path("sx13_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)

    # Always attach phage Moriniere k-mers (needed for cross_term arm and harmless otherwise).
    if "cross_term" in arms or "marginal" in arms:
        attach_moriniere_slot(context, args.phage_kmer_slot_path)
    if "marginal" in arms or "cross_term" in arms:
        attach_csv_slot(
            context=context,
            features_path=args.host_kmer_slot_path,
            slot_name=HOST_OMP_KMER_SLOT,
            feature_prefix=HOST_KMER_PREFIX,
            entity_key="bacteria",
        )
    if "path1_cluster" in arms:
        attach_csv_slot(
            context=context,
            features_path=args.host_cluster_slot_path,
            slot_name=HOST_OMP_CLUSTER_SLOT,
            feature_prefix=CLUSTER_PREFIX,
            entity_key="bacteria",
        )

    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    ambiguous = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean_frame), clean_frame["bacteria"].nunique())

    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    fold_assignments = assign_bacteria_folds(sorted(clean_frame["bacteria"].unique()))

    arm_specs = {
        "baseline": {
            "host_slots": ["host_surface", "host_typing", "host_stats", "host_defense"],
            "phage_slots": ["phage_projection", "phage_stats"],
        },
        "marginal": {
            "host_slots": ["host_surface", "host_typing", "host_stats", "host_defense", HOST_OMP_KMER_SLOT],
            "phage_slots": ["phage_projection", "phage_stats"],
        },
        "cross_term": {
            "host_slots": ["host_surface", "host_typing", "host_stats", "host_defense", HOST_OMP_KMER_SLOT],
            "phage_slots": ["phage_projection", "phage_stats", MORINIERE_SLOT_NAME],
        },
        "path1_cluster": {
            "host_slots": ["host_surface", "host_typing", "host_stats", "host_defense", HOST_OMP_CLUSTER_SLOT],
            "phage_slots": ["phage_projection", "phage_stats"],
        },
    }

    args.output_dir.mkdir(parents=True, exist_ok=True)
    arm_summary_rows: list[dict[str, object]] = []
    bootstrap_bundle: dict[str, dict[str, dict[str, float | None]]] = {}

    for arm in arms:
        spec = arm_specs[arm]
        predictions = run_arm(
            arm_id=arm,
            candidate_module=candidate_module,
            context=context,
            clean_frame=clean_frame,
            fold_assignments=fold_assignments,
            mlc_lookup=mlc_lookup,
            device_type=args.device_type,
            host_slots=spec["host_slots"],
            phage_slots=spec["phage_slots"],
        )

        LOGGER.info("=== Bootstrap for arm %s (%d predictions) ===", arm, len(predictions))
        ci_results = bootstrap_spandex_cis(
            predictions, bootstrap_samples=BOOTSTRAP_SAMPLES, bootstrap_random_state=BOOTSTRAP_RANDOM_STATE
        )
        bootstrap_bundle[arm] = {
            metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
            for metric, ci in ci_results.items()
        }
        summary_row: dict[str, object] = {"arm_id": arm}
        for metric, ci in ci_results.items():
            summary_row[metric] = ci.point_estimate
            summary_row[f"{metric}_ci_low"] = ci.ci_low
            summary_row[f"{metric}_ci_high"] = ci.ci_high
        arm_summary_rows.append(summary_row)
        pd.DataFrame(predictions).to_csv(args.output_dir / f"{arm}_predictions.csv", index=False)

    with open(args.output_dir / "bootstrap_results.json", "w", encoding="utf-8") as f:
        json.dump(bootstrap_bundle, f, indent=2)
    pd.DataFrame(arm_summary_rows).to_csv(args.output_dir / "arm_comparison.csv", index=False)

    elapsed = (datetime.now(timezone.utc) - start).total_seconds()
    LOGGER.info("SX13 eval complete in %.0fs", elapsed)


if __name__ == "__main__":
    main()
