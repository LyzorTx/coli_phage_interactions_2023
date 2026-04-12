#!/usr/bin/env python3
"""SX01: Graded evaluation framework + clean-label baseline.

Implements the SPANDEX evaluation suite:
  - nDCG with graded MLC 0-4 relevance (primary ranking metric)
  - mAP with binary relevance (secondary ranking metric)
  - AUC and Brier (secondary discrimination/calibration metrics)
  - 10-fold bacteria-stratified CV with bootstrap CIs
  - Clean-label training (exclude ambiguous 'n' pairs)

Pre-flight gate: checks whether existing GT03 predictions separate MLC
grades. If Spearman rho < 0.1 among positives, graded nDCG is cosmetic.

Usage:
    python -m lyzortx.pipeline.autoresearch.sx01_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
from scipy.stats import spearmanr

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    BootstrapMetricCI,
    build_st03_training_frame,
    load_module_from_path,
    load_st03_holdout_frame,
    safe_round,
)
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.spandex_metrics import (
    SPANDEX_METRIC_NAMES,
    evaluate_holdout_rows,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx01_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")
INTERACTION_MATRIX_PATH = Path("data/interactions/interaction_matrix.csv")

N_FOLDS = 10
FOLD_SALT = "spandex_v1"
SEEDS = [7, 42, 123]
BOOTSTRAP_SAMPLES = 1000
BOOTSTRAP_RANDOM_STATE = 42


def load_mlc_scores() -> pd.DataFrame:
    """Load MLC 0-4 scores from interaction matrix, return long-form DataFrame."""
    matrix = pd.read_csv(INTERACTION_MATRIX_PATH, sep=";", index_col=0)
    rows = []
    for bacteria in matrix.index:
        for phage in matrix.columns:
            rows.append(
                {
                    "bacteria": str(bacteria),
                    "phage": str(phage),
                    "mlc_score": float(matrix.loc[bacteria, phage]),
                }
            )
    return pd.DataFrame(rows)


def assign_bacteria_folds(
    bacteria_list: Sequence[str], n_folds: int = N_FOLDS, salt: str = FOLD_SALT
) -> dict[str, int]:
    """Deterministic fold assignment via hashing."""
    assignments = {}
    for bacteria in bacteria_list:
        h = hashlib.sha256(f"{salt}:{bacteria}".encode()).hexdigest()
        assignments[bacteria] = int(h, 16) % n_folds
    return assignments


def run_preflight(output_dir: Path) -> dict[str, object]:
    """Pre-flight: check MLC grade separability in existing GT03 predictions."""
    gt03_path = Path("lyzortx/generated_outputs/gt03_eval/aggregated_predictions.csv")
    if not gt03_path.exists():
        LOGGER.warning("GT03 predictions not found at %s — skipping pre-flight", gt03_path)
        return {"status": "skipped", "reason": "GT03 predictions not found"}

    preds = pd.read_csv(gt03_path)
    best_arm = preds[preds["arm_id"] == "all_gates_rfe"].copy()
    mlc = load_mlc_scores()
    merged = best_arm.merge(mlc, on=["bacteria", "phage"], how="left")

    positives = merged[merged["mlc_score"] >= 1]
    mean_by_grade = {}
    for grade in [0, 1, 2, 3, 4]:
        subset = merged[merged["mlc_score"] == grade]
        if len(subset) > 0:
            mean_by_grade[str(grade)] = round(float(subset["predicted_probability"].mean()), 4)

    rho, pval = spearmanr(positives["mlc_score"], positives["predicted_probability"])
    rho = float(rho)

    result = {
        "status": "pass" if rho >= 0.1 else "fail",
        "spearman_rho_positives": round(rho, 4),
        "spearman_pvalue": float(pval),
        "mean_pred_by_mlc_grade": mean_by_grade,
        "n_positives": len(positives),
        "monotonic": all(mean_by_grade.get(str(i), 0) <= mean_by_grade.get(str(i + 1), 0) for i in range(1, 4)),
    }

    LOGGER.info("Pre-flight: Spearman rho=%.4f (threshold 0.1), status=%s", rho, result["status"])
    for grade, mean_p in sorted(mean_by_grade.items()):
        LOGGER.info("  MLC=%s: mean P(lysis)=%.4f", grade, mean_p)

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "preflight_results.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)

    return result


def enrich_rows_with_mlc(
    rows: Sequence[dict[str, object]],
    mlc_lookup: dict[tuple[str, str], float],
) -> list[dict[str, object]]:
    """Add mlc_score and label_binary to prediction rows."""
    enriched = []
    for row in rows:
        key = (str(row["bacteria"]), str(row["phage"]))
        mlc = mlc_lookup.get(key)
        enriched.append(
            {
                **row,
                "mlc_score": mlc,
                "label_binary": int(mlc > 0) if mlc is not None else None,
            }
        )
    return enriched


def train_and_predict_fold(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    seed: int,
    device_type: str,
) -> list[dict[str, object]]:
    """Train with RFE + per-phage blending, predict on holdout_frame."""
    from lyzortx.autoresearch.per_phage_model import fit_per_phage_models, predict_per_phage
    from lyzortx.pipeline.autoresearch.candidate_replay import temporary_module_attribute
    from lyzortx.pipeline.autoresearch.derive_pairwise_depo_capsule_features import (
        compute_pairwise_depo_capsule_features,
    )
    from lyzortx.pipeline.autoresearch.derive_pairwise_receptor_omp_features import (
        compute_pairwise_receptor_omp_features,
    )
    from lyzortx.pipeline.autoresearch.gt03_eval import apply_rfe

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

    # Build design matrices.
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

    # Identify feature columns.
    prefixes = tuple(f"{s}__" for s in host_slots + phage_slots) + (
        "pair_depo_capsule__",
        "pair_receptor_omp__",
    )
    feature_columns = [col for col in train_design.columns if col.startswith(prefixes)]
    categorical_columns = [col for col in (host_categorical + phage_categorical) if col in feature_columns]

    # Apply RFE.
    y_train = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    LOGGER.info(
        "Fold training: %d features after RFE (from %d), %d train, %d holdout",
        len(rfe_features),
        len(feature_columns),
        len(train_design),
        len(holdout_design),
    )

    # Train all-pairs model.
    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[rfe_features],
        y_train,
        sample_weight=sample_weight,
        categorical_feature=rfe_categorical,
    )
    all_pairs_predictions = estimator.predict_proba(holdout_design[rfe_features])[:, 1]

    # Per-phage blending.
    host_feature_columns = [
        col
        for col in rfe_features
        if col.startswith(("host_surface__", "host_typing__", "host_stats__", "host_defense__"))
    ]
    per_phage_models = fit_per_phage_models(
        train_design, host_feature_columns, device_type=device_type, random_state=seed
    )
    predictions, _ = predict_per_phage(
        per_phage_models, holdout_design, host_feature_columns, all_pairs_predictions, blend_alpha=0.5
    )

    # Build result rows.
    rows = []
    for row, probability in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        predictions,
    ):
        rows.append(
            {
                "arm_id": "spandex_baseline",
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(probability)),
            }
        )
    return rows


def run_kfold_evaluation(
    *,
    candidate_module: ModuleType,
    context: Any,
    full_frame: pd.DataFrame,
    mlc_lookup: dict[tuple[str, str], float],
    device_type: str,
    output_dir: Path,
) -> dict[str, object]:
    """Run 10-fold bacteria-stratified CV with SPANDEX metrics."""
    all_bacteria = sorted(full_frame["bacteria"].unique())
    fold_assignments = assign_bacteria_folds(all_bacteria)

    LOGGER.info("k-fold CV: %d bacteria across %d folds", len(all_bacteria), N_FOLDS)
    fold_sizes = defaultdict(int)
    for fold in fold_assignments.values():
        fold_sizes[fold] += 1
    for fold_id in range(N_FOLDS):
        LOGGER.info("  Fold %d: %d bacteria", fold_id, fold_sizes[fold_id])

    all_predictions: list[dict[str, object]] = []
    fold_metrics: list[dict[str, object]] = []

    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}

        holdout_frame = full_frame[full_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_frame = full_frame[full_frame["bacteria"].isin(train_bacteria)].copy()

        LOGGER.info(
            "=== Fold %d: %d train bacteria (%d pairs), %d holdout bacteria (%d pairs) ===",
            fold_id,
            len(train_bacteria),
            len(training_frame),
            len(holdout_bacteria),
            len(holdout_frame),
        )

        fold_seed_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_and_predict_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=training_frame,
                holdout_frame=holdout_frame,
                seed=seed,
                device_type=device_type,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_seed_rows.extend(rows)

        # Aggregate across seeds.
        df = pd.DataFrame(fold_seed_rows)
        aggregated = (
            df.groupby(["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
                "predicted_probability"
            ]
            .mean()
            .sort_values(["bacteria", "phage"])
        )
        agg_rows = aggregated.to_dict(orient="records")
        enriched = enrich_rows_with_mlc(agg_rows, mlc_lookup)
        all_predictions.extend(enriched)

        metrics = evaluate_holdout_rows(enriched)
        metrics["fold_id"] = fold_id
        metrics["n_holdout_bacteria"] = len(holdout_bacteria)
        metrics["n_holdout_pairs"] = len(holdout_frame)
        fold_metrics.append(metrics)

        LOGGER.info(
            "Fold %d: nDCG=%.4f, mAP=%.4f, AUC=%.4f, Brier=%.4f",
            fold_id,
            metrics.get("holdout_ndcg") or 0,
            metrics.get("holdout_map") or 0,
            metrics.get("holdout_roc_auc") or 0,
            metrics.get("holdout_brier_score") or 0,
        )

    # Overall metrics across all folds.
    overall_metrics = evaluate_holdout_rows(all_predictions)
    LOGGER.info("=" * 60)
    LOGGER.info("Overall k-fold CV metrics:")
    for name, value in overall_metrics.items():
        LOGGER.info("  %s: %.4f", name, value or 0)

    # Bootstrap CIs.
    LOGGER.info("Running bootstrap CIs (%d resamples)...", BOOTSTRAP_SAMPLES)
    bootstrap_results = bootstrap_spandex_cis(
        all_predictions,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    LOGGER.info("Bootstrap CIs:")
    for metric_name, ci in bootstrap_results.items():
        LOGGER.info(
            "  %s: %.4f [%.4f, %.4f]",
            metric_name,
            ci.point_estimate or 0,
            ci.ci_low or 0,
            ci.ci_high or 0,
        )

    # Write outputs.
    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_predictions).to_csv(output_dir / "kfold_predictions.csv", index=False)
    pd.DataFrame(fold_metrics).to_csv(output_dir / "fold_metrics.csv", index=False)

    bootstrap_json = {
        metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
        for metric, ci in bootstrap_results.items()
    }
    with open(output_dir / "bootstrap_results.json", "w", encoding="utf-8") as f:
        json.dump(bootstrap_json, f, indent=2)

    return {
        "overall_metrics": overall_metrics,
        "fold_metrics": fold_metrics,
        "bootstrap_results": bootstrap_results,
    }


def bootstrap_spandex_cis(
    rows: Sequence[dict[str, object]],
    *,
    bootstrap_samples: int,
    bootstrap_random_state: int,
) -> dict[str, BootstrapMetricCI]:
    """Bootstrap CIs for SPANDEX metrics, resampling at bacterium level."""
    by_bacterium: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in rows:
        by_bacterium[str(row["bacteria"])].append(row)

    bacteria_ids = tuple(sorted(by_bacterium.keys()))
    n_bacteria = len(bacteria_ids)
    rng = np.random.default_rng(bootstrap_random_state)

    metric_samples: dict[str, list[float]] = {name: [] for name in SPANDEX_METRIC_NAMES}

    progress_interval = max(1, bootstrap_samples // 5)
    for i in range(bootstrap_samples):
        if i == 0 or (i + 1) % progress_interval == 0 or i + 1 == bootstrap_samples:
            LOGGER.info("Bootstrap progress: %d/%d", i + 1, bootstrap_samples)

        sampled_indices = rng.integers(0, n_bacteria, size=n_bacteria)
        sampled_rows: list[dict[str, object]] = []
        for idx in sampled_indices.tolist():
            sampled_rows.extend(by_bacterium[bacteria_ids[idx]])

        metrics = evaluate_holdout_rows(sampled_rows)
        for name in SPANDEX_METRIC_NAMES:
            val = metrics.get(name)
            if val is not None:
                metric_samples[name].append(val)

    actual_metrics = evaluate_holdout_rows(list(rows))

    def _ci(values: Sequence[float]) -> tuple[Optional[float], Optional[float], int]:
        if not values:
            return None, None, 0
        arr = np.asarray(values, dtype=float)
        low, high = np.quantile(arr, [0.025, 0.975])
        return safe_round(float(low)), safe_round(float(high)), len(values)

    results = {}
    for name in SPANDEX_METRIC_NAMES:
        ci_low, ci_high, used = _ci(metric_samples[name])
        results[name] = BootstrapMetricCI(
            point_estimate=actual_metrics.get(name),
            ci_low=ci_low,
            ci_high=ci_high,
            bootstrap_samples_requested=bootstrap_samples,
            bootstrap_samples_used=used,
        )
    return results


def run_sx01_eval(
    *,
    candidate_module: ModuleType,
    context: Any,
    device_type: str,
    output_dir: Path,
) -> None:
    start_time = datetime.now(timezone.utc)

    # Pre-flight gate.
    preflight = run_preflight(output_dir)
    if preflight["status"] == "fail":
        LOGGER.error(
            "PRE-FLIGHT FAILED: Spearman rho=%.4f < 0.1 — graded nDCG is cosmetic. Consider cancelling SX04.",
            preflight.get("spearman_rho_positives", 0),
        )
        # Continue anyway — mAP and k-fold CV are still valuable.

    # Load full frame and apply clean-label filter.
    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    LOGGER.info("Full frame: %d pairs, %d bacteria", len(full_frame), full_frame["bacteria"].nunique())

    ambiguous_pairs = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous_pairs)].copy()
    LOGGER.info(
        "Clean-label frame: %d pairs (%d excluded), %d bacteria",
        len(clean_frame),
        len(full_frame) - len(clean_frame),
        clean_frame["bacteria"].nunique(),
    )

    # Load MLC scores for graded evaluation.
    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}
    LOGGER.info("MLC lookup: %d pairs", len(mlc_lookup))

    # Run k-fold CV on clean data.
    run_kfold_evaluation(
        candidate_module=candidate_module,
        context=context,
        full_frame=clean_frame,
        mlc_lookup=mlc_lookup,
        device_type=device_type,
        output_dir=output_dir,
    )

    elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
    LOGGER.info("SX01 evaluation completed in %.1f seconds", elapsed)


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
    LOGGER.info("SX01 eval starting at %s", datetime.now(timezone.utc).isoformat())

    candidate_module = load_module_from_path("sx01_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)

    run_sx01_eval(
        candidate_module=candidate_module,
        context=context,
        device_type=args.device_type,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
