#!/usr/bin/env python3
"""CH03: Safety-net SX10 rerun on the row-expanded training frame.

Loads the row-expanded raw-observation frame (see `ch03_row_expansion`), collapses
rows back to pair level with any_lysis semantics, and runs the SX10 canonical
configuration (GT03 all_gates_rfe + AX02 per-phage blending, any_lysis labels,
SX10 feature bundle) under the CH02-fixed cv_group folds. Must reproduce the
CH02 revalidated numbers (AUC 0.8521, Brier 0.1317) within 0.005 — otherwise the
row-expansion plumbing has introduced a bug that would contaminate CH04's per-row
training.

Integrity guarantee: because all rows of a given bacterium are routed to the
same fold by the CH02 cv_group hash, no (pair, concentration, replicate) raw
observation can appear in both the train and test partition of a fold.

Usage:
    PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch03_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    load_module_from_path,
    safe_round,
)
from lyzortx.pipeline.autoresearch.ch03_row_expansion import (
    collapse_to_pair_level_for_training,
    load_row_expanded_frame,
    verify_rollup_matches_st02,
)
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.sx01_eval import (
    RAW_INTERACTIONS_PATH,
    bacteria_to_cv_group_map,
    load_mlc_scores,
    run_kfold_evaluation,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch03_row_expansion")

CH02_REVALIDATED_METRICS_PATH = Path("lyzortx/generated_outputs/ch02_cv_group_fix/ch02_sx10_revalidated_metrics.json")
REGRESSION_TOLERANCE = 0.005


def check_row_level_fold_integrity(row_frame: pd.DataFrame, pair_level_frame: pd.DataFrame) -> None:
    """Verify no (pair, log_dilution, replicate) row appears in both train and test of any fold.

    Partitions the row-expanded frame into per-fold train / test subsets using the CH02
    cv_group hash and asserts disjointness on the raw-observation key. This is the
    load-bearing CH03 invariant check: the acceptance criterion `no (pair, concentration,
    replicate) row appears in both train and test` must hold for every fold, otherwise
    CH04's per-row training would leak observations across the split.
    """
    from lyzortx.pipeline.autoresearch.sx01_eval import N_FOLDS, assign_bacteria_folds

    mapping = bacteria_to_cv_group_map(pair_level_frame)
    fold_assignments = assign_bacteria_folds(mapping)

    rows_with_fold = row_frame.assign(fold_id=row_frame["bacteria"].map(fold_assignments))
    if rows_with_fold["fold_id"].isna().any():
        unassigned = rows_with_fold[rows_with_fold["fold_id"].isna()]["bacteria"].unique()
        raise ValueError(f"{len(unassigned)} bacteria missing fold assignment: {unassigned[:5].tolist()}")

    observation_key_cols = ["pair_id", "log_dilution", "replicate"]
    row_keys = pd.MultiIndex.from_frame(rows_with_fold[observation_key_cols])
    if not row_keys.is_unique:
        duplicates = rows_with_fold[row_keys.duplicated(keep=False)]
        raise ValueError(
            f"raw frame has duplicate (pair_id, log_dilution, replicate) rows — "
            f"{len(duplicates)} rows with shared keys (example: {duplicates.head(3).to_dict(orient='records')})"
        )

    for fold_id in range(N_FOLDS):
        fold_mask = rows_with_fold["fold_id"] == fold_id
        train_keys = pd.MultiIndex.from_frame(rows_with_fold.loc[~fold_mask, observation_key_cols])
        test_keys = pd.MultiIndex.from_frame(rows_with_fold.loc[fold_mask, observation_key_cols])
        overlap = train_keys.intersection(test_keys)
        if len(overlap) > 0:
            raise ValueError(
                f"fold {fold_id}: {len(overlap)} (pair, log_dilution, replicate) rows "
                f"appear in both train and test partitions (example: {overlap[:3].tolist()})"
            )

    LOGGER.info(
        "Fold integrity verified: %d raw rows / %d pairs / %d bacteria route to %d folds, "
        "train and test disjoint on (pair_id, log_dilution, replicate) across all %d folds",
        len(row_frame),
        row_frame["pair_id"].nunique(),
        rows_with_fold["bacteria"].nunique(),
        rows_with_fold["fold_id"].nunique(),
        N_FOLDS,
    )


def load_ch02_baseline_metrics() -> dict[str, float]:
    if not CH02_REVALIDATED_METRICS_PATH.exists():
        raise FileNotFoundError(
            f"CH02 baseline metrics not found: {CH02_REVALIDATED_METRICS_PATH}. Run CH02 first (see track_CHISEL.md)."
        )
    with open(CH02_REVALIDATED_METRICS_PATH, encoding="utf-8") as f:
        data = json.load(f)
    return {
        "auc": data["chisel_fixed_folds"]["holdout_roc_auc"]["point_estimate"],
        "brier": data["chisel_fixed_folds"]["holdout_brier_score"]["point_estimate"],
    }


def run_ch03_eval(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.now(timezone.utc)

    LOGGER.info("CH03 evaluation starting at %s", start_time.isoformat())
    row_frame = load_row_expanded_frame()

    rollup_summary = verify_rollup_matches_st02(row_frame)
    LOGGER.info("Rollup-vs-ST02 check: %s", rollup_summary)
    if rollup_summary["n_mismatches"] != 0:
        raise AssertionError(
            f"ST01B rollup disagrees with ST02 on {rollup_summary['n_mismatches']} pairs — "
            "row-expansion does not preserve any_lysis semantics"
        )

    pair_frame = collapse_to_pair_level_for_training(row_frame)
    check_row_level_fold_integrity(row_frame, pair_frame)

    # Persist the row-expanded frame artifact (CH04 reuses this). CSV keeps the
    # environment dependency-free — pyarrow/fastparquet not required.
    expanded_path = output_dir / "ch03_expanded_training_frame.csv"
    row_frame.to_csv(expanded_path, index=False)
    LOGGER.info("Wrote row-expanded frame: %s (%d rows)", expanded_path, len(row_frame))

    # Mirror sx01_eval's clean-label + split filtering on the collapsed pair frame.
    pair_frame["label_any_lysis"] = pair_frame["label_hard_any_lysis"]
    training = pair_frame[
        (pair_frame["split_holdout"] == "train_non_holdout") & (pair_frame["is_hard_trainable"] == "1")
    ].copy()
    holdout = pair_frame[
        (pair_frame["split_holdout"] == "holdout_test") & (pair_frame["is_hard_trainable"] == "1")
    ].copy()
    full_frame = pd.concat([training, holdout], ignore_index=True)
    LOGGER.info(
        "CH03 pair-level frame: %d training, %d holdout, %d bacteria",
        len(training),
        len(holdout),
        full_frame["bacteria"].nunique(),
    )

    ambiguous_pairs = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous_pairs)].copy()
    LOGGER.info(
        "Clean-label CH03 frame: %d pairs (%d excluded), %d bacteria",
        len(clean_frame),
        len(full_frame) - len(clean_frame),
        clean_frame["bacteria"].nunique(),
    )

    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    candidate_module = load_module_from_path("ch03_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)

    kfold_output_dir = output_dir / "sx10_on_row_expanded"
    run_kfold_evaluation(
        candidate_module=candidate_module,
        context=context,
        full_frame=clean_frame,
        mlc_lookup=mlc_lookup,
        device_type=device_type,
        output_dir=kfold_output_dir,
    )

    with open(kfold_output_dir / "bootstrap_results.json", encoding="utf-8") as f:
        bootstrap = json.load(f)
    ch03_auc = float(bootstrap["holdout_roc_auc"]["point_estimate"])
    ch03_brier = float(bootstrap["holdout_brier_score"]["point_estimate"])

    baseline = load_ch02_baseline_metrics()
    delta_auc = ch03_auc - baseline["auc"]
    delta_brier = ch03_brier - baseline["brier"]
    passed = abs(delta_auc) < REGRESSION_TOLERANCE and abs(delta_brier) < REGRESSION_TOLERANCE

    regression = {
        "task_id": "CH03",
        "baseline_source": str(CH02_REVALIDATED_METRICS_PATH),
        "ch02_revalidated": baseline,
        "ch03_row_expanded": {"auc": safe_round(ch03_auc), "brier": safe_round(ch03_brier)},
        "delta": {"auc": safe_round(delta_auc), "brier": safe_round(delta_brier)},
        "tolerance": REGRESSION_TOLERANCE,
        "passed": passed,
        "rollup_vs_st02": rollup_summary,
        "n_raw_observations": int(len(row_frame)),
        "n_pairs": int(row_frame["pair_id"].nunique()),
        "elapsed_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 1),
    }
    regression_path = output_dir / "ch03_regression_check.json"
    with open(regression_path, "w", encoding="utf-8") as f:
        json.dump(regression, f, indent=2)
    LOGGER.info("Regression check: %s (wrote %s)", "PASS" if passed else "FAIL", regression_path)

    if not passed:
        raise AssertionError(
            f"Regression tolerance exceeded: delta_auc={delta_auc:.4f}, "
            f"delta_brier={delta_brier:.4f} (tolerance={REGRESSION_TOLERANCE})"
        )
    return regression


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
    run_ch03_eval(
        device_type=args.device_type,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
    )


if __name__ == "__main__":
    main()
