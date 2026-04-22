#!/usr/bin/env python3
"""Compute per-pair SHAP values for the CH05 unified model and persist as Parquet.

Replays CH05's two-axis fold loops (bacteria-axis cv_group folds + phage-axis ICTV
StratifiedKFold), trains LightGBM seeds per fold, and runs shap.TreeExplainer on the
max-concentration pair frame for each held-out fold. SHAP values are averaged across
seeds per fold, then concatenated across folds so every pair appears exactly once per
axis. Writes three Parquet tables per axis:

- `shap_values.parquet` — wide form, one row per pair, one column per RFE-retained
  feature (plus `pair_id`, `bacteria`, `phage`, `source`). Features absent from a fold's
  RFE selection are encoded as NaN for pairs in that fold.
- `shap_base_values.parquet` — one row per pair with `base_value` (fold-level mean
  training-set prediction). Waterfall decomposition is `base_value + sum(shap) =
  predicted_logit`.
- `feature_values.parquet` — the raw feature matrix for the same pair set, same columns.

This script runs once per refresh of the CH05 canonical — it is not part of CH05 itself
because SHAP compute is optional (only the explainability UI consumes it). Expected
wallclock: 45-60 min on a 2023-vintage MacBook with `--num-workers 3`.
"""

from __future__ import annotations

import argparse
import logging
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Optional

import numpy as np
import pandas as pd
import shap
from lightgbm import LGBMClassifier

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    load_module_from_path,
    temporary_module_attribute,
)
from lyzortx.pipeline.autoresearch.ch04_eval import (
    build_clean_row_training_frame,
)
from lyzortx.pipeline.autoresearch.ch04_parallel import (
    prepare_fold_design_matrices,
    select_rfe_features,
)
from lyzortx.pipeline.autoresearch.ch05_eval import (
    BASEL_LOG10_PFU_ML,
    assign_phage_folds,
    load_unified_phage_family_map,
    load_unified_row_frame,
    patch_context_with_extended_slots,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch05_shap")


def _fit_booster_for_fold(
    *,
    candidate_module: ModuleType,
    train_design: pd.DataFrame,
    rfe_features: list[str],
    rfe_categorical: list[str],
    device_type: str,
    seed: int,
) -> LGBMClassifier:
    """Replay CH04's single-seed fit. Returns the trained estimator (not discarded)."""
    y_train = train_design["label_row_binary"].astype(int).to_numpy(dtype=int)
    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[rfe_features],
        y_train,
        categorical_feature=rfe_categorical,
    )
    return estimator


def _shap_for_holdout(
    estimator: LGBMClassifier,
    holdout_design: pd.DataFrame,
    rfe_features: list[str],
) -> tuple[np.ndarray, float]:
    """Compute SHAP values + the model's base value (log-odds) for a fold's holdout.

    `shap.TreeExplainer(model).shap_values(X)` returns the contribution of each feature
    to the model's log-odds output for the positive class. The base value is the mean
    model output over the training set (LightGBM stores this as `init_score`-like
    quantity; shap exposes it via `explainer.expected_value`).

    For binary classifiers shap >=0.43 returns a single 2D array of shape (n_rows, n_features)
    for the positive class; older versions may return a list of two arrays. We coerce
    to the single-array form here so downstream code is version-robust.
    """
    explainer = shap.TreeExplainer(estimator)
    shap_values = explainer.shap_values(holdout_design[rfe_features])
    if isinstance(shap_values, list):
        # Binary classifier, older shap: [neg_class_shap, pos_class_shap]. Use positive.
        shap_values = shap_values[1]
    base_value = explainer.expected_value
    if isinstance(base_value, (list, np.ndarray)):
        base_value = float(np.asarray(base_value).flatten()[-1])
    return np.asarray(shap_values, dtype=np.float32), float(base_value)


def _aggregate_fold_shap_to_pair_max_conc(
    *,
    holdout_design: pd.DataFrame,
    shap_matrix: np.ndarray,
    rfe_features: list[str],
    pair_source: pd.Series,
) -> pd.DataFrame:
    """Convert per-row SHAP (same row count as holdout_design) into per-pair SHAP.

    Strategy: each row in holdout_design is a (pair, log_dilution, replicate) triple.
    For the UI waterfall we want one SHAP vector per pair at its max-concentration
    observation — the same semantics as `select_pair_max_concentration_rows`. So we
    attach the SHAP matrix as extra columns on holdout_design, aggregate rows by
    (pair_id, max log10_pfu_ml), and for ties (multiple replicates at the same max
    concentration) take the mean SHAP across replicates.
    """
    shap_df = pd.DataFrame(
        shap_matrix,
        columns=[f"shap__{feat}" for feat in rfe_features],
        index=holdout_design.index,
    )
    merged = pd.concat(
        [
            holdout_design[["pair_id", "bacteria", "phage", "log10_pfu_ml"]].reset_index(drop=True),
            shap_df.reset_index(drop=True),
        ],
        axis=1,
    )
    max_conc_per_pair = merged.groupby("pair_id")["log10_pfu_ml"].transform("max")
    at_max = merged[merged["log10_pfu_ml"] == max_conc_per_pair]
    shap_cols = [c for c in at_max.columns if c.startswith("shap__")]
    agg = at_max.groupby(["pair_id", "bacteria", "phage"], as_index=False)[shap_cols].mean()
    agg["source"] = agg["pair_id"].map(pair_source)
    return agg


def _aggregate_fold_feature_values_to_pair_max_conc(
    *,
    holdout_design: pd.DataFrame,
    rfe_features: list[str],
    pair_source: pd.Series,
) -> pd.DataFrame:
    """Same shape as SHAP aggregator, but carrying raw feature values (not contributions)."""
    value_df = holdout_design[rfe_features + ["pair_id", "bacteria", "phage", "log10_pfu_ml"]].copy()
    max_conc_per_pair = value_df.groupby("pair_id")["log10_pfu_ml"].transform("max")
    at_max = value_df[value_df["log10_pfu_ml"] == max_conc_per_pair]
    # Use .agg with `first` for feature values (not mean — raw features don't average meaningfully
    # across replicates for categorical slots; all replicates of a pair share identical features
    # anyway since features are pair-level, so `first` is stable).
    agg = at_max.groupby(["pair_id", "bacteria", "phage"], as_index=False).first()
    agg = agg.drop(columns=["log10_pfu_ml"])
    agg["source"] = agg["pair_id"].map(pair_source)
    return agg


def run_axis(
    *,
    axis_name: str,
    fold_assignments: dict[str, int],
    entity_col: str,
    clean_rows: pd.DataFrame,
    candidate_module: ModuleType,
    context: Any,
    device_type: str,
    pair_source: pd.Series,
    seeds: tuple[int, ...],
    max_folds: Optional[int],
) -> dict[str, pd.DataFrame]:
    """Run CH05's axis loop with SHAP capture. Returns three DataFrames per axis."""
    folds_to_run = N_FOLDS if max_folds is None else min(max_folds, N_FOLDS)
    axis_shap_frames: list[pd.DataFrame] = []
    axis_values_frames: list[pd.DataFrame] = []
    axis_base_values: list[dict[str, object]] = []

    for fold_id in range(folds_to_run):
        holdout_ids = {k for k, f in fold_assignments.items() if f == fold_id}
        train_ids = {k for k, f in fold_assignments.items() if f != fold_id}
        fold_train = clean_rows[clean_rows[entity_col].isin(train_ids)].copy()
        fold_holdout = clean_rows[clean_rows[entity_col].isin(holdout_ids)].copy()
        LOGGER.info(
            "=== %s fold %d/%d: %d train %s (%d rows), %d holdout %s (%d rows) ===",
            axis_name,
            fold_id,
            folds_to_run - 1,
            len(train_ids),
            entity_col,
            len(fold_train),
            len(holdout_ids),
            entity_col,
            len(fold_holdout),
        )
        train_design, holdout_design, feature_columns, categorical_columns = prepare_fold_design_matrices(
            candidate_module=candidate_module,
            context=context,
            training_frame=fold_train,
            holdout_frame=fold_holdout,
        )
        rfe_features, rfe_categorical = select_rfe_features(
            train_design=train_design,
            feature_columns=feature_columns,
            categorical_columns=categorical_columns,
        )
        LOGGER.info("%s fold %d: %d RFE features", axis_name, fold_id, len(rfe_features))

        seed_shap_frames: list[np.ndarray] = []
        seed_base_values: list[float] = []
        for seed in seeds:
            estimator = _fit_booster_for_fold(
                candidate_module=candidate_module,
                train_design=train_design,
                rfe_features=rfe_features,
                rfe_categorical=rfe_categorical,
                device_type=device_type,
                seed=seed,
            )
            shap_matrix, base_value = _shap_for_holdout(estimator, holdout_design, rfe_features)
            seed_shap_frames.append(shap_matrix)
            seed_base_values.append(base_value)
            LOGGER.info(
                "%s fold %d seed %d: SHAP matrix shape=%s, base_value=%.4f",
                axis_name,
                fold_id,
                seed,
                shap_matrix.shape,
                base_value,
            )
        mean_shap = np.mean(np.stack(seed_shap_frames, axis=0), axis=0).astype(np.float32)
        mean_base = float(np.mean(seed_base_values))

        fold_shap = _aggregate_fold_shap_to_pair_max_conc(
            holdout_design=holdout_design,
            shap_matrix=mean_shap,
            rfe_features=rfe_features,
            pair_source=pair_source,
        )
        fold_shap["fold_id"] = fold_id
        axis_shap_frames.append(fold_shap)

        fold_values = _aggregate_fold_feature_values_to_pair_max_conc(
            holdout_design=holdout_design,
            rfe_features=rfe_features,
            pair_source=pair_source,
        )
        fold_values["fold_id"] = fold_id
        axis_values_frames.append(fold_values)

        for pair_id in fold_shap["pair_id"]:
            axis_base_values.append({"pair_id": pair_id, "base_value": mean_base, "fold_id": fold_id})

    shap_df = pd.concat(axis_shap_frames, ignore_index=True)
    values_df = pd.concat(axis_values_frames, ignore_index=True)
    base_df = pd.DataFrame(axis_base_values)
    LOGGER.info(
        "%s axis totals: %d pairs with SHAP, %d features (union across folds)",
        axis_name,
        len(shap_df),
        len([c for c in shap_df.columns if c.startswith("shap__")]),
    )
    return {"shap": shap_df, "values": values_df, "base": base_df}


def write_axis_outputs(axis_name: str, frames: dict[str, pd.DataFrame], output_dir: Path) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    prefix = {"bacteria_axis": "bacteria_axis", "phage_axis": "phage_axis"}[axis_name]
    for kind, df in frames.items():
        suffix = {"shap": "shap_values", "values": "feature_values", "base": "shap_base_values"}[kind]
        path = output_dir / f"{prefix}_{suffix}.parquet"
        df.to_parquet(path, index=False, compression="snappy")
        LOGGER.info("Wrote %s (%d rows, %.1f KB)", path, len(df), path.stat().st_size / 1024)


def run_shap_snapshot(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
    seeds: tuple[int, ...] = SEEDS,
    max_folds: Optional[int] = None,
    basel_log10_pfu_ml: float = BASEL_LOG10_PFU_ML,
) -> dict[str, dict[str, pd.DataFrame]]:
    output_dir.mkdir(parents=True, exist_ok=True)
    start = datetime.now(timezone.utc)
    LOGGER.info(
        "derive_shap_snapshot starting at %s (out=%s, seeds=%s, max_folds=%s)",
        start.isoformat(),
        output_dir,
        seeds,
        max_folds if max_folds is not None else "all",
    )

    unified = load_unified_row_frame(basel_log10_pfu_ml=basel_log10_pfu_ml)
    clean_rows = build_clean_row_training_frame(unified)
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]
    phage_family = load_unified_phage_family_map()

    candidate_module = load_module_from_path("ch05_shap_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context)

    LOGGER.info("Starting bacteria-axis pass")
    bact_mapping = bacteria_to_cv_group_map(clean_rows)
    bact_folds = assign_bacteria_folds(bact_mapping)
    bact_frames = run_axis(
        axis_name="bacteria_axis",
        fold_assignments=bact_folds,
        entity_col="bacteria",
        clean_rows=clean_rows,
        candidate_module=candidate_module,
        context=context,
        device_type=device_type,
        pair_source=pair_source,
        seeds=seeds,
        max_folds=max_folds,
    )
    write_axis_outputs("bacteria_axis", bact_frames, output_dir)

    LOGGER.info("Starting phage-axis pass")
    phages = sorted(clean_rows["phage"].unique())
    phage_folds = assign_phage_folds(phages, phage_family)
    phage_frames = run_axis(
        axis_name="phage_axis",
        fold_assignments=phage_folds,
        entity_col="phage",
        clean_rows=clean_rows,
        candidate_module=candidate_module,
        context=context,
        device_type=device_type,
        pair_source=pair_source,
        seeds=seeds,
        max_folds=max_folds,
    )
    write_axis_outputs("phage_axis", phage_frames, output_dir)

    elapsed = (datetime.now(timezone.utc) - start).total_seconds()
    LOGGER.info("derive_shap_snapshot complete in %.1f s", elapsed)
    return {"bacteria_axis": bact_frames, "phage_axis": phage_frames}


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", default="cpu")
    parser.add_argument("--out-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument(
        "--max-folds",
        type=int,
        default=None,
        help="Limit number of folds per axis (for smoke testing). Defaults to all 10.",
    )
    parser.add_argument(
        "--seeds",
        type=int,
        nargs="+",
        default=None,
        help="Override seed list (default: SEEDS from sx01_eval.py, usually 3 seeds).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    seeds = tuple(args.seeds) if args.seeds is not None else SEEDS
    run_shap_snapshot(
        device_type=args.device_type,
        output_dir=args.out_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
        seeds=seeds,
        max_folds=args.max_folds,
    )


if __name__ == "__main__":
    main()
