"""CH04/CH05 parallel training helpers: split fold prep from per-seed fit.

Refactors the CH04/CH05 training loop so the per-fold preparation (design matrix
build + RFE feature selection) runs once per fold, and the three seed fits can be
dispatched to a `multiprocessing.Pool` instead of running sequentially. This is
CH06's ENGINEERING PRE-FLIGHT — the 4-arm sweep inherits and multiplies whatever
runtime the baseline has, so we fix the loop first.

Determinism contract: sequential (num_workers=1) and parallel (num_workers>1) must
produce bit-identical per-row predictions under identical inputs. RFE uses
`seed=42` internally (hardcoded RFE_SEED; independent of the per-fold SEEDS from
`sx01_eval`), so RFE output is the same across all three seeds in a fold;
computing it once per fold is the other half of the speedup.

Workers are spawned via `multiprocessing.get_context("spawn")` so they do not
inherit the parent's loaded `candidate_module`. Each worker reloads the module
freshly from `candidate_dir`. Design matrices are pickled through the Pool's IPC
(~0.5 GB per matrix, negligible vs the ~60-90 s fit time).
"""

from __future__ import annotations

import logging
import multiprocessing as mp
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd

from lyzortx.pipeline.autoresearch.candidate_replay import (
    load_module_from_path,
    safe_round,
    temporary_module_attribute,
)
from lyzortx.pipeline.autoresearch.gt03_eval import apply_rfe

LOGGER = logging.getLogger(__name__)

# CH04 seed-independent RFE seed (matches the prior inline hardcoded 42).
RFE_SEED = 42

HOST_SLOTS = ("host_surface", "host_typing", "host_stats", "host_defense")
PHAGE_SLOTS = ("phage_projection", "phage_stats")
# Allowlist of design-matrix column prefixes that reach RFE + LightGBM. Any slot
# whose columns are not covered here is silently filtered out of the feature set
# (columns remain in the design matrix for downstream joins, but no model ever
# sees them). CH08 caught this: earlier wave-2 ablation tickets attached
# `host_omp_kmer` and `phage_moriniere_kmer` slots to the context but their
# features were silently dropped because the prefixes were missing here. When
# adding a new slot family whose columns carry a distinct prefix, register it
# below.
FEATURE_COLUMN_PREFIXES = (
    "host_surface__",
    "host_typing__",
    "host_stats__",
    "host_defense__",
    "host_omp_kmer__",
    "host_omp_cluster__",
    "phage_projection__",
    "phage_stats__",
    "phage_moriniere_kmer__",
    "pair_depo_capsule__",
    "pair_receptor_omp__",
    "pair_concentration__",
)
ROW_LEVEL_COLUMNS = (
    "pair_concentration__log10_pfu_ml",
    "label_row_binary",
    "log10_pfu_ml",
    "log_dilution",
    "replicate",
)


def prepare_fold_design_matrices(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    host_slots: tuple[str, ...] = HOST_SLOTS,
    phage_slots: tuple[str, ...] = PHAGE_SLOTS,
) -> tuple[pd.DataFrame, pd.DataFrame, list[str], list[str]]:
    """Build train and holdout design matrices for one fold. Deterministic, no seed.

    Returns `(train_design, holdout_design, feature_columns, categorical_columns)`.
    The returned matrices carry the row-level columns (`label_row_binary`,
    `log10_pfu_ml`, etc.) appended after the feature prefixes, so callers can
    select `feature_columns` for fitting and still read the row metadata for
    writing prediction rows.
    """
    from lyzortx.pipeline.autoresearch.derive_pairwise_depo_capsule_features import (
        compute_pairwise_depo_capsule_features,
    )
    from lyzortx.pipeline.autoresearch.derive_pairwise_receptor_omp_features import (
        compute_pairwise_receptor_omp_features,
    )

    host_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=list(host_slots), entity_key="bacteria"
    )
    phage_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=list(phage_slots), entity_key="phage"
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

    for column in ROW_LEVEL_COLUMNS:
        train_design[column] = training_frame[column].to_numpy()
        holdout_design[column] = holdout_frame[column].to_numpy()

    feature_columns = [col for col in train_design.columns if col.startswith(FEATURE_COLUMN_PREFIXES)]
    categorical_columns = [col for col in (host_categorical + phage_categorical) if col in feature_columns]
    return train_design, holdout_design, feature_columns, categorical_columns


def select_rfe_features(
    *,
    train_design: pd.DataFrame,
    feature_columns: list[str],
    categorical_columns: list[str],
) -> tuple[list[str], list[str]]:
    """RFE feature selection once per fold. Deterministic via hardcoded RFE_SEED=42.

    Returns `(rfe_features, rfe_categorical)`. Because RFE_SEED is independent of
    the per-seed SEEDS, RFE output is identical across all three seeds in a fold —
    this is why computing it once is a pure speedup, not an approximation.
    """
    y_train = train_design["label_row_binary"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train, seed=RFE_SEED)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    return rfe_features, rfe_categorical


def _fit_worker(args: tuple) -> tuple[int, list[dict[str, object]], pd.DataFrame]:
    """Pool worker: fit LightGBM for one seed on a prepared fold, return (seed, rows, feat_importance).

    Top-level module function so Pool can pickle it. Each worker re-loads
    `candidate_module` (Python module objects are not picklable) from
    `candidate_dir`. Design matrices arrive via Pool IPC; feature columns and
    categorical lists arrive alongside.
    """
    (
        seed,
        candidate_dir,
        train_design,
        holdout_design,
        rfe_features,
        rfe_categorical,
        device_type,
        arm_id,
    ) = args
    candidate_module = load_module_from_path("ch04_candidate", Path(candidate_dir) / "train.py")
    return _fit_one_seed(
        seed=seed,
        candidate_module=candidate_module,
        train_design=train_design,
        holdout_design=holdout_design,
        rfe_features=rfe_features,
        rfe_categorical=rfe_categorical,
        device_type=device_type,
        arm_id=arm_id,
    )


def _fit_one_seed(
    *,
    seed: int,
    candidate_module: ModuleType,
    train_design: pd.DataFrame,
    holdout_design: pd.DataFrame,
    rfe_features: list[str],
    rfe_categorical: list[str],
    device_type: str,
    arm_id: str,
) -> tuple[int, list[dict[str, object]], pd.DataFrame]:
    """Pure fit+predict for one seed. Shared between sequential path and worker path."""
    y_train = train_design["label_row_binary"].astype(int).to_numpy(dtype=int)
    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[rfe_features],
        y_train,
        categorical_feature=rfe_categorical,
    )
    predictions = estimator.predict_proba(holdout_design[rfe_features])[:, 1]

    rows: list[dict[str, object]] = []
    holdout_meta = holdout_design[
        ["pair_id", "bacteria", "phage", "label_row_binary", "log_dilution", "log10_pfu_ml", "replicate"]
    ]
    for (_, row), probability in zip(holdout_meta.iterrows(), predictions):
        rows.append(
            {
                "arm_id": arm_id,
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "log_dilution": int(row["log_dilution"]),
                "log10_pfu_ml": float(row["log10_pfu_ml"]),
                "replicate": int(row["replicate"]),
                "label_row_binary": int(row["label_row_binary"]),
                "predicted_probability": safe_round(float(probability)),
            }
        )

    from lyzortx.pipeline.autoresearch.ch04_eval import (
        CONCENTRATION_FEATURE_COLUMN,
    )  # avoid circular import at module load

    feature_importance = pd.DataFrame(
        {
            "feature": rfe_features,
            "importance": estimator.feature_importances_,
        }
    )
    feature_importance["is_concentration_feature"] = feature_importance["feature"] == CONCENTRATION_FEATURE_COLUMN
    return seed, rows, feature_importance


def fit_seeds(
    *,
    seeds: tuple[int, ...],
    candidate_module: ModuleType,
    candidate_dir: Path,
    train_design: pd.DataFrame,
    holdout_design: pd.DataFrame,
    rfe_features: list[str],
    rfe_categorical: list[str],
    device_type: str,
    num_workers: int,
    arm_id: str = "chisel_baseline",
) -> list[tuple[int, list[dict[str, object]], pd.DataFrame]]:
    """Fit LightGBM for each seed. Sequential when num_workers<=1, pooled otherwise.

    Result order matches `seeds` argument regardless of worker completion order
    (`Pool.map` preserves input order). Determinism is a contract: sequential
    and parallel paths must yield bit-identical per-row predictions.
    """
    if num_workers <= 1:
        LOGGER.info("fit_seeds: sequential (num_workers=%d, n_seeds=%d)", num_workers, len(seeds))
        return [
            _fit_one_seed(
                seed=seed,
                candidate_module=candidate_module,
                train_design=train_design,
                holdout_design=holdout_design,
                rfe_features=rfe_features,
                rfe_categorical=rfe_categorical,
                device_type=device_type,
                arm_id=arm_id,
            )
            for seed in seeds
        ]

    LOGGER.info("fit_seeds: pooled (num_workers=%d, n_seeds=%d)", num_workers, len(seeds))
    tasks = [
        (
            seed,
            str(candidate_dir),
            train_design,
            holdout_design,
            rfe_features,
            rfe_categorical,
            device_type,
            arm_id,
        )
        for seed in seeds
    ]
    ctx = mp.get_context("spawn")
    with ctx.Pool(processes=num_workers) as pool:
        return list(pool.map(_fit_worker, tasks))
