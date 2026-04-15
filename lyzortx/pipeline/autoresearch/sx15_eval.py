#!/usr/bin/env python3
"""SX15: Unified Guelin+BASEL k-fold evaluation framework (bacteria + phage axes).

Re-evaluates the wave-2 consolidated baseline (= SX10 unchanged since all wave-2 arms failed
their gates) under a unified Guelin+BASEL panel with stratified reporting per SX14's framework.

Two axes:
  - BACTERIA-AXIS: 10-fold StratifiedKFold on ~390 unique bacteria (369 Guelin + 25 BASEL minus
    overlap), stratified on (source, phylogroup). Every fold must contain >=2 BASEL and >=20 Guelin.
  - PHAGE-AXIS: 10-fold StratifiedKFold on 148 unique phages, stratified on phage family. Tests
    deployability on genuinely unseen phages.

Label unification:
  - Guelin: MLC 0-3 as-is.
  - BASEL: positive -> MLC=2 (default, Option B), negative -> MLC=0.
  - Sensitivity: rerun under Option A (BASEL+ -> MLC=1) and Option C (BASEL+ -> MLC=3).

Stratified evaluation (mandatory per SX14): every bootstrap reports aggregate + per-stratum
metrics across SX14 strata (within_family, cross_family, narrow_host_phage, phylogroup_orphan)
plus SX15-specific strata (source_basel, source_guelin; holdout_phage_source on phage axis).

Phage-axis limitation: held-out phages have zero training pairs, so per-phage blending cannot
apply. The phage-axis evaluation measures the all-pairs model only, not full SX10.

Usage:
    python -m lyzortx.pipeline.autoresearch.sx15_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import logging
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any

import numpy as np
import pandas as pd
from sklearn.model_selection import StratifiedKFold

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
from lyzortx.pipeline.autoresearch.sx03_eval import load_basel_interactions, patch_context_with_extended_slots
from lyzortx.pipeline.autoresearch.sx14_eval import (
    NARROW_HOST_LYSIS_THRESHOLD,
    PHYLOGROUP_ORPHAN_MAX_SIBLINGS,
    WITHIN_FAMILY_MIN_POSITIVES,
    bootstrap_per_stratum,
    compute_phage_lysis_rates,
    load_host_cv_group_map,
    load_host_phylogroup_map,
    load_interaction_matrix,
    load_phage_family_map,
    pivot_for_notebook,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    N_FOLDS,
    SEEDS,
    load_mlc_scores,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx15_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")

# Label sensitivity options for BASEL positives.
BASEL_POSITIVE_OPTIONS = {"A": 1, "B": 2, "C": 3}
BASEL_NEGATIVE_LABEL = 0
DEFAULT_OPTION = "B"  # BASEL+ -> MLC=2 (middle of positive range)
BASEL_PHAGE_METADATA = Path(".scratch/basel/S2_Table_phage_genome_data.xlsx")


def load_unified_phage_family_map() -> dict[str, str]:
    """Combined Guelin (from guelin_collection.csv) + BASEL (from S2 Table) phage→family mapping."""
    guelin = load_phage_family_map()
    basel_df = pd.read_excel(BASEL_PHAGE_METADATA, header=6)
    basel_df = basel_df[["Bas##", "ICTV family"]].dropna()
    basel = {str(r["Bas##"]): str(r["ICTV family"]) for _, r in basel_df.iterrows()}
    combined = {**guelin, **basel}
    LOGGER.info("Combined phage family map: %d Guelin + %d BASEL = %d phages", len(guelin), len(basel), len(combined))
    return combined


@dataclass
class AxisResult:
    """Container for one axis × option evaluation result."""

    axis: str  # "bacteria" or "phage"
    option: str  # "A", "B", or "C"
    predictions: pd.DataFrame  # per-pair with stratum labels
    metrics: pd.DataFrame  # stratified metric table


def load_basel_as_training_rows(basel_positive_mlc: int) -> pd.DataFrame:
    """Load BASEL pairs as SX01-compatible training rows with MLC unified per option.

    Returns columns: pair_id, bacteria, phage, label_any_lysis, mlc_score, source,
    training_weight_v3 (fixed to 1.0 per SX03 convention).
    """
    basel = load_basel_interactions()
    basel["label_any_lysis"] = (basel["interaction"] == 1).astype(int)
    basel["mlc_score"] = basel["label_any_lysis"].map({0: BASEL_NEGATIVE_LABEL, 1: basel_positive_mlc})
    basel["training_weight_v3"] = 1.0
    # Keep only confirmed observed cells (drop rows where interaction is missing).
    basel = basel.dropna(subset=["interaction"]).copy()
    return basel[["pair_id", "bacteria", "phage", "label_any_lysis", "mlc_score", "source", "training_weight_v3"]]


def build_unified_frame(basel_positive_mlc: int) -> pd.DataFrame:
    """Merge Guelin (MLC 0-3) + BASEL (mapped) into one training-ready DataFrame."""
    holdout = load_st03_holdout_frame()
    training = build_st03_training_frame()
    guelin = pd.concat([training, holdout], ignore_index=True)
    ambiguous = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    guelin = guelin[~guelin["pair_id"].isin(ambiguous)].copy()
    guelin["source"] = "guelin"

    # Attach MLC scores to Guelin pairs.
    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): float(r["mlc_score"]) for _, r in mlc_df.iterrows()}
    guelin["mlc_score"] = guelin.apply(lambda r: mlc_lookup.get((r["bacteria"], r["phage"]), 0.0), axis=1)

    basel = load_basel_as_training_rows(basel_positive_mlc)

    # Concat (BASEL shared bacteria contribute duplicate rows; we keep both since pair_id is unique).
    unified = pd.concat([guelin, basel], ignore_index=True)
    LOGGER.info(
        "Unified frame (BASEL+ -> MLC=%d): %d Guelin + %d BASEL = %d total pairs, %d unique bacteria, %d unique phages",
        basel_positive_mlc,
        len(guelin),
        len(basel),
        len(unified),
        unified["bacteria"].nunique(),
        unified["phage"].nunique(),
    )
    return unified


def stratified_bacteria_folds(
    unified: pd.DataFrame,
    host_phylogroup: dict[str, str],
    n_splits: int = N_FOLDS,
    random_state: int = 42,
) -> dict[str, int]:
    """StratifiedKFold on unique bacteria with (source, phylogroup) composite key."""
    bact_df = unified[["bacteria", "source"]].drop_duplicates().reset_index(drop=True)
    # If a bacterium appears with multiple sources (shared ECOR), pick a deterministic priority
    bact_df["source"] = bact_df.groupby("bacteria")["source"].transform(lambda g: "both" if len(g) > 1 else g.iloc[0])
    bact_df = bact_df.drop_duplicates(subset=["bacteria"]).reset_index(drop=True)
    bact_df["phylogroup"] = bact_df["bacteria"].map(host_phylogroup).fillna("UNKNOWN")
    bact_df["strat_key"] = bact_df["source"] + "__" + bact_df["phylogroup"]

    # Collapse rare strat keys (n<2) into a single bucket so StratifiedKFold can handle them.
    counts = bact_df["strat_key"].value_counts()
    bact_df.loc[bact_df["strat_key"].map(counts) < n_splits, "strat_key"] = "other"

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    fold_map: dict[str, int] = {}
    for fold_id, (_, test_idx) in enumerate(skf.split(bact_df["bacteria"].values, bact_df["strat_key"].values)):
        for b in bact_df.iloc[test_idx]["bacteria"].values:
            fold_map[str(b)] = fold_id

    # Sanity report
    per_fold = {i: [] for i in range(n_splits)}
    for b, f in fold_map.items():
        per_fold[f].append(b)
    for f, bacts in per_fold.items():
        n_basel = sum(1 for b in bacts if (unified[unified["bacteria"] == b]["source"] == "basel").any())
        n_guelin = len(bacts) - n_basel
        LOGGER.info("bacteria-axis fold %d: %d bacteria (%d guelin, %d basel)", f, len(bacts), n_guelin, n_basel)
    return fold_map


def stratified_phage_folds(
    unified: pd.DataFrame,
    phage_family: dict[str, str],
    n_splits: int = N_FOLDS,
    random_state: int = 42,
) -> dict[str, int]:
    """StratifiedKFold on unique phages with family stratification."""
    phages = sorted(unified["phage"].unique())
    phage_df = pd.DataFrame({"phage": phages})
    phage_df["family"] = phage_df["phage"].map(phage_family).fillna("UNKNOWN")

    # Collapse rare families
    counts = phage_df["family"].value_counts()
    phage_df.loc[phage_df["family"].map(counts) < n_splits, "family"] = "other"

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    fold_map: dict[str, int] = {}
    for fold_id, (_, test_idx) in enumerate(skf.split(phage_df["phage"].values, phage_df["family"].values)):
        for p in phage_df.iloc[test_idx]["phage"].values:
            fold_map[str(p)] = fold_id

    per_fold = {i: 0 for i in range(n_splits)}
    for _, f in fold_map.items():
        per_fold[f] += 1
    LOGGER.info("phage-axis fold sizes: %s", dict(per_fold))
    return fold_map


def train_and_predict_all_pairs_only(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    seed: int,
    device_type: str,
) -> list[dict[str, object]]:
    """All-pairs model prediction only (no per-phage blending). Used on phage-axis splits where
    held-out phages have zero training data.

    Mirrors the first half of train_and_predict_fold but skips the per-phage blending step.
    """
    from lyzortx.pipeline.autoresearch.candidate_replay import temporary_module_attribute

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

    y_train = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(train_design[rfe_features], y_train, sample_weight=sample_weight, categorical_feature=rfe_categorical)
    preds = estimator.predict_proba(holdout_design[rfe_features])[:, 1]

    rows = []
    for row, p in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        preds,
    ):
        rows.append(
            {
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(p)),
            }
        )
    return rows


def run_bacteria_axis(
    *,
    candidate_module: ModuleType,
    context: Any,
    unified: pd.DataFrame,
    host_phylogroup: dict[str, str],
    device_type: str,
) -> pd.DataFrame:
    """Bacteria-axis 10-fold CV with per-phage blending (full SX10 config)."""
    from lyzortx.pipeline.autoresearch.sx01_eval import train_and_predict_fold

    fold_map = stratified_bacteria_folds(unified, host_phylogroup)
    all_predictions: list[dict[str, object]] = []

    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_map.items() if f == fold_id}
        train_bacteria = set(fold_map.keys()) - holdout_bacteria
        holdout_frame = unified[unified["bacteria"].isin(holdout_bacteria)].copy()
        training_frame = unified[unified["bacteria"].isin(train_bacteria)].copy()
        LOGGER.info(
            "bacteria-axis fold %d: train=%d pairs, holdout=%d pairs (%d holdout bacteria)",
            fold_id,
            len(training_frame),
            len(holdout_frame),
            len(holdout_bacteria),
        )

        fold_rows: list[dict[str, object]] = []
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
            fold_rows.extend(rows)

        df = pd.DataFrame(fold_rows)
        agg = df.groupby(["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
            "predicted_probability"
        ].mean()
        all_predictions.extend(agg.to_dict(orient="records"))

    return pd.DataFrame(all_predictions)


def run_phage_axis(
    *,
    candidate_module: ModuleType,
    context: Any,
    unified: pd.DataFrame,
    phage_family: dict[str, str],
    device_type: str,
) -> pd.DataFrame:
    """Phage-axis 10-fold CV — all-pairs model only (no per-phage blending on held-out phages)."""
    fold_map = stratified_phage_folds(unified, phage_family)
    all_predictions: list[dict[str, object]] = []

    for fold_id in range(N_FOLDS):
        holdout_phages = {p for p, f in fold_map.items() if f == fold_id}
        train_phages = set(fold_map.keys()) - holdout_phages
        holdout_frame = unified[unified["phage"].isin(holdout_phages)].copy()
        training_frame = unified[unified["phage"].isin(train_phages)].copy()
        LOGGER.info(
            "phage-axis fold %d: train=%d pairs, holdout=%d pairs (%d holdout phages)",
            fold_id,
            len(training_frame),
            len(holdout_frame),
            len(holdout_phages),
        )

        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_and_predict_all_pairs_only(
                candidate_module=candidate_module,
                context=context,
                training_frame=training_frame,
                holdout_frame=holdout_frame,
                seed=seed,
                device_type=device_type,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_rows.extend(rows)

        df = pd.DataFrame(fold_rows)
        agg = df.groupby(["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
            "predicted_probability"
        ].mean()
        all_predictions.extend(agg.to_dict(orient="records"))

    return pd.DataFrame(all_predictions)


def attach_sx15_stratum_labels(
    predictions: pd.DataFrame,
    *,
    unified: pd.DataFrame,
    phage_family: dict[str, str],
    host_phylogroup: dict[str, str],
    host_cv_group: dict[str, int],
    phage_lysis: dict[str, float],
    axis: str,
) -> pd.DataFrame:
    """Attach SX14 strata + SX15-specific source strata.

    Includes: within_family, cross_family, narrow_host_phage, phylogroup_orphan, source_basel,
    source_guelin. On phage-axis, also adds holdout_phage_source_basel /
    holdout_phage_source_guelin.
    """
    out = predictions.copy()
    out["bacteria"] = out["bacteria"].astype(str)
    out["phage"] = out["phage"].astype(str)

    # Source: attach from unified frame (a pair_id uniquely maps to source).
    pair_source = unified.drop_duplicates(subset=["pair_id"])[["pair_id", "source"]].set_index("pair_id")["source"]
    out["source"] = out["pair_id"].map(pair_source).fillna("unknown")

    # mlc_score: take from unified (reflects the option's BASEL mapping for BASEL pairs)
    pair_mlc = unified.drop_duplicates(subset=["pair_id"])[["pair_id", "mlc_score"]].set_index("pair_id")["mlc_score"]
    out["mlc_score"] = out["pair_id"].map(pair_mlc).fillna(0.0)

    out["phage_family"] = out["phage"].map(phage_family).fillna("UNKNOWN")
    out["host_phylogroup"] = out["bacteria"].map(host_phylogroup).fillna("UNKNOWN")
    out["host_cv_group"] = out["bacteria"].map(host_cv_group).fillna(-1).astype(int)
    out["phage_lysis_rate"] = out["phage"].map(phage_lysis).fillna(0.0)

    out["stratum_narrow_host_phage"] = out["phage_lysis_rate"] < NARROW_HOST_LYSIS_THRESHOLD
    out["stratum_source_basel"] = out["source"] == "basel"
    out["stratum_source_guelin"] = out["source"] == "guelin"

    # Training-positive counts per (family, cv_group) per fold (using the unified frame).
    uni = unified.copy()
    uni["bacteria"] = uni["bacteria"].astype(str)
    uni["phage"] = uni["phage"].astype(str)
    uni["any_lysis"] = (uni["mlc_score"] >= 1).astype(int)
    uni["phage_family"] = uni["phage"].map(phage_family).fillna("UNKNOWN")
    uni["host_cv_group"] = uni["bacteria"].map(host_cv_group).fillna(-1).astype(int)

    fold_to_holdout: dict[int, set[str]]
    if axis == "bacteria":
        fold_to_holdout = out.groupby("fold_id")["bacteria"].agg(lambda x: set(x.astype(str))).to_dict()
        for fold_id, holdout_bact in fold_to_holdout.items():
            train_mask = ~uni["bacteria"].isin(holdout_bact) & (uni["any_lysis"] == 1)
            train_long = uni.loc[train_mask]
            counts = train_long.groupby(["phage_family", "host_cv_group"]).size().to_dict()
            fold_rows_mask = out["fold_id"] == fold_id

            def _classify(row, c=counts):
                key = (row["phage_family"], int(row["host_cv_group"]))
                n = c.get(key, 0)
                return pd.Series([n >= WITHIN_FAMILY_MIN_POSITIVES, n == 0])

            fold_part = out.loc[fold_rows_mask].apply(_classify, axis=1)
            fold_part.columns = ["stratum_within_family", "stratum_cross_family"]
            out.loc[fold_rows_mask, "stratum_within_family"] = fold_part["stratum_within_family"]
            out.loc[fold_rows_mask, "stratum_cross_family"] = fold_part["stratum_cross_family"]
    else:
        # phage-axis: training set is pairs whose PHAGE is not in this fold's holdout.
        fold_to_holdout = out.groupby("fold_id")["phage"].agg(lambda x: set(x.astype(str))).to_dict()
        for fold_id, holdout_phages in fold_to_holdout.items():
            train_mask = ~uni["phage"].isin(holdout_phages) & (uni["any_lysis"] == 1)
            train_long = uni.loc[train_mask]
            counts = train_long.groupby(["phage_family", "host_cv_group"]).size().to_dict()
            fold_rows_mask = out["fold_id"] == fold_id

            def _classify_p(row, c=counts):
                key = (row["phage_family"], int(row["host_cv_group"]))
                n = c.get(key, 0)
                return pd.Series([n >= WITHIN_FAMILY_MIN_POSITIVES, n == 0])

            fold_part = out.loc[fold_rows_mask].apply(_classify_p, axis=1)
            fold_part.columns = ["stratum_within_family", "stratum_cross_family"]
            out.loc[fold_rows_mask, "stratum_within_family"] = fold_part["stratum_within_family"]
            out.loc[fold_rows_mask, "stratum_cross_family"] = fold_part["stratum_cross_family"]

    out["stratum_within_family"] = out["stratum_within_family"].fillna(False).astype(bool)
    out["stratum_cross_family"] = out["stratum_cross_family"].fillna(False).astype(bool)

    # Phylogroup-orphan: count training phylogroup siblings per fold.
    bact_pgrp = {b: host_phylogroup.get(b, "UNKNOWN") for b in out["bacteria"].unique()}
    if axis == "bacteria":
        for fold_id, holdout_bact in fold_to_holdout.items():
            all_bact = set(unified["bacteria"].astype(str).unique())
            train_bact = all_bact - holdout_bact
            pgrp_counts: dict[str, int] = {}
            for b in train_bact:
                pgrp = bact_pgrp.get(b, "UNKNOWN")
                pgrp_counts[pgrp] = pgrp_counts.get(pgrp, 0) + 1
            fmask = out["fold_id"] == fold_id
            out.loc[fmask, "__train_pgrp_siblings"] = (
                out.loc[fmask, "host_phylogroup"].map(pgrp_counts).fillna(0).astype(int)
            )
    else:
        # phage-axis: all bacteria are in training for every fold, so pgrp-orphan is constant.
        all_bact = set(unified["bacteria"].astype(str).unique())
        pgrp_counts = {}
        for b in all_bact:
            pgrp = bact_pgrp.get(b, "UNKNOWN")
            pgrp_counts[pgrp] = pgrp_counts.get(pgrp, 0) + 1
        out["__train_pgrp_siblings"] = out["host_phylogroup"].map(pgrp_counts).fillna(0).astype(int)
    out["stratum_phylogroup_orphan"] = out["__train_pgrp_siblings"] <= PHYLOGROUP_ORPHAN_MAX_SIBLINGS
    out = out.drop(columns=["__train_pgrp_siblings"])

    # Phage-axis only: split holdout phages by source.
    if axis == "phage":
        phage_source = unified.drop_duplicates(subset=["phage"])[["phage", "source"]].set_index("phage")["source"]
        out["holdout_phage_source"] = out["phage"].map(phage_source).fillna("unknown")
        out["stratum_holdout_phage_guelin"] = out["holdout_phage_source"] == "guelin"
        out["stratum_holdout_phage_basel"] = out["holdout_phage_source"] == "basel"

    return out


def run_axis_with_sensitivity(
    *,
    axis: str,
    candidate_module: ModuleType,
    context: Any,
    host_phylogroup: dict[str, str],
    phage_family: dict[str, str],
    host_cv_group: dict[str, int],
    phage_lysis: dict[str, float],
    device_type: str,
    output_dir: Path,
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Run one axis across all 3 sensitivity options. Returns (all_predictions_df, all_metrics_df)."""
    all_preds: list[pd.DataFrame] = []
    all_metrics: list[dict[str, object]] = []

    stratum_cols = [
        "stratum_within_family",
        "stratum_cross_family",
        "stratum_narrow_host_phage",
        "stratum_phylogroup_orphan",
        "stratum_source_basel",
        "stratum_source_guelin",
    ]
    if axis == "phage":
        stratum_cols = stratum_cols + ["stratum_holdout_phage_guelin", "stratum_holdout_phage_basel"]

    for option_label, basel_plus_mlc in BASEL_POSITIVE_OPTIONS.items():
        LOGGER.info("=== SX15 %s-axis  Option %s  (BASEL+ -> MLC=%d) ===", axis.upper(), option_label, basel_plus_mlc)
        unified = build_unified_frame(basel_plus_mlc)

        if axis == "bacteria":
            preds = run_bacteria_axis(
                candidate_module=candidate_module,
                context=context,
                unified=unified,
                host_phylogroup=host_phylogroup,
                device_type=device_type,
            )
        else:
            preds = run_phage_axis(
                candidate_module=candidate_module,
                context=context,
                unified=unified,
                phage_family=phage_family,
                device_type=device_type,
            )

        enriched = attach_sx15_stratum_labels(
            preds,
            unified=unified,
            phage_family=phage_family,
            host_phylogroup=host_phylogroup,
            host_cv_group=host_cv_group,
            phage_lysis=phage_lysis,
            axis=axis,
        )
        enriched.insert(0, "option", option_label)
        enriched.insert(0, "axis", axis)
        all_preds.append(enriched)

        # Bootstrap per stratum
        rows = bootstrap_per_stratum(enriched, arm_id=f"sx10_{axis}_axis_opt{option_label}", stratum_cols=stratum_cols)
        for r in rows:
            r["axis"] = axis
            r["option"] = option_label
        all_metrics.extend(rows)

    preds_df = pd.concat(all_preds, ignore_index=True)
    metrics_df = pd.DataFrame(all_metrics)
    return preds_df, metrics_df


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    p.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    p.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    p.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    p.add_argument("--axes", type=str, default="bacteria,phage", help="Comma-separated: bacteria,phage")
    p.add_argument("--options", type=str, default="A,B,C", help="Sensitivity options to run")
    return p.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    start = datetime.now(timezone.utc)
    LOGGER.info("SX15 eval starting at %s", start.isoformat())

    global BASEL_POSITIVE_OPTIONS  # noqa: PLW0603
    # Filter sensitivity options
    requested_opts = [o.strip() for o in args.options.split(",") if o.strip()]
    for o in requested_opts:
        if o not in BASEL_POSITIVE_OPTIONS:
            raise ValueError(f"Unknown option {o}; valid: {list(BASEL_POSITIVE_OPTIONS)}")
    BASEL_POSITIVE_OPTIONS = {k: v for k, v in BASEL_POSITIVE_OPTIONS.items() if k in requested_opts}

    candidate_module = load_module_from_path("sx15_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)
    # Patch BASEL phage features (projection/stats/rbp_struct) into slot artifacts — required for unified eval.
    patch_context_with_extended_slots(context)

    phage_family = load_unified_phage_family_map()
    host_phylogroup = load_host_phylogroup_map()
    host_cv_group = load_host_cv_group_map()
    # Unified phage lysis rates: Guelin from interaction_matrix.csv + BASEL from BASEL interactions.
    long_matrix = load_interaction_matrix()
    phage_lysis = compute_phage_lysis_rates(long_matrix)
    basel = load_basel_interactions()
    basel_rates = basel.assign(any_lysis=(basel["interaction"] == 1).astype(int)).groupby("phage")["any_lysis"].mean()
    for p, r in basel_rates.items():
        phage_lysis[str(p)] = float(r)
    LOGGER.info(
        "Unified phage lysis rates: %d phages (min=%.3f, median=%.3f, max=%.3f)",
        len(phage_lysis),
        min(phage_lysis.values()),
        float(np.median(list(phage_lysis.values()))),
        max(phage_lysis.values()),
    )

    axes = [a.strip() for a in args.axes.split(",") if a.strip()]
    for axis in axes:
        if axis not in ("bacteria", "phage"):
            raise ValueError(f"Unknown axis {axis}")

    for axis in axes:
        preds_df, metrics_df = run_axis_with_sensitivity(
            axis=axis,
            candidate_module=candidate_module,
            context=context,
            host_phylogroup=host_phylogroup,
            phage_family=phage_family,
            host_cv_group=host_cv_group,
            phage_lysis=phage_lysis,
            device_type=args.device_type,
            output_dir=args.output_dir,
        )
        preds_path = args.output_dir / f"sx15_{axis}_axis_predictions.csv"
        metrics_path = args.output_dir / f"sx15_{axis}_axis_stratified_metrics.csv"
        preds_df.to_csv(preds_path, index=False)
        metrics_df.to_csv(metrics_path, index=False)
        LOGGER.info("Wrote %d predictions and %d metric rows for %s-axis", len(preds_df), len(metrics_df), axis)

    # Comparison table
    if all((args.output_dir / f"sx15_{a}_axis_stratified_metrics.csv").exists() for a in axes):
        all_metrics = pd.concat(
            [pd.read_csv(args.output_dir / f"sx15_{a}_axis_stratified_metrics.csv") for a in axes],
            ignore_index=True,
        )
        # Use default option B only for the headline comparison
        default_metrics = all_metrics[all_metrics["option"] == "B"].copy()
        default_metrics["arm_id"] = default_metrics["axis"] + "_axis (SX15 default BASEL+→MLC=2)"
        md_table = pivot_for_notebook(default_metrics)
        (args.output_dir / "sx15_comparison_table.md").write_text(md_table, encoding="utf-8")
        LOGGER.info("Wrote comparison table to %s/sx15_comparison_table.md", args.output_dir)

    elapsed = (datetime.now(timezone.utc) - start).total_seconds()
    LOGGER.info("SX15 complete in %.0f sec", elapsed)


if __name__ == "__main__":
    main()
