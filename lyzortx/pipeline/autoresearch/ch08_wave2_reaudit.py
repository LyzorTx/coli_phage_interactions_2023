#!/usr/bin/env python3
"""CH08: SX12 + SX13 wave-2 feature-family re-audit under the CHISEL frame.

Verifies that the SPANDEX wave-2 nulls (SX12 Moriniere 815-kmer phage-side
features; SX13 host OMP 5-mer features) still hold under the CHISEL training
frame (per-row binary labels, `pair_concentration__log10_pfu_ml` feature,
all-pairs only, neat-only filter).

SX11 (loss-function ablation) is explicitly out of scope per plan.yml — it
chased MLC-graded potency which no longer exists in the CHISEL label frame.

Arms evaluated (each is a full CH04 10-fold bacteria-axis CV on the
369-bacterium Guelin panel):
  1. baseline   — CH04 canonical (reused from lyzortx/generated_outputs/ch04_chisel_baseline/)
  2. sx12       — baseline + `phage_moriniere_kmer` slot (815 k-mers)
  3. sx13       — baseline + `host_omp_kmer` slot (marginal arm only)

Deltas are computed by paired bacterium-level bootstrap: the same 1000 resamples
(same random state, same unit IDs) are applied to both baseline and +slot
predictions so the AUC/Brier deltas are within-resample, preserving the pairing.

Artifacts under `lyzortx/generated_outputs/ch08_wave2_reaudit/`:
  - ch08_sx12_delta.json — SX12 arm: AUC/Brier with CIs, paired delta CI
  - ch08_sx13_delta.json — SX13 arm: same, marginal arm only
  - ch08_summary.csv     — baseline vs +sx12 vs +sx13 one-row-per-arm table
  - ch08_<arm>_predictions.csv — per-arm pair-level predictions (for audit)

Expected outcome: null (within ~0.5 pp of CH04 baseline) per the SX12/SX13
knowledge units `kmer-receptor-expansion-neutral` and
`host-omp-variation-unpredictive`. If CI on the delta is disjoint from zero,
the wave-2 null is reopened and the corresponding knowledge unit must be
reverted.

Usage:
    python -m lyzortx.pipeline.autoresearch.ch08_wave2_reaudit --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    load_module_from_path,
    safe_round,
)
from lyzortx.pipeline.autoresearch.ch03_row_expansion import load_row_expanded_frame
from lyzortx.pipeline.autoresearch.ch04_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    build_clean_row_training_frame,
    compute_aggregate_auc_brier,
    select_pair_max_concentration_rows,
)
from lyzortx.pipeline.autoresearch.ch04_parallel import (
    HOST_SLOTS as BASE_HOST_SLOTS,
    PHAGE_SLOTS as BASE_PHAGE_SLOTS,
    fit_seeds,
    prepare_fold_design_matrices,
    select_rfe_features,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)
from lyzortx.pipeline.autoresearch.sx12_eval import (
    DEFAULT_KMER_SLOT_PATH as DEFAULT_PHAGE_KMER_SLOT_PATH,
    MORINIERE_SLOT_NAME,
    attach_moriniere_slot,
)
from lyzortx.pipeline.autoresearch.sx13_eval import (
    HOST_OMP_KMER_SLOT,
    attach_csv_slot,
)
from lyzortx.pipeline.autoresearch.build_host_omp_kmer_slot import (
    DEFAULT_SLOT_OUTPUT_PATH as DEFAULT_HOST_OMP_KMER_SLOT_PATH,
    FEATURE_PREFIX as HOST_OMP_KMER_PREFIX,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch08_wave2_reaudit")
CH04_CANONICAL_DIR = Path("lyzortx/generated_outputs/ch04_chisel_baseline")
CH04_CANONICAL_PREDICTIONS = CH04_CANONICAL_DIR / "ch04_predictions.csv"
CH04_CANONICAL_METRICS = CH04_CANONICAL_DIR / "ch04_aggregate_metrics.json"

ARM_BASELINE = "baseline"
ARM_SX12 = "sx12"
ARM_SX13 = "sx13"

# Kmer slots have 815 (phage Moriniere) and 5546 (host OMP) binary features. With
# RFECV step=0.1, cv=5, and 281K-row training frames the full-feature RFE runs
# at >60min/fold (20 folds = 20h+) — infeasible per plan.yml's "estimate before
# running and confirm wallclock is feasible" criterion. We cap each kmer slot
# at the top-K features by binary variance (maximally discriminative across
# entities). K=100 chosen empirically to land RFE fold time in the CH04
# baseline range (~2-3 min/fold). Top-K by variance preserves:
#  * sparse receptor-class-specific signals (e.g. kmers in ~10% of phages —
#    high variance, kept)
#  * balanced OMP 5mers (support near 50% — maximum variance, kept)
# Features with near-constant presence (variance → 0) are dropped.
#
# Rank is computed entity-level once per slot before any fold split — leakage-
# free. Threshold K is recorded in each arm's delta JSON alongside the delta
# estimates.
KMER_TOPK_PER_SLOT = 100


def _aggregate_fold_rows_to_pairs(fold_rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(fold_rows)
    aggregated = (
        df.groupby(
            [
                "fold_id",
                "pair_id",
                "bacteria",
                "phage",
                "log_dilution",
                "log10_pfu_ml",
                "replicate",
                "label_row_binary",
            ],
            as_index=False,
        )["predicted_probability"]
        .mean()
        .sort_values(["bacteria", "phage", "log10_pfu_ml", "replicate"])
    )
    return aggregated


def run_ch04_variant(
    *,
    arm_id: str,
    phage_slots: tuple[str, ...],
    host_slots: tuple[str, ...],
    slot_attach_fns: list,
    device_type: str,
    cache_dir: Path,
    candidate_dir: Path,
    output_dir: Path,
    num_workers: int,
    drop_high_titer_only_positives: bool,
) -> pd.DataFrame:
    """Run one CH04 bacteria-axis 10-fold variant with an extended slot bundle.

    `slot_attach_fns` are callables taking the loaded context and mutating it in place
    (e.g., `attach_moriniere_slot`, `attach_csv_slot`). `phage_slots` / `host_slots`
    are the slot name tuples passed to `prepare_fold_design_matrices` so the new slot
    participates in feature table construction.
    """
    LOGGER.info(
        "=== CH08 variant '%s' starting: phage_slots=%s, host_slots=%s ===",
        arm_id,
        list(phage_slots),
        list(host_slots),
    )
    row_frame = load_row_expanded_frame()
    clean_rows = build_clean_row_training_frame(
        row_frame, drop_high_titer_only_positives=drop_high_titer_only_positives
    )
    training = clean_rows[
        (clean_rows["split_holdout"] == "train_non_holdout") & (clean_rows["is_hard_trainable"] == "1")
    ].copy()
    holdout = clean_rows[
        (clean_rows["split_holdout"] == "holdout_test") & (clean_rows["is_hard_trainable"] == "1")
    ].copy()
    full_frame = pd.concat([training, holdout], ignore_index=True)
    mapping = bacteria_to_cv_group_map(full_frame)
    fold_assignments = assign_bacteria_folds(mapping)

    candidate_module = load_module_from_path(f"ch08_{arm_id}_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    for fn in slot_attach_fns:
        fn(context)

    all_per_row: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_b = {b for b, f in fold_assignments.items() if f == fold_id}
        train_b = {b for b, f in fold_assignments.items() if f != fold_id}
        fold_train = full_frame[full_frame["bacteria"].isin(train_b)].copy()
        fold_holdout = full_frame[full_frame["bacteria"].isin(holdout_b)].copy()
        LOGGER.info(
            "=== CH08/%s Fold %d: %d train bacteria (%d rows), %d holdout bacteria (%d rows) ===",
            arm_id,
            fold_id,
            len(train_b),
            len(fold_train),
            len(holdout_b),
            len(fold_holdout),
        )
        train_design, holdout_design, feature_columns, categorical_columns = prepare_fold_design_matrices(
            candidate_module=candidate_module,
            context=context,
            training_frame=fold_train,
            holdout_frame=fold_holdout,
            host_slots=host_slots,
            phage_slots=phage_slots,
        )
        rfe_features, rfe_categorical = select_rfe_features(
            train_design=train_design,
            feature_columns=feature_columns,
            categorical_columns=categorical_columns,
        )
        seed_results = fit_seeds(
            seeds=SEEDS,
            candidate_module=candidate_module,
            candidate_dir=candidate_dir,
            train_design=train_design,
            holdout_design=holdout_design,
            rfe_features=rfe_features,
            rfe_categorical=rfe_categorical,
            device_type=device_type,
            num_workers=num_workers,
        )
        fold_seed_rows: list[dict[str, object]] = []
        for _seed, rows, _fi in seed_results:
            for r in rows:
                r["fold_id"] = fold_id
            fold_seed_rows.extend(rows)
        aggregated = _aggregate_fold_rows_to_pairs(fold_seed_rows)
        all_per_row.extend(aggregated.to_dict(orient="records"))
    per_row_df = pd.DataFrame(all_per_row)
    pair_pred = select_pair_max_concentration_rows(per_row_df)
    pair_pred.to_csv(output_dir / f"ch08_{arm_id}_predictions.csv", index=False)
    LOGGER.info("CH08 variant '%s' done: %d pair-level predictions", arm_id, len(pair_pred))
    return pair_pred


def _bootstrap_paired_by_bacterium(
    baseline_rows: Sequence[dict[str, object]],
    variant_rows: Sequence[dict[str, object]],
    *,
    bootstrap_samples: int = BOOTSTRAP_SAMPLES,
    bootstrap_random_state: int = BOOTSTRAP_RANDOM_STATE,
) -> dict[str, object]:
    """Paired bacterium-level bootstrap: baseline and variant use same resample indices.

    Returns point AUC/Brier for baseline and variant plus the delta distribution with CIs.
    Requires pair_id sets to match — we align on pair_id before resampling.
    """
    base_by_pair = {str(r["pair_id"]): r for r in baseline_rows}
    var_by_pair = {str(r["pair_id"]): r for r in variant_rows}
    common_pairs = sorted(set(base_by_pair.keys()) & set(var_by_pair.keys()))
    if len(common_pairs) < max(len(base_by_pair), len(var_by_pair)) * 0.95:
        LOGGER.warning(
            "Baseline and variant share only %d pairs (base=%d, variant=%d) — delta may be biased",
            len(common_pairs),
            len(base_by_pair),
            len(var_by_pair),
        )
    base_rows = [base_by_pair[p] for p in common_pairs]
    var_rows = [var_by_pair[p] for p in common_pairs]

    grouped_base: dict[str, list[int]] = defaultdict(list)
    for i, r in enumerate(base_rows):
        grouped_base[str(r["bacteria"])].append(i)
    unit_ids = tuple(sorted(grouped_base.keys()))
    rng = np.random.default_rng(bootstrap_random_state)

    base_aucs: list[float] = []
    base_briers: list[float] = []
    var_aucs: list[float] = []
    var_briers: list[float] = []
    delta_aucs: list[float] = []
    delta_briers: list[float] = []

    progress_interval = max(1, bootstrap_samples // 5)
    for i in range(bootstrap_samples):
        if i == 0 or (i + 1) % progress_interval == 0 or i + 1 == bootstrap_samples:
            LOGGER.info("Paired bootstrap progress: %d/%d", i + 1, bootstrap_samples)
        idx = rng.integers(0, len(unit_ids), size=len(unit_ids))
        row_indices: list[int] = []
        for j in idx.tolist():
            row_indices.extend(grouped_base[unit_ids[j]])
        if not row_indices:
            continue
        y_base = np.array([int(base_rows[k]["label_row_binary"]) for k in row_indices])
        y_var = np.array([int(var_rows[k]["label_row_binary"]) for k in row_indices])
        assert np.array_equal(y_base, y_var), "Baseline and variant label vectors diverged on pair-id alignment"
        p_base = np.array([float(base_rows[k]["predicted_probability"]) for k in row_indices])
        p_var = np.array([float(var_rows[k]["predicted_probability"]) for k in row_indices])
        if len(np.unique(y_base)) < 2:
            continue
        auc_b = float(roc_auc_score(y_base, p_base))
        auc_v = float(roc_auc_score(y_var, p_var))
        brier_b = float(brier_score_loss(y_base, p_base))
        brier_v = float(brier_score_loss(y_var, p_var))
        base_aucs.append(auc_b)
        var_aucs.append(auc_v)
        base_briers.append(brier_b)
        var_briers.append(brier_v)
        delta_aucs.append(auc_v - auc_b)
        delta_briers.append(brier_v - brier_b)

    point_base = compute_aggregate_auc_brier(base_rows)
    point_var = compute_aggregate_auc_brier(var_rows)

    def _ci(values: Sequence[float]) -> dict[str, Optional[float]]:
        if not values:
            return {"low": None, "high": None, "mean": None, "n_resamples": 0}
        arr = np.asarray(values, dtype=float)
        low, high = np.quantile(arr, [0.025, 0.975])
        return {
            "low": safe_round(float(low)),
            "high": safe_round(float(high)),
            "mean": safe_round(float(arr.mean())),
            "n_resamples": len(values),
        }

    return {
        "n_common_pairs": len(common_pairs),
        "baseline": {
            "auc_point": safe_round(point_base["auc"]),
            "brier_point": safe_round(point_base["brier"]),
            "auc_ci": _ci(base_aucs),
            "brier_ci": _ci(base_briers),
        },
        "variant": {
            "auc_point": safe_round(point_var["auc"]),
            "brier_point": safe_round(point_var["brier"]),
            "auc_ci": _ci(var_aucs),
            "brier_ci": _ci(var_briers),
        },
        "delta": {
            "auc_point": safe_round(point_var["auc"] - point_base["auc"]),
            "brier_point": safe_round(point_var["brier"] - point_base["brier"]),
            "auc_ci": _ci(delta_aucs),
            "brier_ci": _ci(delta_briers),
        },
    }


def prefilter_slot_topk_by_variance(
    context: Any,
    *,
    slot_name: str,
    feature_prefix: str,
    top_k: int,
) -> dict[str, int]:
    """Keep the top-K features ranked by entity-level variance; drop the rest.

    Rank-by-variance preserves both max-balanced features (support near 50% →
    Var=0.25) and narrow-support features concentrated in small cliques
    (support 10/148 → Var=0.063 — still informative). Drops near-constant
    features (support 1 of 369 → Var=0.003) which carry no discriminative
    signal for the task and just slow down RFE.

    Deterministic (sort is stable on feature name for ties), leakage-free
    (variance is over entities, not over training/test pairs), ablation-
    consistent (top_k is recorded alongside delta estimates). Returns
    `{"before": N, "after": K, "top_k": K}` for logging.
    """
    from lyzortx.autoresearch.train import SlotArtifact

    if slot_name not in context.slot_artifacts:
        raise KeyError(f"Slot {slot_name!r} not attached to context — call attach first")
    art = context.slot_artifacts[slot_name]
    frame = art.frame
    n_entities = len(frame)
    feature_cols = [c for c in frame.columns if c.startswith(feature_prefix)]
    before = len(feature_cols)
    if before <= top_k:
        LOGGER.info("Slot %s has %d features ≤ top_k=%d; keeping all", slot_name, before, top_k)
        return {"before": before, "after": before, "top_k": top_k, "n_entities": n_entities}
    variances = frame[feature_cols].var(axis=0).sort_values(ascending=False, kind="stable")
    kept = variances.head(top_k).index.tolist()
    after = len(kept)
    non_feature_cols = [c for c in frame.columns if c not in feature_cols]
    new_frame = frame[non_feature_cols + kept].copy()
    new_artifact = SlotArtifact(
        slot_name=slot_name,
        entity_key=art.entity_key,
        feature_columns=tuple(kept),
        frame=new_frame,
    )
    context.slot_artifacts[slot_name] = new_artifact
    LOGGER.info(
        "Pre-filtered slot %s by top-%d variance (n_entities=%d): %d → %d features",
        slot_name,
        top_k,
        n_entities,
        before,
        after,
    )
    return {"before": before, "after": after, "top_k": top_k, "n_entities": n_entities}


def load_canonical_baseline_predictions() -> pd.DataFrame:
    if not CH04_CANONICAL_PREDICTIONS.exists():
        raise FileNotFoundError(
            f"CH04 canonical predictions not found at {CH04_CANONICAL_PREDICTIONS}; "
            f"CH08 reuses these as the baseline arm. Run CH04 first."
        )
    df = pd.read_csv(CH04_CANONICAL_PREDICTIONS)
    LOGGER.info("Loaded CH04 canonical baseline predictions: %d rows", len(df))
    return df


def run_ch08_eval(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
    phage_kmer_slot_path: Path = DEFAULT_PHAGE_KMER_SLOT_PATH,
    host_omp_kmer_slot_path: Path = DEFAULT_HOST_OMP_KMER_SLOT_PATH,
    num_workers: int = 3,
    drop_high_titer_only_positives: bool = False,
    max_folds: Optional[int] = None,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.now(timezone.utc)
    LOGGER.info(
        "CH08 wave-2 re-audit starting at %s (num_workers=%d, drop_high_titer_only=%s)",
        start_time.isoformat(),
        num_workers,
        drop_high_titer_only_positives,
    )

    # Baseline reused from CH04 canonical.
    baseline_df = load_canonical_baseline_predictions()
    baseline_rows = baseline_df.to_dict(orient="records")

    # Arm SX12 — extend phage_slots with MORINIERE_SLOT_NAME.
    sx12_phage_slots = tuple(BASE_PHAGE_SLOTS) + (MORINIERE_SLOT_NAME,)

    def _attach_sx12_slot(ctx: Any) -> None:
        attach_moriniere_slot(ctx, phage_kmer_slot_path)
        prefilter_slot_topk_by_variance(
            ctx,
            slot_name=MORINIERE_SLOT_NAME,
            feature_prefix="phage_moriniere_kmer__",
            top_k=KMER_TOPK_PER_SLOT,
        )

    sx12_pair_preds = run_ch04_variant(
        arm_id=ARM_SX12,
        phage_slots=sx12_phage_slots,
        host_slots=tuple(BASE_HOST_SLOTS),
        slot_attach_fns=[_attach_sx12_slot],
        device_type=device_type,
        cache_dir=cache_dir,
        candidate_dir=candidate_dir,
        output_dir=output_dir,
        num_workers=num_workers,
        drop_high_titer_only_positives=drop_high_titer_only_positives,
    )
    if max_folds is not None:
        LOGGER.warning("max_folds not passed through — CH08 always runs all 10 folds per arm")

    sx12_report = _bootstrap_paired_by_bacterium(baseline_rows, sx12_pair_preds.to_dict(orient="records"))
    with open(output_dir / "ch08_sx12_delta.json", "w", encoding="utf-8") as f:
        json.dump({"arm": ARM_SX12, "extra_slot": MORINIERE_SLOT_NAME, **sx12_report}, f, indent=2)
    LOGGER.info(
        "SX12 delta: AUC %s (CI %s, %s), Brier %s (CI %s, %s), n_common=%d",
        sx12_report["delta"]["auc_point"],
        sx12_report["delta"]["auc_ci"]["low"],
        sx12_report["delta"]["auc_ci"]["high"],
        sx12_report["delta"]["brier_point"],
        sx12_report["delta"]["brier_ci"]["low"],
        sx12_report["delta"]["brier_ci"]["high"],
        sx12_report["n_common_pairs"],
    )

    # Arm SX13 marginal — extend host_slots with host_omp_kmer. Per plan.yml: "the
    # marginal arm only (cross_term, path1_cluster added no additional lift in SX13 and
    # can be skipped)".
    sx13_host_slots = tuple(BASE_HOST_SLOTS) + (HOST_OMP_KMER_SLOT,)

    def _attach_sx13_slot(ctx: Any) -> None:
        attach_csv_slot(
            context=ctx,
            features_path=host_omp_kmer_slot_path,
            slot_name=HOST_OMP_KMER_SLOT,
            feature_prefix=HOST_OMP_KMER_PREFIX,
            entity_key="bacteria",
        )
        prefilter_slot_topk_by_variance(
            ctx,
            slot_name=HOST_OMP_KMER_SLOT,
            feature_prefix=HOST_OMP_KMER_PREFIX,
            top_k=KMER_TOPK_PER_SLOT,
        )

    sx13_pair_preds = run_ch04_variant(
        arm_id=ARM_SX13,
        phage_slots=tuple(BASE_PHAGE_SLOTS),
        host_slots=sx13_host_slots,
        slot_attach_fns=[_attach_sx13_slot],
        device_type=device_type,
        cache_dir=cache_dir,
        candidate_dir=candidate_dir,
        output_dir=output_dir,
        num_workers=num_workers,
        drop_high_titer_only_positives=drop_high_titer_only_positives,
    )
    sx13_report = _bootstrap_paired_by_bacterium(baseline_rows, sx13_pair_preds.to_dict(orient="records"))
    with open(output_dir / "ch08_sx13_delta.json", "w", encoding="utf-8") as f:
        json.dump({"arm": ARM_SX13, "extra_slot": HOST_OMP_KMER_SLOT, **sx13_report}, f, indent=2)
    LOGGER.info(
        "SX13 delta: AUC %s (CI %s, %s), Brier %s (CI %s, %s), n_common=%d",
        sx13_report["delta"]["auc_point"],
        sx13_report["delta"]["auc_ci"]["low"],
        sx13_report["delta"]["auc_ci"]["high"],
        sx13_report["delta"]["brier_point"],
        sx13_report["delta"]["brier_ci"]["low"],
        sx13_report["delta"]["brier_ci"]["high"],
        sx13_report["n_common_pairs"],
    )

    summary = pd.DataFrame(
        [
            {
                "arm": ARM_BASELINE,
                "extra_slot": None,
                "auc_point": sx12_report["baseline"]["auc_point"],
                "brier_point": sx12_report["baseline"]["brier_point"],
                "auc_low": sx12_report["baseline"]["auc_ci"]["low"],
                "auc_high": sx12_report["baseline"]["auc_ci"]["high"],
                "brier_low": sx12_report["baseline"]["brier_ci"]["low"],
                "brier_high": sx12_report["baseline"]["brier_ci"]["high"],
                "delta_auc": None,
                "delta_auc_low": None,
                "delta_auc_high": None,
                "delta_brier": None,
                "delta_brier_low": None,
                "delta_brier_high": None,
            },
            {
                "arm": ARM_SX12,
                "extra_slot": MORINIERE_SLOT_NAME,
                "auc_point": sx12_report["variant"]["auc_point"],
                "brier_point": sx12_report["variant"]["brier_point"],
                "auc_low": sx12_report["variant"]["auc_ci"]["low"],
                "auc_high": sx12_report["variant"]["auc_ci"]["high"],
                "brier_low": sx12_report["variant"]["brier_ci"]["low"],
                "brier_high": sx12_report["variant"]["brier_ci"]["high"],
                "delta_auc": sx12_report["delta"]["auc_point"],
                "delta_auc_low": sx12_report["delta"]["auc_ci"]["low"],
                "delta_auc_high": sx12_report["delta"]["auc_ci"]["high"],
                "delta_brier": sx12_report["delta"]["brier_point"],
                "delta_brier_low": sx12_report["delta"]["brier_ci"]["low"],
                "delta_brier_high": sx12_report["delta"]["brier_ci"]["high"],
            },
            {
                "arm": ARM_SX13,
                "extra_slot": HOST_OMP_KMER_SLOT,
                "auc_point": sx13_report["variant"]["auc_point"],
                "brier_point": sx13_report["variant"]["brier_point"],
                "auc_low": sx13_report["variant"]["auc_ci"]["low"],
                "auc_high": sx13_report["variant"]["auc_ci"]["high"],
                "brier_low": sx13_report["variant"]["brier_ci"]["low"],
                "brier_high": sx13_report["variant"]["brier_ci"]["high"],
                "delta_auc": sx13_report["delta"]["auc_point"],
                "delta_auc_low": sx13_report["delta"]["auc_ci"]["low"],
                "delta_auc_high": sx13_report["delta"]["auc_ci"]["high"],
                "delta_brier": sx13_report["delta"]["brier_point"],
                "delta_brier_low": sx13_report["delta"]["brier_ci"]["low"],
                "delta_brier_high": sx13_report["delta"]["brier_ci"]["high"],
            },
        ]
    )
    summary.to_csv(output_dir / "ch08_summary.csv", index=False)
    LOGGER.info("CH08 summary:\n%s", summary.to_string(index=False))

    combined = {
        "task_id": "CH08",
        "scorecard": "AUC + Brier (paired bacterium-level bootstrap)",
        "drop_high_titer_only_positives": drop_high_titer_only_positives,
        "sx12": sx12_report,
        "sx13": sx13_report,
        "elapsed_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 1),
    }
    with open(output_dir / "ch08_combined_summary.json", "w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2)
    return combined


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--phage-kmer-slot-path", type=Path, default=DEFAULT_PHAGE_KMER_SLOT_PATH)
    parser.add_argument("--host-omp-kmer-slot-path", type=Path, default=DEFAULT_HOST_OMP_KMER_SLOT_PATH)
    parser.add_argument("--num-workers", type=int, default=3)
    parser.add_argument(
        "--drop-high-titer-only-positives",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Opt-in neat-only positive filter (CH10 demoted; default OFF).",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    run_ch08_eval(
        device_type=args.device_type,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
        phage_kmer_slot_path=args.phage_kmer_slot_path,
        host_omp_kmer_slot_path=args.host_omp_kmer_slot_path,
        num_workers=args.num_workers,
        drop_high_titer_only_positives=args.drop_high_titer_only_positives,
    )


if __name__ == "__main__":
    main()


__all__ = [
    "run_ch08_eval",
]
