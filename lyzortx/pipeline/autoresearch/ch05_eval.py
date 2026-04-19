#!/usr/bin/env python3
"""CH05: Unified Guelin+BASEL CHISEL evaluation — both axes, cross-source breakdown.

Successor to SX15 under the CHISEL frame. Runs CH04's per-row all-pairs pipeline on the
unified 148-phage × 369-bacteria panel:

- **Bacteria-axis** (10-fold via CH02 cv_group hash): all 148 phages remain in training;
  bacteria are held out. Measures generalization to unseen bacteria.
- **Phage-axis** (10-fold StratifiedKFold by ICTV family): all 369 bacteria remain in
  training; phages are held out. Measures generalization to unseen phages — the
  deployment question per `deployment-goal`.

**Per-phage blending is retired on both axes.** `per-phage-retired-under-chisel` applies
track-wide (not just to bacteria-axis). On phage-axis per-phage is structurally impossible
(held-out phages have zero training rows); on bacteria-axis the retirement rationale from
CH04 still holds — it is not deployable, inflates bacteria-axis AUC non-transferably, and
an intermediate CH04 diagnostic showed it deflates positive predictions under per-row
training. This deviates from the CH05 ticket's "BACTERIA-AXIS: per-phage blending enabled"
directive, justified per AGENTS.md "Question the requirement" and documented in
`per-phage-retired-under-chisel`. The ticket's "BLENDING TAX" is consequently retired too —
under all-pairs-only evaluation, the bacteria-vs-phage-axis gap is the plain phage-axis
generalization gap, not a blending tax.

**Cross-source sanity check**: on phage-axis, AUC is computed separately for held-out
Guelin phages (96) and held-out BASEL phages (52). SX15 under nDCG produced a BASEL-vs-
Guelin artifact (11 pp gap) that was a metric artifact, not a real generalization
difference. Under CHISEL AUC+Brier, the two cohorts should land within ~1 pp of each
other.

Usage:
    PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch05_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score
from sklearn.model_selection import StratifiedKFold

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    BootstrapMetricCI,
    load_module_from_path,
    safe_round,
)
from lyzortx.pipeline.autoresearch.ch03_row_expansion import (
    RAW_SCORE_NEGATIVE,
    RAW_SCORE_POSITIVE,
    RAW_SCORE_UNINTERPRETABLE,
    load_row_expanded_frame,
)
from lyzortx.pipeline.autoresearch.ch04_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    CONCENTRATION_FEATURE_COLUMN,
    build_clean_row_training_frame,
    compute_aggregate_auc_brier,
    select_pair_max_concentration_rows,
    train_and_predict_per_row_fold,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)
from lyzortx.pipeline.autoresearch.sx03_eval import (
    load_basel_interactions,
    patch_context_with_extended_slots,
)
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch05_unified_kfold")
CH04_AGGREGATE_METRICS_PATH = Path("lyzortx/generated_outputs/ch04_chisel_baseline/ch04_aggregate_metrics.json")

SOURCE_GUELIN = "guelin"
SOURCE_BASEL = "basel"
BASEL_REPLICATE_ID = 1
BASEL_LOG_DILUTION = 0  # BASEL's single spot test approximates Guelin's neat concentration.
PHAGE_AXIS_RANDOM_STATE = 42


def load_basel_as_row_frame(bacteria_cv_group: dict[str, str]) -> pd.DataFrame:
    """Materialise BASEL pairs as row-level training rows.

    BASEL is a single spot test (>10⁹ pfu/ml) per pair, not a dilution series, so each
    pair contributes exactly one training row. We encode it as `log_dilution=0` (matching
    Guelin's neat concentration) and `replicate=1`, with `score` copied from the binary
    interaction outcome.

    BASEL bacteria are a subset of Guelin bacteria (all 25 BASEL ECOR strains appear in
    the 369-bacterium Guelin panel), so `cv_group` is inherited from the Guelin mapping.
    Rows whose bacterium is absent from the Guelin cv_group map are dropped (should be
    zero in practice — the BASEL ECOR panel is fully contained in Guelin).
    """
    basel = load_basel_interactions().copy()
    basel = basel.dropna(subset=["interaction"])
    basel["score"] = basel["interaction"].astype(int).astype(str)
    basel["log_dilution"] = BASEL_LOG_DILUTION
    basel["replicate"] = BASEL_REPLICATE_ID
    basel["source"] = SOURCE_BASEL

    basel["cv_group"] = basel["bacteria"].map(bacteria_cv_group)
    missing = basel[basel["cv_group"].isna()]
    if not missing.empty:
        LOGGER.warning(
            "BASEL rows with no Guelin cv_group (dropped): %d (example bacteria: %s)",
            len(missing),
            missing["bacteria"].unique()[:5].tolist(),
        )
        basel = basel.dropna(subset=["cv_group"])

    # Downstream per-row training frame filters on split_holdout and is_hard_trainable.
    # Treat every BASEL observation as trainable — no ST03-style holdout split applies to
    # BASEL (the CHISEL unified evaluation handles splits via k-fold CV, not ST03).
    basel["split_holdout"] = "train_non_holdout"
    basel["is_hard_trainable"] = "1"

    LOGGER.info(
        "BASEL row frame: %d rows (%d positive), %d pairs, %d bacteria, %d phages",
        len(basel),
        int((basel["score"] == RAW_SCORE_POSITIVE).sum()),
        basel["pair_id"].nunique(),
        basel["bacteria"].nunique(),
        basel["phage"].nunique(),
    )
    return basel[
        [
            "pair_id",
            "bacteria",
            "phage",
            "log_dilution",
            "replicate",
            "score",
            "source",
            "cv_group",
            "split_holdout",
            "is_hard_trainable",
        ]
    ].copy()


def load_unified_row_frame() -> pd.DataFrame:
    """Concatenate Guelin (318K rows, 9 replicates × pair) + BASEL (~1.2K rows, 1 row × pair).

    Both sides carry `source ∈ {"guelin", "basel"}` and the minimal column set CH04's
    per-row pipeline consumes. Guelin rows inherit split_holdout / is_hard_trainable from
    ST02+ST03; BASEL rows are marked universally trainable.
    """
    guelin = load_row_expanded_frame().copy()
    guelin["source"] = SOURCE_GUELIN

    # Extract Guelin bacteria → cv_group so BASEL rows can inherit cv_group.
    guelin_map = bacteria_to_cv_group_map(guelin)
    basel = load_basel_as_row_frame(guelin_map)

    # Reduce both frames to the shared column set before concatenation. Guelin has many
    # additional pair-level columns (host/phage metadata, training_weight_v3, etc.) that
    # BASEL lacks — downstream code only reads the columns below plus anything derived
    # from the feature-slot cache, so keep it minimal.
    keep = [
        "pair_id",
        "bacteria",
        "phage",
        "log_dilution",
        "replicate",
        "score",
        "source",
        "cv_group",
        "split_holdout",
        "is_hard_trainable",
    ]
    unified = pd.concat([guelin[keep], basel[keep]], ignore_index=True)
    LOGGER.info(
        "Unified row frame: %d total rows (%d guelin, %d basel); %d pairs, %d bacteria, %d phages",
        len(unified),
        (unified["source"] == SOURCE_GUELIN).sum(),
        (unified["source"] == SOURCE_BASEL).sum(),
        unified["pair_id"].nunique(),
        unified["bacteria"].nunique(),
        unified["phage"].nunique(),
    )
    return unified


def assign_phage_folds(
    phages: Sequence[str],
    phage_family: dict[str, str],
    n_splits: int = N_FOLDS,
    random_state: int = PHAGE_AXIS_RANDOM_STATE,
) -> dict[str, int]:
    """10-fold StratifiedKFold on the phage list, stratified by ICTV family.

    Rare families (fewer than `n_splits` phages) are collapsed into an "other" bucket so
    StratifiedKFold can split them. Assignment is deterministic given `random_state`.
    Phages with no family annotation are mapped to "UNKNOWN" (and fall into the "other"
    bucket if that stratum has fewer than n_splits members).
    """
    phage_list = sorted(set(str(p) for p in phages))
    df = pd.DataFrame({"phage": phage_list})
    df["family"] = df["phage"].map(phage_family).fillna("UNKNOWN")
    counts = df["family"].value_counts()
    rare = counts[counts < n_splits].index
    df.loc[df["family"].isin(rare), "family"] = "other"

    skf = StratifiedKFold(n_splits=n_splits, shuffle=True, random_state=random_state)
    fold_map: dict[str, int] = {}
    for fold_id, (_, test_idx) in enumerate(skf.split(df["phage"].to_numpy(), df["family"].to_numpy())):
        for phage in df.iloc[test_idx]["phage"].to_list():
            fold_map[phage] = fold_id
    per_fold_counts = defaultdict(int)
    for fold in fold_map.values():
        per_fold_counts[fold] += 1
    LOGGER.info(
        "phage-axis folds (n=%d): sizes=%s (stratified by %d unique families after collapsing rares)",
        len(phage_list),
        dict(sorted(per_fold_counts.items())),
        df["family"].nunique(),
    )
    return fold_map


def _aggregate_fold_rows_to_pairs(fold_rows: list[dict[str, object]]) -> pd.DataFrame:
    """Collapse per-seed row predictions to mean-over-seeds then pick pair-max-conc rows."""
    df = pd.DataFrame(fold_rows)
    aggregated = (
        df.groupby(
            [
                "fold_id",
                "pair_id",
                "bacteria",
                "phage",
                "log_dilution",
                "replicate",
                "label_row_binary",
            ],
            as_index=False,
        )["predicted_probability"]
        .mean()
        .sort_values(["bacteria", "phage", "log_dilution", "replicate"])
    )
    return aggregated


def _bootstrap_by_unit(
    pair_rows: Sequence[dict[str, object]],
    *,
    unit_key: str,
    bootstrap_samples: int,
    bootstrap_random_state: int,
) -> dict[str, BootstrapMetricCI]:
    """Resample `unit_key` (bacterium or phage) with replacement and compute AUC+Brier per resample."""
    grouped: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        grouped[str(row[unit_key])].append(row)
    unit_ids = tuple(sorted(grouped.keys()))
    n_units = len(unit_ids)
    rng = np.random.default_rng(bootstrap_random_state)

    aucs: list[float] = []
    briers: list[float] = []
    progress_interval = max(1, bootstrap_samples // 5)
    for i in range(bootstrap_samples):
        if i == 0 or (i + 1) % progress_interval == 0 or i + 1 == bootstrap_samples:
            LOGGER.info("Bootstrap progress (%s): %d/%d", unit_key, i + 1, bootstrap_samples)
        idx = rng.integers(0, n_units, size=n_units)
        sampled: list[dict[str, object]] = []
        for j in idx.tolist():
            sampled.extend(grouped[unit_ids[j]])
        labels = np.array([int(r["label_row_binary"]) for r in sampled])
        preds = np.array([float(r["predicted_probability"]) for r in sampled])
        if len(np.unique(labels)) < 2:
            continue
        aucs.append(float(roc_auc_score(labels, preds)))
        briers.append(float(brier_score_loss(labels, preds)))

    point = compute_aggregate_auc_brier(pair_rows)

    def _ci(values: Sequence[float]) -> tuple[Optional[float], Optional[float], int]:
        if not values:
            return None, None, 0
        arr = np.asarray(values, dtype=float)
        low, high = np.quantile(arr, [0.025, 0.975])
        return safe_round(float(low)), safe_round(float(high)), len(values)

    auc_low, auc_high, auc_used = _ci(aucs)
    brier_low, brier_high, brier_used = _ci(briers)
    return {
        "holdout_roc_auc": BootstrapMetricCI(
            point_estimate=point["auc"],
            ci_low=auc_low,
            ci_high=auc_high,
            bootstrap_samples_requested=bootstrap_samples,
            bootstrap_samples_used=auc_used,
        ),
        "holdout_brier_score": BootstrapMetricCI(
            point_estimate=point["brier"],
            ci_low=brier_low,
            ci_high=brier_high,
            bootstrap_samples_requested=bootstrap_samples,
            bootstrap_samples_used=brier_used,
        ),
    }


def _ci_to_dict(ci: BootstrapMetricCI) -> dict[str, Optional[float]]:
    return {
        "point_estimate": safe_round(ci.point_estimate) if ci.point_estimate is not None else None,
        "ci_low": ci.ci_low,
        "ci_high": ci.ci_high,
    }


def run_bacteria_axis(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_rows: pd.DataFrame,
    device_type: str,
    output_dir: Path,
) -> pd.DataFrame:
    """10-fold bacteria-axis CV via CH02 cv_group hash (all 148 phages in training)."""
    LOGGER.info("=== Bacteria-axis evaluation (CH02 cv_group folds, all-pairs only) ===")
    mapping = bacteria_to_cv_group_map(clean_rows)
    fold_assignments = assign_bacteria_folds(mapping)
    all_rows: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_bac = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bac = {b for b, f in fold_assignments.items() if f != fold_id}
        fold_train = clean_rows[clean_rows["bacteria"].isin(train_bac)].copy()
        fold_holdout = clean_rows[clean_rows["bacteria"].isin(holdout_bac)].copy()
        LOGGER.info(
            "=== bacteria-axis Fold %d: %d train bacteria (%d rows), %d holdout bacteria (%d rows) ===",
            fold_id,
            len(train_bac),
            len(fold_train),
            len(holdout_bac),
            len(fold_holdout),
        )
        fold_seed_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows, _ = train_and_predict_per_row_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=fold_train,
                holdout_frame=fold_holdout,
                seed=seed,
                device_type=device_type,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_seed_rows.extend(rows)
        aggregated = _aggregate_fold_rows_to_pairs(fold_seed_rows)
        pair_pred = select_pair_max_concentration_rows(aggregated)
        fold_point = compute_aggregate_auc_brier(pair_pred.to_dict(orient="records"))
        LOGGER.info(
            "bacteria-axis Fold %d: AUC=%.4f, Brier=%.4f, n_pairs=%d",
            fold_id,
            fold_point["auc"],
            fold_point["brier"],
            len(pair_pred),
        )
        all_rows.extend(aggregated.to_dict(orient="records"))
    per_row_df = pd.DataFrame(all_rows)
    per_row_df.to_csv(output_dir / "ch05_bacteria_axis_per_row_predictions.csv", index=False)
    return per_row_df


def run_phage_axis(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_rows: pd.DataFrame,
    phage_family: dict[str, str],
    device_type: str,
    output_dir: Path,
) -> pd.DataFrame:
    """10-fold phage-axis CV via StratifiedKFold (all 369 bacteria in training)."""
    LOGGER.info("=== Phage-axis evaluation (StratifiedKFold by ICTV family, all-pairs only) ===")
    phages = sorted(clean_rows["phage"].unique())
    fold_assignments = assign_phage_folds(phages, phage_family)
    all_rows: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_phages = {p for p, f in fold_assignments.items() if f == fold_id}
        train_phages = {p for p, f in fold_assignments.items() if f != fold_id}
        fold_train = clean_rows[clean_rows["phage"].isin(train_phages)].copy()
        fold_holdout = clean_rows[clean_rows["phage"].isin(holdout_phages)].copy()
        LOGGER.info(
            "=== phage-axis Fold %d: %d train phages (%d rows), %d holdout phages (%d rows) ===",
            fold_id,
            len(train_phages),
            len(fold_train),
            len(holdout_phages),
            len(fold_holdout),
        )
        fold_seed_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows, _ = train_and_predict_per_row_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=fold_train,
                holdout_frame=fold_holdout,
                seed=seed,
                device_type=device_type,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_seed_rows.extend(rows)
        aggregated = _aggregate_fold_rows_to_pairs(fold_seed_rows)
        pair_pred = select_pair_max_concentration_rows(aggregated)
        fold_point = compute_aggregate_auc_brier(pair_pred.to_dict(orient="records"))
        LOGGER.info(
            "phage-axis Fold %d: AUC=%.4f, Brier=%.4f, n_pairs=%d",
            fold_id,
            fold_point["auc"],
            fold_point["brier"],
            len(pair_pred),
        )
        all_rows.extend(aggregated.to_dict(orient="records"))
    per_row_df = pd.DataFrame(all_rows)
    per_row_df.to_csv(output_dir / "ch05_phage_axis_per_row_predictions.csv", index=False)
    return per_row_df


def run_ch05_eval(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.now(timezone.utc)
    LOGGER.info("CH05 evaluation starting at %s", start_time.isoformat())

    unified = load_unified_row_frame()
    clean_rows = build_clean_row_training_frame(unified)
    LOGGER.info(
        "CH05 clean row frame: %d rows, %d pairs, %d bacteria, %d phages",
        len(clean_rows),
        clean_rows["pair_id"].nunique(),
        clean_rows["bacteria"].nunique(),
        clean_rows["phage"].nunique(),
    )
    # Carry source forward so cross-source breakdown works after CH04 pipeline strips columns.
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]

    phage_family = load_unified_phage_family_map()
    candidate_module = load_module_from_path("ch05_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context)

    bacteria_axis_per_row = run_bacteria_axis(
        candidate_module=candidate_module,
        context=context,
        clean_rows=clean_rows,
        device_type=device_type,
        output_dir=output_dir,
    )
    phage_axis_per_row = run_phage_axis(
        candidate_module=candidate_module,
        context=context,
        clean_rows=clean_rows,
        phage_family=phage_family,
        device_type=device_type,
        output_dir=output_dir,
    )

    # Per-axis aggregate + bootstrap.
    bacteria_pairs = select_pair_max_concentration_rows(bacteria_axis_per_row)
    bacteria_pairs["source"] = bacteria_pairs["pair_id"].map(pair_source)
    bacteria_pairs.to_csv(output_dir / "ch05_bacteria_axis_predictions.csv", index=False)
    bacteria_rows = bacteria_pairs.to_dict(orient="records")
    bacteria_cis = _bootstrap_by_unit(
        bacteria_rows,
        unit_key="bacteria",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )
    bacteria_summary = {
        "n_pairs": len(bacteria_pairs),
        "n_bacteria": int(bacteria_pairs["bacteria"].nunique()),
        "n_phages": int(bacteria_pairs["phage"].nunique()),
        "aggregate": {name: _ci_to_dict(ci) for name, ci in bacteria_cis.items()},
    }
    with open(output_dir / "ch05_bacteria_axis_metrics.json", "w", encoding="utf-8") as f:
        json.dump(bacteria_summary, f, indent=2)
    LOGGER.info(
        "Bacteria-axis aggregate: AUC %.4f [%.4f, %.4f], Brier %.4f [%.4f, %.4f]",
        bacteria_cis["holdout_roc_auc"].point_estimate,
        bacteria_cis["holdout_roc_auc"].ci_low or 0,
        bacteria_cis["holdout_roc_auc"].ci_high or 0,
        bacteria_cis["holdout_brier_score"].point_estimate,
        bacteria_cis["holdout_brier_score"].ci_low or 0,
        bacteria_cis["holdout_brier_score"].ci_high or 0,
    )

    phage_pairs = select_pair_max_concentration_rows(phage_axis_per_row)
    phage_pairs["source"] = phage_pairs["pair_id"].map(pair_source)
    phage_pairs.to_csv(output_dir / "ch05_phage_axis_predictions.csv", index=False)
    phage_rows_all = phage_pairs.to_dict(orient="records")
    phage_cis_all = _bootstrap_by_unit(
        phage_rows_all,
        unit_key="phage",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    cross_source_rows: list[dict[str, object]] = []
    for source_label in (SOURCE_GUELIN, SOURCE_BASEL):
        subset = [r for r in phage_rows_all if r["source"] == source_label]
        LOGGER.info(
            "Phage-axis cross-source subset %s: %d pairs, %d phages",
            source_label,
            len(subset),
            len({r["phage"] for r in subset}),
        )
        if not subset:
            continue
        subset_cis = _bootstrap_by_unit(
            subset,
            unit_key="phage",
            bootstrap_samples=BOOTSTRAP_SAMPLES,
            bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
        )
        cross_source_rows.append(
            {
                "source": source_label,
                "n_pairs": len(subset),
                "n_phages": len({r["phage"] for r in subset}),
                "auc_point": safe_round(subset_cis["holdout_roc_auc"].point_estimate),
                "auc_low": subset_cis["holdout_roc_auc"].ci_low,
                "auc_high": subset_cis["holdout_roc_auc"].ci_high,
                "brier_point": safe_round(subset_cis["holdout_brier_score"].point_estimate),
                "brier_low": subset_cis["holdout_brier_score"].ci_low,
                "brier_high": subset_cis["holdout_brier_score"].ci_high,
            }
        )
    cross_source_df = pd.DataFrame(cross_source_rows)
    cross_source_df.to_csv(output_dir / "ch05_cross_source_breakdown.csv", index=False)
    LOGGER.info("Cross-source breakdown (phage-axis):\n%s", cross_source_df.to_string(index=False))

    phage_summary = {
        "n_pairs": len(phage_pairs),
        "n_bacteria": int(phage_pairs["bacteria"].nunique()),
        "n_phages": int(phage_pairs["phage"].nunique()),
        "aggregate": {name: _ci_to_dict(ci) for name, ci in phage_cis_all.items()},
        "cross_source": cross_source_rows,
    }
    with open(output_dir / "ch05_phage_axis_metrics.json", "w", encoding="utf-8") as f:
        json.dump(phage_summary, f, indent=2)

    # "Blending tax" is retired; under all-pairs-only this becomes the plain phage-axis
    # generalization gap.
    phage_axis_gap = safe_round(
        bacteria_cis["holdout_roc_auc"].point_estimate - phage_cis_all["holdout_roc_auc"].point_estimate
    )
    cross_source_gap = None
    if cross_source_df.shape[0] == 2:
        guelin_auc = cross_source_df.loc[cross_source_df["source"] == SOURCE_GUELIN, "auc_point"].iloc[0]
        basel_auc = cross_source_df.loc[cross_source_df["source"] == SOURCE_BASEL, "auc_point"].iloc[0]
        cross_source_gap = safe_round(abs(guelin_auc - basel_auc))

    ch04_auc = None
    if CH04_AGGREGATE_METRICS_PATH.exists():
        with open(CH04_AGGREGATE_METRICS_PATH, encoding="utf-8") as f:
            ch04 = json.load(f)
        ch04_auc = ch04["chisel_baseline"]["holdout_roc_auc"]["point_estimate"]

    combined = {
        "task_id": "CH05",
        "scorecard": "AUC + Brier (nDCG/mAP/top-k retired)",
        "per_phage_blending": "retired (see per-phage-retired-under-chisel)",
        "bacteria_axis": bacteria_summary,
        "phage_axis": phage_summary,
        "phage_axis_generalization_gap_auc": phage_axis_gap,
        "cross_source_auc_delta": cross_source_gap,
        "ch04_baseline_auc": ch04_auc,
        "elapsed_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 1),
    }
    with open(output_dir / "ch05_combined_summary.json", "w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2)
    LOGGER.info(
        "Phage-axis generalization gap (bacteria_auc − phage_auc) = %.4f; cross-source |ΔAUC| = %s",
        phage_axis_gap,
        cross_source_gap,
    )
    return combined


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
    run_ch05_eval(
        device_type=args.device_type,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
    )


if __name__ == "__main__":
    main()


__all__ = [
    "assign_phage_folds",
    "load_basel_as_row_frame",
    "load_unified_row_frame",
    "run_ch05_eval",
    "CONCENTRATION_FEATURE_COLUMN",
    "RAW_SCORE_NEGATIVE",
    "RAW_SCORE_POSITIVE",
    "RAW_SCORE_UNINTERPRETABLE",
]
