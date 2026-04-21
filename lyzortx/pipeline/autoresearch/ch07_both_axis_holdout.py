#!/usr/bin/env python3
"""CH07: Both-axis 10×10 double cross-validation — bacterium AND phage held out.

The strongest generalization test in the CHISEL track: in each of 100 train/test
cells, both the held-out bacteria fold AND the held-out phage fold are removed
from training simultaneously. Addresses the external-critic concern that the
"panel-size ceiling" (see `panel-size-ceiling` knowledge unit) has so far only
been probed on single-axis CV (CH04 bacteria-axis, CH05 phage-axis).

Design:
  - 10 bacteria folds (CH02 cv_group hash, same as CH04/CH05 bacteria-axis)
  - 10 phage folds (StratifiedKFold by ICTV family, same as CH05 phage-axis)
  - 100 cells: each trained on (all-but-this-bact-fold) × (all-but-this-phage-fold)
    and evaluated on the held-out bact × held-out phage cell

Phage-side feature slot: **Arm 3 (Moriniere per-receptor k-mer fractions)**. The
plan.yml direction — "the winning phage-side feature bundle from CH06, not the
TL17 baseline that CH05 diagnosed as broken". Arm 3 was validated in CH06 and
promoted in the `moriniere-receptor-fractions-validated` knowledge unit.

All-pairs only (per-phage blending retired, see `per-phage-retired-under-chisel`).

Artifacts (under `lyzortx/generated_outputs/ch07_both_axis_holdout/`):
  - ch07_cell_metrics.csv — per-cell AUC, Brier, n_pairs, n_bacteria, n_phages
  - ch07_per_row_predictions.csv — pooled per-row predictions, all 100 cells
  - ch07_pair_predictions.csv — pooled pair-level (max-conc) predictions
  - ch07_aggregate.json — pooled AUC + Brier with pair-level bootstrap CIs
  - ch07_cell_distribution.png — histogram of per-cell AUC

Usage:
    PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch07_both_axis_holdout --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    BootstrapMetricCI,
    load_module_from_path,
    safe_round,
)
from lyzortx.pipeline.autoresearch.ch04_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    build_clean_row_training_frame,
    compute_aggregate_auc_brier,
    select_pair_max_concentration_rows,
)
from lyzortx.pipeline.autoresearch.ch04_parallel import (
    fit_seeds,
    prepare_fold_design_matrices,
    select_rfe_features,
)
from lyzortx.pipeline.autoresearch.ch05_eval import (
    BASEL_LOG10_PFU_ML,
    PHAGE_AXIS_RANDOM_STATE,
    assign_phage_folds,
    load_unified_row_frame,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)
from lyzortx.pipeline.autoresearch.sx03_eval import patch_context_with_extended_slots
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch07_both_axis_holdout")
# Arm 3 slot directory — Moriniere per-receptor k-mer fractions (13-dim), validated
# in CH06 as the panel-independent replacement for the Guelin-derived TL17
# phage_projection slot. See moriniere-receptor-fractions-validated.
ARM3_SLOTS_DIR = Path(".scratch/basel/feature_slots_arm3")

CH05_COMBINED_METRICS_PATH = Path("lyzortx/generated_outputs/ch05_unified_kfold/ch05_combined_summary.json")


def _cells_plan(n_bact_folds: int, n_phage_folds: int) -> list[tuple[int, int]]:
    return [(b, p) for b in range(n_bact_folds) for p in range(n_phage_folds)]


def _bootstrap_pair_level(
    pair_rows: Sequence[dict[str, object]],
    *,
    bootstrap_samples: int,
    bootstrap_random_state: int,
) -> dict[str, BootstrapMetricCI]:
    """Pair-level bootstrap: resample held-out (bacterium, phage) pairs with replacement.

    Pair-level resampling matches the CH07 design decision in plan.yml: per-cell pair
    counts are too small to bootstrap at the cell level, so the 100-cell aggregate
    is bootstrapped over the full pooled pair list. Each draw picks `n_pairs` rows
    with replacement and computes AUC + Brier on that resample.
    """
    rows = list(pair_rows)
    n = len(rows)
    rng = np.random.default_rng(bootstrap_random_state)
    aucs: list[float] = []
    briers: list[float] = []
    progress_interval = max(1, bootstrap_samples // 5)
    for i in range(bootstrap_samples):
        if i == 0 or (i + 1) % progress_interval == 0 or i + 1 == bootstrap_samples:
            LOGGER.info("Pair-level bootstrap progress: %d/%d", i + 1, bootstrap_samples)
        idx = rng.integers(0, n, size=n)
        sampled = [rows[j] for j in idx.tolist()]
        labels = np.array([int(r["label_row_binary"]) for r in sampled])
        preds = np.array([float(r["predicted_probability"]) for r in sampled])
        if len(np.unique(labels)) < 2:
            continue
        aucs.append(float(roc_auc_score(labels, preds)))
        briers.append(float(brier_score_loss(labels, preds)))

    point = compute_aggregate_auc_brier(rows)

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
        "bootstrap_samples_requested": ci.bootstrap_samples_requested,
        "bootstrap_samples_used": ci.bootstrap_samples_used,
    }


def _aggregate_cell_rows_to_pairs(cell_rows: list[dict[str, object]]) -> pd.DataFrame:
    df = pd.DataFrame(cell_rows)
    aggregated = (
        df.groupby(
            [
                "bacteria_fold",
                "phage_fold",
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


def run_ch07_eval(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
    phage_slots_dir: Path = ARM3_SLOTS_DIR,
    basel_log10_pfu_ml: float = BASEL_LOG10_PFU_ML,
    num_workers: int = 3,
    drop_high_titer_only_positives: bool = False,
    max_cells: Optional[int] = None,
) -> dict[str, object]:
    """Run CH07 100-cell double cross-validation.

    `max_cells` limits the number of cells (subset iteration for engineering). The
    full 10×10 run evaluates 100 cells; smaller values iterate row-major
    (bacteria fold outer, phage fold inner).
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.now(timezone.utc)
    LOGGER.info(
        "CH07 evaluation starting at %s (phage_slots_dir=%s, basel_log10_pfu_ml=%.2f, "
        "num_workers=%d, drop_high_titer_only=%s, max_cells=%s)",
        start_time.isoformat(),
        phage_slots_dir,
        basel_log10_pfu_ml,
        num_workers,
        drop_high_titer_only_positives,
        max_cells if max_cells is not None else "100",
    )

    unified = load_unified_row_frame(basel_log10_pfu_ml=basel_log10_pfu_ml)
    clean_rows = build_clean_row_training_frame(unified, drop_high_titer_only_positives=drop_high_titer_only_positives)
    LOGGER.info(
        "CH07 clean row frame: %d rows, %d pairs, %d bacteria, %d phages",
        len(clean_rows),
        clean_rows["pair_id"].nunique(),
        clean_rows["bacteria"].nunique(),
        clean_rows["phage"].nunique(),
    )
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]

    # Bacteria folds (CH02 hashing) + phage folds (StratifiedKFold by family).
    bac_mapping = bacteria_to_cv_group_map(clean_rows)
    bacteria_fold = assign_bacteria_folds(bac_mapping)
    phage_family = load_unified_phage_family_map()
    phages = sorted(clean_rows["phage"].unique())
    phage_fold = assign_phage_folds(phages, phage_family, random_state=PHAGE_AXIS_RANDOM_STATE)

    clean_rows = clean_rows.copy()
    clean_rows["bacteria_fold"] = clean_rows["bacteria"].map(bacteria_fold).astype(int)
    clean_rows["phage_fold"] = clean_rows["phage"].map(phage_fold).astype(int)

    # Candidate module + context, patched with Arm 3 phage slot.
    candidate_module = load_module_from_path("ch07_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context, slots_dir=phage_slots_dir)

    cell_plan = _cells_plan(N_FOLDS, N_FOLDS)
    if max_cells is not None:
        cell_plan = cell_plan[:max_cells]
    total_cells = len(cell_plan)

    per_cell_metrics: list[dict[str, object]] = []
    all_per_row: list[dict[str, object]] = []
    all_pair: list[dict[str, object]] = []

    cell_start = datetime.now(timezone.utc)
    for cell_idx, (bact_fold_id, phage_fold_id) in enumerate(cell_plan, start=1):
        train = clean_rows[
            (clean_rows["bacteria_fold"] != bact_fold_id) & (clean_rows["phage_fold"] != phage_fold_id)
        ].copy()
        test = clean_rows[
            (clean_rows["bacteria_fold"] == bact_fold_id) & (clean_rows["phage_fold"] == phage_fold_id)
        ].copy()
        n_holdout_bact = test["bacteria"].nunique()
        n_holdout_phage = test["phage"].nunique()
        LOGGER.info(
            "=== Cell %d/%d (bact_fold=%d, phage_fold=%d): train=%d rows, test=%d rows (%d bacteria × %d phages) ===",
            cell_idx,
            total_cells,
            bact_fold_id,
            phage_fold_id,
            len(train),
            len(test),
            n_holdout_bact,
            n_holdout_phage,
        )
        if len(test) == 0:
            LOGGER.warning("Cell has zero holdout rows — skipping (no held-out phages × held-out bacteria)")
            continue

        train_design, holdout_design, feature_columns, categorical_columns = prepare_fold_design_matrices(
            candidate_module=candidate_module,
            context=context,
            training_frame=train,
            holdout_frame=test,
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
        cell_seed_rows: list[dict[str, object]] = []
        for _seed, rows, _fi in seed_results:
            for r in rows:
                r["bacteria_fold"] = bact_fold_id
                r["phage_fold"] = phage_fold_id
            cell_seed_rows.extend(rows)
        aggregated = _aggregate_cell_rows_to_pairs(cell_seed_rows)
        pair_pred = select_pair_max_concentration_rows(aggregated)
        pair_pred["bacteria_fold"] = bact_fold_id
        pair_pred["phage_fold"] = phage_fold_id
        cell_point = compute_aggregate_auc_brier(pair_pred.to_dict(orient="records"))
        per_cell_metrics.append(
            {
                "bacteria_fold": bact_fold_id,
                "phage_fold": phage_fold_id,
                "n_pairs": len(pair_pred),
                "n_bacteria": n_holdout_bact,
                "n_phages": n_holdout_phage,
                "positive_rate": safe_round(float(pair_pred["label_row_binary"].mean())),
                "auc": safe_round(cell_point["auc"]) if cell_point["auc"] is not None else None,
                "brier": safe_round(cell_point["brier"]),
            }
        )
        all_per_row.extend(aggregated.to_dict(orient="records"))
        all_pair.extend(pair_pred.to_dict(orient="records"))

        elapsed = (datetime.now(timezone.utc) - cell_start).total_seconds()
        cells_per_sec = cell_idx / elapsed if elapsed > 0 else 0.0
        remaining = (total_cells - cell_idx) / cells_per_sec if cells_per_sec > 0 else float("inf")
        LOGGER.info(
            "Cell %d/%d done: AUC=%s, Brier=%.4f, n_pairs=%d. Elapsed=%.0fs, ETA=%.0fs",
            cell_idx,
            total_cells,
            "nan" if cell_point["auc"] is None else f"{cell_point['auc']:.4f}",
            cell_point["brier"],
            len(pair_pred),
            elapsed,
            remaining,
        )

    # Persist per-cell and per-row / per-pair artifacts.
    cell_df = pd.DataFrame(per_cell_metrics)
    cell_df.to_csv(output_dir / "ch07_cell_metrics.csv", index=False)
    per_row_df = pd.DataFrame(all_per_row)
    per_row_df.to_csv(output_dir / "ch07_per_row_predictions.csv", index=False)
    pair_df = pd.DataFrame(all_pair)
    pair_df["source"] = pair_df["pair_id"].map(pair_source)
    pair_df.to_csv(output_dir / "ch07_pair_predictions.csv", index=False)

    # Pooled aggregate + pair-level bootstrap CI.
    pair_rows_list = pair_df.to_dict(orient="records")
    cis = _bootstrap_pair_level(
        pair_rows_list,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    # Cross-source decomposition on the pooled pair predictions.
    cross_source_rows: list[dict[str, object]] = []
    for source_label in ("guelin", "basel"):
        subset = [r for r in pair_rows_list if r.get("source") == source_label]
        if not subset:
            LOGGER.warning("Cross-source subset %s is empty — skipping", source_label)
            continue
        sub_cis = _bootstrap_pair_level(
            subset,
            bootstrap_samples=BOOTSTRAP_SAMPLES,
            bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
        )
        cross_source_rows.append(
            {
                "source": source_label,
                "n_pairs": len(subset),
                "n_phages": len({r["phage"] for r in subset}),
                "n_bacteria": len({r["bacteria"] for r in subset}),
                "auc_point": safe_round(sub_cis["holdout_roc_auc"].point_estimate),
                "auc_low": sub_cis["holdout_roc_auc"].ci_low,
                "auc_high": sub_cis["holdout_roc_auc"].ci_high,
                "brier_point": safe_round(sub_cis["holdout_brier_score"].point_estimate),
                "brier_low": sub_cis["holdout_brier_score"].ci_low,
                "brier_high": sub_cis["holdout_brier_score"].ci_high,
            }
        )
    cross_df = pd.DataFrame(cross_source_rows)
    cross_df.to_csv(output_dir / "ch07_cross_source_breakdown.csv", index=False)

    ch05_auc = None
    if CH05_COMBINED_METRICS_PATH.exists():
        with open(CH05_COMBINED_METRICS_PATH, encoding="utf-8") as f:
            ch05 = json.load(f)
        ch05_auc = ch05["bacteria_axis"]["aggregate"]["holdout_roc_auc"]["point_estimate"]

    # Cell-AUC distribution summary.
    valid_cell_aucs = [m["auc"] for m in per_cell_metrics if m["auc"] is not None]
    cell_stats = {
        "n_cells": len(per_cell_metrics),
        "n_cells_with_auc": len(valid_cell_aucs),
        "cell_auc_mean": safe_round(float(np.mean(valid_cell_aucs))) if valid_cell_aucs else None,
        "cell_auc_median": safe_round(float(np.median(valid_cell_aucs))) if valid_cell_aucs else None,
        "cell_auc_std": safe_round(float(np.std(valid_cell_aucs))) if valid_cell_aucs else None,
        "cell_auc_min": safe_round(float(np.min(valid_cell_aucs))) if valid_cell_aucs else None,
        "cell_auc_max": safe_round(float(np.max(valid_cell_aucs))) if valid_cell_aucs else None,
    }

    aggregate = {
        "task_id": "CH07",
        "scorecard": "AUC + Brier (pair-level bootstrap)",
        "phage_slots_dir": str(phage_slots_dir),
        "drop_high_titer_only_positives": drop_high_titer_only_positives,
        "n_pairs": len(pair_df),
        "n_bacteria": int(pair_df["bacteria"].nunique()),
        "n_phages": int(pair_df["phage"].nunique()),
        "aggregate": {name: _ci_to_dict(ci) for name, ci in cis.items()},
        "cross_source": cross_source_rows,
        "cell_distribution": cell_stats,
        "ch05_bacteria_axis_baseline_auc": ch05_auc,
        "elapsed_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 1),
    }
    with open(output_dir / "ch07_aggregate.json", "w", encoding="utf-8") as f:
        json.dump(aggregate, f, indent=2)

    # Per-cell AUC histogram.
    _write_cell_distribution_plot(per_cell_metrics, output_dir / "ch07_cell_distribution.png")

    LOGGER.info(
        "CH07 aggregate: AUC %.4f [%.4f, %.4f], Brier %.4f [%.4f, %.4f], %d cells",
        cis["holdout_roc_auc"].point_estimate,
        cis["holdout_roc_auc"].ci_low or 0,
        cis["holdout_roc_auc"].ci_high or 0,
        cis["holdout_brier_score"].point_estimate,
        cis["holdout_brier_score"].ci_low or 0,
        cis["holdout_brier_score"].ci_high or 0,
        len(per_cell_metrics),
    )
    if cross_source_rows:
        LOGGER.info("Cross-source breakdown:\n%s", cross_df.to_string(index=False))
    LOGGER.info(
        "Per-cell AUC distribution: mean=%s, median=%s, std=%s, min=%s, max=%s",
        cell_stats["cell_auc_mean"],
        cell_stats["cell_auc_median"],
        cell_stats["cell_auc_std"],
        cell_stats["cell_auc_min"],
        cell_stats["cell_auc_max"],
    )
    return aggregate


def _write_cell_distribution_plot(
    per_cell_metrics: list[dict[str, object]],
    output_path: Path,
) -> None:
    try:
        import matplotlib

        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        LOGGER.warning("matplotlib not installed — skipping cell distribution plot")
        return
    aucs = [m["auc"] for m in per_cell_metrics if m["auc"] is not None]
    if not aucs:
        LOGGER.warning("No valid per-cell AUCs — skipping distribution plot")
        return
    fig, ax = plt.subplots(figsize=(7, 4))
    ax.hist(aucs, bins=np.linspace(0.4, 1.0, 25), edgecolor="black", color="#4c72b0")
    ax.axvline(float(np.mean(aucs)), color="red", linestyle="--", label=f"mean={np.mean(aucs):.3f}")
    ax.axvline(float(np.median(aucs)), color="orange", linestyle="--", label=f"median={np.median(aucs):.3f}")
    ax.set_xlabel("Per-cell AUC")
    ax.set_ylabel("Count")
    ax.set_title(f"CH07 per-cell AUC distribution ({len(aucs)} cells)")
    ax.legend()
    fig.tight_layout()
    fig.savefig(output_path, dpi=140)
    plt.close(fig)
    LOGGER.info("Wrote %s", output_path)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--phage-slots-dir",
        type=Path,
        default=ARM3_SLOTS_DIR,
        help=(
            "Directory containing the phage-side slot override. Default is the CH06 Arm 3 "
            "(Moriniere per-receptor fractions) slot at .scratch/basel/feature_slots_arm3. "
            "Pass the canonical extended slots directory for the pre-Arm-3 baseline."
        ),
    )
    parser.add_argument(
        "--basel-log10-pfu-ml",
        type=float,
        default=BASEL_LOG10_PFU_ML,
    )
    parser.add_argument(
        "--num-workers",
        type=int,
        default=3,
        help="Seed-level parallelism (1 = sequential; 3 matches len(SEEDS)).",
    )
    parser.add_argument(
        "--drop-high-titer-only-positives",
        action=argparse.BooleanOptionalAction,
        default=False,
        help="Opt-in neat-only positive filter (CH10 demoted; default OFF).",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=None,
        help="Limit number of (bact_fold, phage_fold) cells evaluated — for subset iteration.",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    run_ch07_eval(
        device_type=args.device_type,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
        phage_slots_dir=args.phage_slots_dir,
        basel_log10_pfu_ml=args.basel_log10_pfu_ml,
        num_workers=args.num_workers,
        drop_high_titer_only_positives=args.drop_high_titer_only_positives,
        max_cells=args.max_cells,
    )


if __name__ == "__main__":
    main()


__all__ = [
    "run_ch07_eval",
    "ARM3_SLOTS_DIR",
]
