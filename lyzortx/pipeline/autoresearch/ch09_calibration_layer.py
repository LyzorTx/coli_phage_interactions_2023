#!/usr/bin/env python3
"""CH09: Post-hoc calibration layer (Arms 1 + 2).

Productionises the CH05 isotonic diagnostic into a deployable calibration layer:
  - ARM 1: fit a single isotonic regressor on ALL Guelin CHISEL training-fold
    predictions (bacteria-axis, pooled), persist as a pipeline artifact under
    `lyzortx/generated_outputs/ch09_calibration_layer/ch09_calibrator.pkl`. Also
    report the leave-one-fold-out (LOOF) ECE that would be obtained at deployment
    time, so the reported Guelin ECE is unbiased by in-sample fit.
  - ARM 2: apply the fitted Guelin calibrator to BASEL predictions (bacteria-axis
    and phage-axis) and report ECE closure. The CH05 diagnostic pre-registered
    ~34-37% closure — CH09 confirms it on the production layer.

ARM 3 (label-threshold sensitivity — retraining CH04 with ambiguous "clearing-
only-at-high-titer" pairs excluded) is deferred to a follow-up commit on this
branch. It requires a CH04 rerun (~25 min on the parallel loop from CH06
pre-flight) and a label-policy change.

Outputs:
  - ch09_calibrator.pkl          (joblib-serialised IsotonicRegression)
  - ch09_calibration_report.json (raw vs calibrated ECE/Brier/AUC, axis x source)
  - ch09_reliability_tables.csv  (per-decile reliability, raw vs calibrated,
                                   all four axis x source subsets)

Reads existing CH05 prediction artifacts under
`lyzortx/generated_outputs/ch05_unified_kfold/` — no model retraining.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Optional

import joblib
import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import (
    SOURCE_BASEL,
    SOURCE_GUELIN,
    assign_phage_folds,
    load_unified_row_frame,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

DEFAULT_CH05_DIR = Path("lyzortx/generated_outputs/ch05_unified_kfold")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch09_calibration_layer")
CALIBRATOR_FILENAME = "ch09_calibrator.pkl"
REPORT_FILENAME = "ch09_calibration_report.json"
RELIABILITY_FILENAME = "ch09_reliability_tables.csv"

ECE_BINS = 10


def expected_calibration_error(
    labels: np.ndarray,
    predictions: np.ndarray,
    *,
    n_bins: int = ECE_BINS,
) -> float:
    """Weighted mean of per-decile |observed − predicted|. Matches CH05 reliability tables."""
    if len(labels) == 0:
        return float("nan")
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_idx = np.clip(np.digitize(predictions, edges) - 1, 0, n_bins - 1)
    total = 0.0
    total_n = 0
    for b in range(n_bins):
        mask = bin_idx == b
        bin_count = int(mask.sum())
        if bin_count == 0:
            continue
        gap = abs(labels[mask].mean() - predictions[mask].mean())
        total += bin_count * gap
        total_n += bin_count
    return float(total / total_n) if total_n else float("nan")


def reliability_table(
    preds: np.ndarray,
    labels: np.ndarray,
    *,
    variant: str,
    n_bins: int = ECE_BINS,
) -> list[dict[str, object]]:
    """Per-decile observed rate vs mean predicted probability (reliability diagram row form)."""
    edges = np.linspace(0.0, 1.0, n_bins + 1)
    bin_idx = np.clip(np.digitize(preds, edges) - 1, 0, n_bins - 1)
    rows: list[dict[str, object]] = []
    for b in range(n_bins):
        mask = bin_idx == b
        n = int(mask.sum())
        if n == 0:
            continue
        rows.append(
            {
                "variant": variant,
                "bin_lo": float(edges[b]),
                "bin_hi": float(edges[b + 1]),
                "n_pairs": n,
                "mean_predicted": float(preds[mask].mean()),
                "observed_rate": float(labels[mask].mean()),
                "reliability_gap": float(labels[mask].mean() - preds[mask].mean()),
            }
        )
    return rows


def subset_metrics(preds: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    """AUC (invariant under monotone isotonic), Brier, ECE for one subset."""
    if len(labels) == 0:
        return {"auc": float("nan"), "brier": float("nan"), "ece": float("nan"), "n_pairs": 0}
    has_both = len(np.unique(labels)) > 1
    return {
        "auc": float(roc_auc_score(labels, preds)) if has_both else float("nan"),
        "brier": float(brier_score_loss(labels, preds)),
        "ece": float(expected_calibration_error(labels, preds)),
        "n_pairs": int(len(labels)),
    }


def fit_production_calibrator(
    guelin_preds: np.ndarray,
    guelin_labels: np.ndarray,
) -> IsotonicRegression:
    """Fit a SINGLE isotonic regressor on ALL Guelin training-fold predictions.

    This is the DEPLOYMENT artifact — it uses every Guelin data point. For
    UNBIASED evaluation (LOOF), see `apply_loof_calibration`.
    """
    iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso.fit(guelin_preds.astype(float), guelin_labels.astype(float))
    LOGGER.info(
        "Production isotonic calibrator fitted on %d Guelin pairs; fitted monotone map has %d change points.",
        len(guelin_preds),
        len(iso.X_thresholds_),
    )
    return iso


def apply_loof_calibration(
    guelin_preds: np.ndarray,
    guelin_labels: np.ndarray,
    fold_ids: np.ndarray,
) -> np.ndarray:
    """Leave-one-fold-out: for each fold, fit isotonic on OTHER folds, predict on held-out fold.

    Gives an unbiased deployment-time ECE estimate (as if each fold were the true
    cold-start test). Used for the ECE < 0.02 acceptance check on Guelin.
    """
    calibrated = np.empty(len(guelin_preds), dtype=float)
    for fold_id in sorted(set(fold_ids.tolist())):
        if fold_id < 0:  # unassigned fold — leave raw
            held_mask = fold_ids == fold_id
            calibrated[held_mask] = guelin_preds[held_mask]
            continue
        held_mask = fold_ids == fold_id
        fit_mask = ~held_mask
        if held_mask.sum() == 0 or fit_mask.sum() == 0 or len(np.unique(guelin_labels[fit_mask])) < 2:
            calibrated[held_mask] = guelin_preds[held_mask]
            continue
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        iso.fit(guelin_preds[fit_mask].astype(float), guelin_labels[fit_mask].astype(float))
        calibrated[held_mask] = iso.predict(guelin_preds[held_mask].astype(float))
    return calibrated


def attach_fold_ids(
    predictions: pd.DataFrame,
    *,
    axis: str,
    unified_row_frame: pd.DataFrame,
    phage_family: dict[str, str],
) -> pd.DataFrame:
    """Attach `fold_id` column to a CH05 predictions DataFrame. `axis` ∈ {bacteria, phage}."""
    out = predictions.copy()
    if axis == "bacteria":
        mapping = bacteria_to_cv_group_map(unified_row_frame)
        fold_map = assign_bacteria_folds(mapping)
        out["fold_id"] = out["bacteria"].map(fold_map)
    elif axis == "phage":
        phages = sorted(predictions["phage"].unique())
        fold_map = assign_phage_folds(phages, phage_family)
        out["fold_id"] = out["phage"].map(fold_map)
    else:
        raise ValueError(f"Unknown axis: {axis}")
    missing = out["fold_id"].isna().sum()
    if missing:
        LOGGER.warning("axis=%s: %d prediction rows missing fold_id (marking as -1)", axis, missing)
        out["fold_id"] = out["fold_id"].fillna(-1).astype(int)
    else:
        out["fold_id"] = out["fold_id"].astype(int)
    return out


def run_ch09(
    *,
    ch05_dir: Path,
    output_dir: Path,
) -> dict[str, object]:
    """CH09 Arms 1 + 2 driver.

    Loads CH05 predictions, fits LOOF and production calibrators on Guelin,
    applies the production calibrator to BASEL, writes artifacts, returns
    a summary dict.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("CH09 calibration layer starting")

    bacteria_preds = pd.read_csv(ch05_dir / "ch05_bacteria_axis_predictions.csv")
    phage_preds = pd.read_csv(ch05_dir / "ch05_phage_axis_predictions.csv")
    unified_row_frame = load_unified_row_frame()
    phage_family = load_unified_phage_family_map()

    bacteria_preds = attach_fold_ids(
        bacteria_preds, axis="bacteria", unified_row_frame=unified_row_frame, phage_family=phage_family
    )
    phage_preds = attach_fold_ids(
        phage_preds, axis="phage", unified_row_frame=unified_row_frame, phage_family=phage_family
    )

    bacteria_guelin = bacteria_preds[bacteria_preds["source"] == SOURCE_GUELIN].reset_index(drop=True)
    bacteria_basel = bacteria_preds[bacteria_preds["source"] == SOURCE_BASEL].reset_index(drop=True)
    phage_guelin = phage_preds[phage_preds["source"] == SOURCE_GUELIN].reset_index(drop=True)
    phage_basel = phage_preds[phage_preds["source"] == SOURCE_BASEL].reset_index(drop=True)

    # ARM 1 — Production calibrator + LOOF evaluation on bacteria-axis Guelin.
    # The production calibrator is fit on ALL Guelin bacteria-axis predictions (the
    # larger, more diverse pool). LOOF ECE is reported for an unbiased deployment
    # ECE estimate (the production calibrator itself is NOT held out when persisted —
    # it's the full-data fit that a deployment pipeline would use).
    bacteria_guelin_raw = bacteria_guelin["predicted_probability"].to_numpy(dtype=float)
    bacteria_guelin_labels = bacteria_guelin["label_row_binary"].to_numpy(dtype=int)
    bacteria_guelin_fold = bacteria_guelin["fold_id"].to_numpy(dtype=int)
    bacteria_guelin_loof = apply_loof_calibration(bacteria_guelin_raw, bacteria_guelin_labels, bacteria_guelin_fold)
    production_calibrator = fit_production_calibrator(bacteria_guelin_raw, bacteria_guelin_labels)
    calibrator_path = output_dir / CALIBRATOR_FILENAME
    joblib.dump(production_calibrator, calibrator_path)
    LOGGER.info("Persisted production calibrator to %s", calibrator_path)

    # ARM 1 — also report LOOF on phage-axis Guelin (different fold definition).
    phage_guelin_raw = phage_guelin["predicted_probability"].to_numpy(dtype=float)
    phage_guelin_labels = phage_guelin["label_row_binary"].to_numpy(dtype=int)
    phage_guelin_fold = phage_guelin["fold_id"].to_numpy(dtype=int)
    phage_guelin_loof = apply_loof_calibration(phage_guelin_raw, phage_guelin_labels, phage_guelin_fold)

    # ARM 2 — apply the production calibrator to BASEL (cross-panel transfer).
    bacteria_basel_raw = bacteria_basel["predicted_probability"].to_numpy(dtype=float)
    bacteria_basel_labels = bacteria_basel["label_row_binary"].to_numpy(dtype=int)
    bacteria_basel_transfer = production_calibrator.predict(bacteria_basel_raw)

    phage_basel_raw = phage_basel["predicted_probability"].to_numpy(dtype=float)
    phage_basel_labels = phage_basel["label_row_binary"].to_numpy(dtype=int)
    phage_basel_transfer = production_calibrator.predict(phage_basel_raw)

    # Metrics.
    axis_source_metrics: dict[str, dict[str, dict[str, float]]] = {
        "bacteria_axis": {
            "guelin_raw": subset_metrics(bacteria_guelin_raw, bacteria_guelin_labels),
            "guelin_loof_calibrated": subset_metrics(bacteria_guelin_loof, bacteria_guelin_labels),
            "basel_raw": subset_metrics(bacteria_basel_raw, bacteria_basel_labels),
            "basel_guelin_calibrator_applied": subset_metrics(bacteria_basel_transfer, bacteria_basel_labels),
        },
        "phage_axis": {
            "guelin_raw": subset_metrics(phage_guelin_raw, phage_guelin_labels),
            "guelin_loof_calibrated": subset_metrics(phage_guelin_loof, phage_guelin_labels),
            "basel_raw": subset_metrics(phage_basel_raw, phage_basel_labels),
            "basel_guelin_calibrator_applied": subset_metrics(phage_basel_transfer, phage_basel_labels),
        },
    }

    # Transfer closure = (raw_ECE - calibrated_ECE) / raw_ECE. Guelin uses LOOF (unbiased);
    # BASEL uses production calibrator applied directly (cross-panel transfer).
    def closure(raw_ece: float, cal_ece: float) -> Optional[float]:
        if raw_ece is None or cal_ece is None or np.isnan(raw_ece) or np.isnan(cal_ece) or raw_ece <= 0:
            return None
        return (raw_ece - cal_ece) / raw_ece

    summary = {
        "task_id": "CH09",
        "scope": "Arm 1 (production calibrator) + Arm 2 (cross-panel transfer); Arm 3 (label-threshold sensitivity) deferred.",
        "calibrator_artifact": str(calibrator_path),
        "calibrator_fit_n": int(len(bacteria_guelin_raw)),
        "calibrator_fit_source": "Guelin bacteria-axis training-fold predictions (pooled)",
        "axis_source_metrics": axis_source_metrics,
        "arm1_guelin_bacteria_axis_loof_ece": axis_source_metrics["bacteria_axis"]["guelin_loof_calibrated"]["ece"],
        "arm1_guelin_phage_axis_loof_ece": axis_source_metrics["phage_axis"]["guelin_loof_calibrated"]["ece"],
        "arm1_acceptance_bacteria_axis_pass": bool(
            axis_source_metrics["bacteria_axis"]["guelin_loof_calibrated"]["ece"] < 0.02
        ),
        "arm1_acceptance_phage_axis_pass": bool(
            axis_source_metrics["phage_axis"]["guelin_loof_calibrated"]["ece"] < 0.02
        ),
        "arm2_basel_bacteria_axis_closure": closure(
            axis_source_metrics["bacteria_axis"]["basel_raw"]["ece"],
            axis_source_metrics["bacteria_axis"]["basel_guelin_calibrator_applied"]["ece"],
        ),
        "arm2_basel_phage_axis_closure": closure(
            axis_source_metrics["phage_axis"]["basel_raw"]["ece"],
            axis_source_metrics["phage_axis"]["basel_guelin_calibrator_applied"]["ece"],
        ),
    }

    with open(output_dir / REPORT_FILENAME, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)

    # Reliability tables (per-decile, for both axes × (guelin raw, guelin calibrated,
    # basel raw, basel transferred)).
    reliability_rows: list[dict[str, object]] = []
    for label, preds, labels in (
        ("bacteria_axis_guelin_raw", bacteria_guelin_raw, bacteria_guelin_labels),
        ("bacteria_axis_guelin_loof_calibrated", bacteria_guelin_loof, bacteria_guelin_labels),
        ("bacteria_axis_basel_raw", bacteria_basel_raw, bacteria_basel_labels),
        ("bacteria_axis_basel_guelin_calibrator_applied", bacteria_basel_transfer, bacteria_basel_labels),
        ("phage_axis_guelin_raw", phage_guelin_raw, phage_guelin_labels),
        ("phage_axis_guelin_loof_calibrated", phage_guelin_loof, phage_guelin_labels),
        ("phage_axis_basel_raw", phage_basel_raw, phage_basel_labels),
        ("phage_axis_basel_guelin_calibrator_applied", phage_basel_transfer, phage_basel_labels),
    ):
        reliability_rows.extend(reliability_table(preds, labels, variant=label))
    pd.DataFrame(reliability_rows).to_csv(output_dir / RELIABILITY_FILENAME, index=False)

    LOGGER.info(
        "CH09 done. Guelin LOOF ECE: bacteria=%.4f, phage=%.4f (target <0.02). "
        "BASEL transfer closure: bacteria=%s, phage=%s.",
        axis_source_metrics["bacteria_axis"]["guelin_loof_calibrated"]["ece"],
        axis_source_metrics["phage_axis"]["guelin_loof_calibrated"]["ece"],
        f"{summary['arm2_basel_bacteria_axis_closure']:.1%}" if summary["arm2_basel_bacteria_axis_closure"] else "n/a",
        f"{summary['arm2_basel_phage_axis_closure']:.1%}" if summary["arm2_basel_phage_axis_closure"] else "n/a",
    )
    return summary


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ch05-dir", type=Path, default=DEFAULT_CH05_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    setup_logging()
    args = parse_args(argv)
    run_ch09(ch05_dir=args.ch05_dir, output_dir=args.output_dir)


if __name__ == "__main__":
    main()
