#!/usr/bin/env python3
"""CH05 post-hoc: isotonic-regression calibration test on Guelin predictions.

Reviewer hypothesis: Guelin's 28-30 pp mid-P miscalibration in CH05 (bacteria-axis
0.7-0.8 bin predicts 0.75 but observes 0.45) is a training-label-vs-deployment-
question mismatch signature. Training treats every score='1' row as a positive,
including "clearing events" at high titer that Gaborieau 2024 Methods explicitly
says can arise from non-productive mechanisms (lysis from without, abortive
infection). The deployment question — "will this phage productively lyse this
strain at therapeutic dose?" — is stricter than the training label. If that's
the driver, the model should have the discrimination (AUC stays) but not the
calibration, and an isotonic fit should close the gap.

Isotonic regression is non-parametric and monotonic, so it preserves ranking
(AUC unchanged) while remapping probability values to match observed frequencies.
If Guelin mid-P gap closes from 28 pp to <10 pp after isotonic, threshold story
is real. If it stays large, the model is genuinely uncertain in mid-P — threshold
isn't the driver.

Per-fold leave-one-out: for each held-out fold's predictions, fit isotonic on the
other 9 folds and apply to this fold. Otherwise the diagnostic is optimistic
(fits and evaluates on the same data).

Reads the existing CH05 predictions. No model rerun. Outputs:
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_isotonic_calibration_test.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.isotonic import IsotonicRegression
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import (
    DEFAULT_OUTPUT_DIR,
    SOURCE_BASEL,
    SOURCE_GUELIN,
    load_unified_row_frame,
)
from lyzortx.pipeline.autoresearch.sx01_eval import N_FOLDS, assign_bacteria_folds, bacteria_to_cv_group_map
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_isotonic_calibration_test.csv"
N_BINS = 10


def _reliability(preds: np.ndarray, labels: np.ndarray, label: str) -> list[dict[str, object]]:
    """Compute per-decile observed rate vs mean predicted probability."""
    edges = np.linspace(0.0, 1.0, N_BINS + 1)
    bin_idx = np.clip(np.digitize(preds, edges[1:-1], right=False), 0, N_BINS - 1)
    rows: list[dict[str, object]] = []
    for b in range(N_BINS):
        mask = bin_idx == b
        n = int(mask.sum())
        if n == 0:
            continue
        mean_pred = float(preds[mask].mean())
        obs_rate = float(labels[mask].mean())
        rows.append(
            {
                "variant": label,
                "bin_lo": round(float(edges[b]), 3),
                "bin_hi": round(float(edges[b + 1]), 3),
                "n_pairs": n,
                "mean_predicted": round(mean_pred, 4),
                "observed_rate": round(obs_rate, 4),
                "reliability_gap": round(obs_rate - mean_pred, 4),
            }
        )
    return rows


def _ece(preds: np.ndarray, labels: np.ndarray) -> float:
    """Expected-calibration-error: weighted mean of per-bin |observed - predicted| gaps."""
    edges = np.linspace(0.0, 1.0, N_BINS + 1)
    bin_idx = np.clip(np.digitize(preds, edges[1:-1], right=False), 0, N_BINS - 1)
    n_total = len(preds)
    ece = 0.0
    for b in range(N_BINS):
        mask = bin_idx == b
        n = int(mask.sum())
        if n == 0:
            continue
        mean_pred = float(preds[mask].mean())
        obs_rate = float(labels[mask].mean())
        ece += (n / n_total) * abs(obs_rate - mean_pred)
    return ece


def _log_metrics(label: str, preds: np.ndarray, labels: np.ndarray) -> dict[str, float]:
    """AUC (should be monotone-invariant), Brier, ECE — sanity checks for the transform."""
    auc = float(roc_auc_score(labels, preds)) if len(set(labels.tolist())) > 1 else float("nan")
    brier = float(brier_score_loss(labels, preds))
    ece = _ece(preds, labels)
    LOGGER.info("%s: AUC=%.4f, Brier=%.4f, ECE=%.4f, n=%d", label, auc, brier, ece, len(preds))
    return {"label": label, "auc": round(auc, 4), "brier": round(brier, 4), "ece": round(ece, 4), "n_pairs": len(preds)}


def main() -> None:
    setup_logging()
    bacteria_preds = pd.read_csv(Path(DEFAULT_OUTPUT_DIR) / "ch05_bacteria_axis_predictions.csv")
    phage_preds = pd.read_csv(Path(DEFAULT_OUTPUT_DIR) / "ch05_phage_axis_predictions.csv")

    family_map = load_unified_phage_family_map()

    bacteria_guelin = bacteria_preds[bacteria_preds["source"] == SOURCE_GUELIN].copy()
    bacteria_basel = bacteria_preds[bacteria_preds["source"] == SOURCE_BASEL].copy()
    phage_guelin = phage_preds[phage_preds["source"] == SOURCE_GUELIN].copy()
    phage_basel = phage_preds[phage_preds["source"] == SOURCE_BASEL].copy()

    rows: list[dict[str, object]] = []
    metrics: list[dict[str, object]] = []

    # Bacteria-axis: need fold_id to leave each fold out. Rebuild fold assignments
    # from the unified row frame's cv_group mapping (the predictions CSV doesn't carry cv_group).
    unified = load_unified_row_frame()
    mapping = bacteria_to_cv_group_map(unified)
    fold_assignments = assign_bacteria_folds(mapping)
    bacteria_guelin["fold_id"] = bacteria_guelin["bacteria"].map(fold_assignments)

    bacteria_calibrated = np.empty(len(bacteria_guelin), dtype=float)
    for fold_id in range(N_FOLDS):
        held_out_mask = bacteria_guelin["fold_id"].to_numpy() == fold_id
        fit_mask = ~held_out_mask
        if held_out_mask.sum() == 0 or fit_mask.sum() == 0:
            bacteria_calibrated[held_out_mask] = bacteria_guelin.loc[held_out_mask, "predicted_probability"].to_numpy()
            continue
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        iso.fit(
            bacteria_guelin.loc[fit_mask, "predicted_probability"].to_numpy(),
            bacteria_guelin.loc[fit_mask, "label_row_binary"].to_numpy().astype(float),
        )
        bacteria_calibrated[held_out_mask] = iso.predict(
            bacteria_guelin.loc[held_out_mask, "predicted_probability"].to_numpy()
        )
    bacteria_labels = bacteria_guelin["label_row_binary"].to_numpy().astype(int)
    bacteria_raw = bacteria_guelin["predicted_probability"].to_numpy().astype(float)

    metrics.append(_log_metrics("bacteria_axis_guelin_raw", bacteria_raw, bacteria_labels))
    metrics.append(_log_metrics("bacteria_axis_guelin_isotonic", bacteria_calibrated, bacteria_labels))
    rows.extend(_reliability(bacteria_raw, bacteria_labels, "bacteria_axis_guelin_raw"))
    rows.extend(_reliability(bacteria_calibrated, bacteria_labels, "bacteria_axis_guelin_isotonic"))

    # Discriminative test: apply a Guelin-fitted calibrator to BASEL predictions.
    # BASEL phages are entirely distinct from Guelin phages (disjoint panels), so fitting on
    # ALL Guelin bacteria-axis predictions and applying to BASEL is leakage-free. If BASEL
    # gap also shrinks → threshold was a unified story; if it stays → TL17-bias is the
    # BASEL-specific driver on top of Guelin's threshold mismatch.
    iso_bacteria_global = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso_bacteria_global.fit(
        bacteria_guelin["predicted_probability"].to_numpy(),
        bacteria_guelin["label_row_binary"].to_numpy().astype(float),
    )
    basel_bacteria_labels = bacteria_basel["label_row_binary"].to_numpy().astype(int)
    basel_bacteria_raw = bacteria_basel["predicted_probability"].to_numpy().astype(float)
    basel_bacteria_iso = iso_bacteria_global.predict(basel_bacteria_raw)
    metrics.append(_log_metrics("bacteria_axis_basel_raw", basel_bacteria_raw, basel_bacteria_labels))
    metrics.append(
        _log_metrics("bacteria_axis_basel_guelin_isotonic_applied", basel_bacteria_iso, basel_bacteria_labels)
    )
    rows.extend(_reliability(basel_bacteria_raw, basel_bacteria_labels, "bacteria_axis_basel_raw"))
    rows.extend(_reliability(basel_bacteria_iso, basel_bacteria_labels, "bacteria_axis_basel_guelin_isotonic_applied"))

    # Phage-axis: held-out phages. Leave each phage out for its calibration fit.
    phage_guelin["family"] = phage_guelin["phage"].map(family_map).fillna("UNKNOWN")
    from lyzortx.pipeline.autoresearch.ch05_eval import assign_phage_folds

    phage_fold_assignments = assign_phage_folds(
        sorted(phage_guelin["phage"].unique().tolist()),
        family_map,
    )
    phage_guelin["fold_id"] = phage_guelin["phage"].map(phage_fold_assignments)

    phage_calibrated = np.empty(len(phage_guelin), dtype=float)
    for fold_id in sorted(phage_guelin["fold_id"].dropna().unique().astype(int)):
        held_out_mask = phage_guelin["fold_id"].to_numpy() == fold_id
        fit_mask = ~held_out_mask
        if held_out_mask.sum() == 0 or fit_mask.sum() == 0:
            phage_calibrated[held_out_mask] = phage_guelin.loc[held_out_mask, "predicted_probability"].to_numpy()
            continue
        iso = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
        iso.fit(
            phage_guelin.loc[fit_mask, "predicted_probability"].to_numpy(),
            phage_guelin.loc[fit_mask, "label_row_binary"].to_numpy().astype(float),
        )
        phage_calibrated[held_out_mask] = iso.predict(
            phage_guelin.loc[held_out_mask, "predicted_probability"].to_numpy()
        )
    phage_labels = phage_guelin["label_row_binary"].to_numpy().astype(int)
    phage_raw = phage_guelin["predicted_probability"].to_numpy().astype(float)

    metrics.append(_log_metrics("phage_axis_guelin_raw", phage_raw, phage_labels))
    metrics.append(_log_metrics("phage_axis_guelin_isotonic", phage_calibrated, phage_labels))
    rows.extend(_reliability(phage_raw, phage_labels, "phage_axis_guelin_raw"))
    rows.extend(_reliability(phage_calibrated, phage_labels, "phage_axis_guelin_isotonic"))

    # Discriminative test (phage-axis): apply Guelin-fitted calibrator to BASEL predictions.
    # Guelin (96) and BASEL (52) phage panels are disjoint, so fitting on all Guelin predictions
    # and applying to BASEL has no leakage.
    iso_phage_global = IsotonicRegression(out_of_bounds="clip", y_min=0.0, y_max=1.0)
    iso_phage_global.fit(
        phage_guelin["predicted_probability"].to_numpy(),
        phage_guelin["label_row_binary"].to_numpy().astype(float),
    )
    basel_phage_labels = phage_basel["label_row_binary"].to_numpy().astype(int)
    basel_phage_raw = phage_basel["predicted_probability"].to_numpy().astype(float)
    basel_phage_iso = iso_phage_global.predict(basel_phage_raw)
    metrics.append(_log_metrics("phage_axis_basel_raw", basel_phage_raw, basel_phage_labels))
    metrics.append(_log_metrics("phage_axis_basel_guelin_isotonic_applied", basel_phage_iso, basel_phage_labels))
    rows.extend(_reliability(basel_phage_raw, basel_phage_labels, "phage_axis_basel_raw"))
    rows.extend(_reliability(basel_phage_iso, basel_phage_labels, "phage_axis_basel_guelin_isotonic_applied"))

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    metrics_df = pd.DataFrame(metrics)
    metrics_df.to_csv(Path(DEFAULT_OUTPUT_DIR) / "ch05_isotonic_calibration_test_metrics.csv", index=False)
    LOGGER.info("Wrote %s (%d rows) + _metrics.csv (%d rows)", OUTPUT_CSV, len(df), len(metrics_df))

    print()
    print("=== Isotonic calibration test ===")
    print("Guelin: leave-one-fold-out (fold-aware, no leakage)")
    print("BASEL: Guelin-fitted calibrator applied to BASEL predictions (disjoint panels, no leakage)")
    print()
    print(metrics_df.to_string(index=False))
    print()
    for axis in ("bacteria_axis", "phage_axis"):
        for source in ("guelin", "basel"):
            prefix = f"{axis}_{source}"
            raw = df[df["variant"] == f"{prefix}_raw"]["reliability_gap"].abs()
            iso_label_suffix = "isotonic" if source == "guelin" else "guelin_isotonic_applied"
            iso = df[df["variant"] == f"{prefix}_{iso_label_suffix}"]["reliability_gap"].abs()
            if raw.empty or iso.empty:
                continue
            print(
                f"{axis} {source}: max |gap| raw={raw.max():.3f}, isotonic={iso.max():.3f}, "
                f"closure={raw.max() - iso.max():+.3f}"
            )


if __name__ == "__main__":
    main()
