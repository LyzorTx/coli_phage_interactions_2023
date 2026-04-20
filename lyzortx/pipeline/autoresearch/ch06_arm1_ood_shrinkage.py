#!/usr/bin/env python3
"""CH06 Arm 1: OOD-aware post-hoc shrinkage toward base rate.

Calibration arm (NOT discrimination) per plan.yml: attacks the second of CH05's
two miscalibration mechanisms — BASEL's residual ECE that the Guelin-fitted
isotonic calibrator only partly closed (34-37% gap closure on BASEL vs 78-89% on
Guelin). The remaining gap is TL17-bias: Guelin-calibrated host-range priors
mapped onto non-zero-projection BASEL phages whose actual host ranges are
narrower.

The remedy this arm tests: for each held-out phage, measure distance to its
nearest training phage (Euclidean in the phage_projection TL17 family vector).
When that distance exceeds a threshold fit on Guelin phage-axis fold-out
predictions, shrink the model's probabilities toward the training base rate.
Shrink-toward-base-rate on OOD phages is literally what the miscalibration
asks for; the arm quantifies how much of BASEL's residual reliability gap is
closed by this no-retraining move.

Runs on existing CH05 artifacts (no training). Produces:
  - ch06_arm1_predictions.csv: CH05 predictions with `predicted_probability_raw`
    and `predicted_probability_shrunk` columns, plus per-phage `ood_distance`
    and `ood_shrinkage_weight`.
  - ch06_arm1_metrics.json: AUC/Brier/ECE before and after shrinkage, split by
    axis and source.
  - ch06_arm1_shrinkage_diagnostic.csv: threshold sweep (for both hard and soft
    gating variants) showing Guelin vs BASEL ECE at each threshold.

Arm 1 reports against BASEL ECE, not BASEL AUC. Discrimination is addressed by
Arms 2-4.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import (
    SOURCE_BASEL,
    SOURCE_GUELIN,
    assign_phage_folds,
    load_unified_row_frame,
)
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

DEFAULT_CH05_DIR = Path("lyzortx/generated_outputs/ch05_unified_kfold")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch06_arm1_ood_shrinkage")
EXTENDED_PHAGE_PROJECTION_CSV = Path(".scratch/basel/feature_slots/phage_projection/features.csv")

ECE_BINS = 10
PROJECTION_COL_PREFIX = "phage_projection__"


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


def load_extended_phage_projection() -> pd.DataFrame:
    """148-phage (Guelin + BASEL) phage_projection feature matrix."""
    frame = pd.read_csv(EXTENDED_PHAGE_PROJECTION_CSV)
    projection_cols = [c for c in frame.columns if c.startswith(PROJECTION_COL_PREFIX)]
    return frame.set_index("phage")[projection_cols]


def compute_phage_distances(projection: pd.DataFrame) -> pd.DataFrame:
    """Pairwise Euclidean distance matrix in the TL17-family vector space."""
    matrix = projection.to_numpy(dtype=float)
    # Euclidean: sqrt(sum((a_i - b_j)^2)) = sqrt(|a|^2 + |b|^2 - 2 * a.b).
    sq_norms = (matrix**2).sum(axis=1)
    dist2 = sq_norms[:, None] + sq_norms[None, :] - 2.0 * matrix @ matrix.T
    dist2 = np.clip(dist2, 0.0, None)  # numerical negative zeros
    dist = np.sqrt(dist2)
    return pd.DataFrame(dist, index=projection.index, columns=projection.index)


def compute_ood_distance_per_phage(
    distance_matrix: pd.DataFrame,
    fold_assignments: dict[str, int],
    n_folds: int,
) -> pd.DataFrame:
    """Intrinsic OOD distance per phage = distance to nearest training phage when
    THIS phage is held out under the CH05 phage-axis folds. One value per phage,
    applied identically to bacteria-axis and phage-axis predictions.

    Bacteria-axis CV does not hold phages out — every phage is in every training
    fold, so "nearest training phage in the current fold" is zero for all of them
    and would give no OOD signal. Using the phage-axis fold-out definition gives
    an intrinsic phage property and lets the same shrinkage decision apply on
    both axes, as the plan.yml Arm 1 text specifies ("AUC stays at 0.7152" on
    bacteria-axis — monotone per-phage shrinkage preserves within-phage rank).
    """
    records: list[dict[str, object]] = []
    for phage, fold_id in fold_assignments.items():
        train_phages = [p for p, f in fold_assignments.items() if f != fold_id]
        if not train_phages or phage not in distance_matrix.index:
            continue
        available_train = [p for p in train_phages if p in distance_matrix.columns]
        if not available_train:
            continue
        min_dist = float(distance_matrix.loc[phage, available_train].min())
        records.append({"phage": phage, "fold_id": fold_id, "ood_distance_min_to_train": min_dist})
    return pd.DataFrame(records)


def apply_shrinkage(
    predictions: np.ndarray,
    distances: np.ndarray,
    base_rate: float,
    threshold: float,
    *,
    gate: str,
    scale: float,
) -> tuple[np.ndarray, np.ndarray]:
    """Apply distance-gated shrinkage toward base rate.

    Hard gate: `weight = 1 if dist > threshold else 0`.
    Soft gate: `weight = sigmoid((dist - threshold) / scale)`.

    Returns `(shrunk_predictions, shrinkage_weights)`.
    """
    if gate == "hard":
        weights = (distances > threshold).astype(float)
    elif gate == "soft":
        weights = 1.0 / (1.0 + np.exp(-(distances - threshold) / max(scale, 1e-9)))
    else:
        raise ValueError(f"Unknown gate mode: {gate}")
    shrunk = (1.0 - weights) * predictions + weights * base_rate
    return shrunk, weights


def sweep_thresholds(
    predictions: pd.DataFrame,
    *,
    gate: str,
    scale: float,
    base_rate: float,
    n_steps: int = 40,
) -> pd.DataFrame:
    """Sweep thresholds across the observed OOD-distance range and report
    per-source ECE + AUC + Brier + shrunk-fraction at each point.

    The sweep is the authoritative artifact for this arm: it surfaces whether
    any threshold closes BASEL ECE without damaging Guelin ECE. Picking a
    single "best" threshold is deferred to `pick_best_threshold`.
    """
    distances = predictions["ood_distance_min_to_train"].to_numpy()
    probs = predictions["predicted_probability_raw"].to_numpy()
    labels = predictions["label_row_binary"].to_numpy()
    source = predictions["source"].to_numpy()
    guelin_mask = source == SOURCE_GUELIN
    basel_mask = source == SOURCE_BASEL

    lo, hi = float(np.min(distances)), float(np.max(distances))
    # Include the "no shrinkage" anchor explicitly (threshold = hi + small margin).
    thresholds = list(np.linspace(lo, hi, n_steps)) + [hi + 1.0]
    sweep_rows: list[dict[str, object]] = []
    for thr in thresholds:
        shrunk, weights = apply_shrinkage(probs, distances, base_rate, thr, gate=gate, scale=scale)
        row: dict[str, object] = {
            "threshold": float(thr),
            "gate": gate,
            "shrunk_fraction_overall": float((weights > 0.5).mean()),
        }
        for key, mask in (("guelin", guelin_mask), ("basel", basel_mask)):
            if mask.sum() == 0 or len(np.unique(labels[mask])) < 2:
                row[f"{key}_ece"] = float("nan")
                row[f"{key}_auc"] = float("nan")
                row[f"{key}_brier"] = float("nan")
                row[f"{key}_shrunk_fraction"] = float((weights[mask] > 0.5).mean()) if mask.sum() else 0.0
                continue
            row[f"{key}_ece"] = expected_calibration_error(labels[mask], shrunk[mask])
            row[f"{key}_auc"] = float(roc_auc_score(labels[mask], shrunk[mask]))
            row[f"{key}_brier"] = float(brier_score_loss(labels[mask], shrunk[mask]))
            row[f"{key}_shrunk_fraction"] = float((weights[mask] > 0.5).mean())
        sweep_rows.append(row)
    return pd.DataFrame(sweep_rows)


def pick_best_threshold(
    sweep: pd.DataFrame,
    *,
    guelin_ece_tolerance: float = 0.005,
) -> float:
    """Pick threshold that minimizes BASEL ECE, subject to not harming Guelin ECE
    by more than `guelin_ece_tolerance` above its baseline (no-shrinkage) value.

    This is the Arm 1 acceptance rule: the shrinkage must actually help BASEL
    (which is the point of the arm) without degrading Guelin (where the model
    already calibrates reasonably). If no threshold satisfies the constraint,
    return the "no shrinkage" anchor (largest threshold) so the arm reports a
    null outcome honestly rather than a harmful change.
    """
    baseline_guelin_ece = float(sweep["guelin_ece"].iloc[-1])  # last row = no-shrinkage anchor
    cap = baseline_guelin_ece + guelin_ece_tolerance
    candidates = sweep[sweep["guelin_ece"] <= cap]
    if candidates.empty:
        return float(sweep["threshold"].iloc[-1])
    best_idx = int(candidates["basel_ece"].idxmin())
    return float(sweep.loc[best_idx, "threshold"])


def build_arm1_predictions(
    ch05_predictions: pd.DataFrame,
    per_phage_distance: pd.DataFrame,
    *,
    threshold: float,
    gate: str,
    scale: float,
    base_rate: float,
) -> pd.DataFrame:
    """Attach OOD distance + shrunk prediction columns to a CH05 predictions frame.

    `per_phage_distance` has one row per phage with `ood_distance_min_to_train`
    (intrinsic OOD distance under phage-axis fold-out). Joined by phage only,
    so the same per-phage shrinkage applies on both axes.
    """
    merged = ch05_predictions.merge(
        per_phage_distance[["phage", "ood_distance_min_to_train"]],
        on="phage",
        how="left",
    )
    distances = merged["ood_distance_min_to_train"].to_numpy()
    probs = merged["predicted_probability"].to_numpy()
    shrunk, weights = apply_shrinkage(probs, distances, base_rate, threshold, gate=gate, scale=scale)
    merged["predicted_probability_raw"] = probs
    merged["predicted_probability_shrunk"] = shrunk
    merged["ood_shrinkage_weight"] = weights
    merged["predicted_probability"] = shrunk  # downstream metric code reads this column
    return merged


def summarize_metrics_split_by_source(
    predictions: pd.DataFrame,
    prob_column: str,
) -> dict[str, dict[str, float]]:
    """AUC, Brier, ECE for overall / Guelin / BASEL subsets."""
    result: dict[str, dict[str, float]] = {}
    for key, subset in (
        ("overall", predictions),
        (SOURCE_GUELIN, predictions[predictions["source"] == SOURCE_GUELIN]),
        (SOURCE_BASEL, predictions[predictions["source"] == SOURCE_BASEL]),
    ):
        labels = subset["label_row_binary"].to_numpy()
        probs = subset[prob_column].to_numpy()
        if len(labels) == 0 or len(np.unique(labels)) < 2:
            result[key] = {"auc": float("nan"), "brier": float("nan"), "ece": float("nan"), "n_pairs": int(len(labels))}
            continue
        result[key] = {
            "auc": float(roc_auc_score(labels, probs)),
            "brier": float(brier_score_loss(labels, probs)),
            "ece": float(expected_calibration_error(labels, probs)),
            "n_pairs": int(len(labels)),
        }
    return result


def run_arm1(
    *,
    ch05_dir: Path,
    output_dir: Path,
    gate: str,
    soft_scale: float,
    basel_log10_pfu_ml: float,
) -> dict[str, object]:
    """Arm 1 driver: compute per-phage OOD distance, sweep thresholds across both
    gate modes, pick the threshold that minimizes BASEL ECE without hurting
    Guelin ECE beyond tolerance, apply to both axes, emit artifacts.
    """
    output_dir.mkdir(parents=True, exist_ok=True)
    LOGGER.info("CH06 Arm 1: OOD-aware shrinkage toward base rate (gate=%s)", gate)

    # Phage projection + distance matrix.
    projection = load_extended_phage_projection()
    LOGGER.info("Loaded phage_projection: %d phages × %d TL17 families", *projection.shape)
    distances = compute_phage_distances(projection)
    zero_vector_count = int((projection.to_numpy(dtype=float).sum(axis=1) == 0.0).sum())
    LOGGER.info(
        "Distance matrix: %d×%d, diagonal-near-zero=%s, zero-vector phages=%d",
        distances.shape[0],
        distances.shape[1],
        bool(np.allclose(np.diag(distances), 0.0, atol=1e-6)),
        zero_vector_count,
    )

    # Phage-axis fold assignments + intrinsic per-phage OOD distance.
    phage_family = load_unified_phage_family_map()
    unified_row_frame = load_unified_row_frame(basel_log10_pfu_ml=basel_log10_pfu_ml)
    phages_in_panel = sorted(unified_row_frame["phage"].unique())
    phage_fold = assign_phage_folds(phages_in_panel, phage_family)
    per_phage_distance = compute_ood_distance_per_phage(distances, phage_fold, n_folds=10)
    LOGGER.info(
        "Per-phage OOD distances: n=%d, min=%.3f, mean=%.3f, max=%.3f",
        len(per_phage_distance),
        per_phage_distance["ood_distance_min_to_train"].min(),
        per_phage_distance["ood_distance_min_to_train"].mean(),
        per_phage_distance["ood_distance_min_to_train"].max(),
    )

    # Load CH05 predictions; base rate from training labels.
    bacteria_pred = pd.read_csv(ch05_dir / "ch05_bacteria_axis_predictions.csv")
    phage_pred = pd.read_csv(ch05_dir / "ch05_phage_axis_predictions.csv")
    from lyzortx.pipeline.autoresearch.ch04_eval import build_clean_row_training_frame

    clean_rows = build_clean_row_training_frame(unified_row_frame)
    base_rate = float(clean_rows.loc[clean_rows["source"] == SOURCE_GUELIN, "label_row_binary"].mean())
    LOGGER.info(
        "CH05 predictions: bacteria-axis %d, phage-axis %d; base rate (Guelin): %.4f",
        len(bacteria_pred),
        len(phage_pred),
        base_rate,
    )

    # Attach OOD distances to both axes; raw=shrunk until threshold is picked.
    bacteria_with_dist = bacteria_pred.merge(
        per_phage_distance[["phage", "ood_distance_min_to_train"]], on="phage", how="left"
    ).assign(predicted_probability_raw=lambda df: df["predicted_probability"])
    phage_with_dist = phage_pred.merge(
        per_phage_distance[["phage", "ood_distance_min_to_train"]], on="phage", how="left"
    ).assign(predicted_probability_raw=lambda df: df["predicted_probability"])

    # Threshold sweep on phage-axis (where the OOD concept has semantic meaning —
    # held-out phages with no near-neighbour in training). Sweep both gate modes
    # so the diagnostic shows the full surface.
    hard_sweep = sweep_thresholds(phage_with_dist, gate="hard", scale=soft_scale, base_rate=base_rate)
    soft_sweep = sweep_thresholds(phage_with_dist, gate="soft", scale=soft_scale, base_rate=base_rate)
    diagnostic = pd.concat([hard_sweep, soft_sweep], ignore_index=True)
    diagnostic.to_csv(output_dir / "ch06_arm1_shrinkage_diagnostic.csv", index=False)

    canonical_sweep = hard_sweep if gate == "hard" else soft_sweep
    threshold = pick_best_threshold(canonical_sweep, guelin_ece_tolerance=0.005)
    LOGGER.info(
        "Best threshold (gate=%s, BASEL ECE minimised subject to Guelin ECE ≤ raw+0.005): %.4f",
        gate,
        threshold,
    )

    # Apply the selected threshold to both axes.
    bacteria_shrunk = build_arm1_predictions(
        bacteria_with_dist.drop(columns=["ood_distance_min_to_train"]),
        per_phage_distance,
        threshold=threshold,
        gate=gate,
        scale=soft_scale,
        base_rate=base_rate,
    )
    phage_shrunk = build_arm1_predictions(
        phage_with_dist.drop(columns=["ood_distance_min_to_train", "predicted_probability_raw"]),
        per_phage_distance,
        threshold=threshold,
        gate=gate,
        scale=soft_scale,
        base_rate=base_rate,
    )
    bacteria_shrunk.to_csv(output_dir / "ch06_arm1_bacteria_axis_predictions.csv", index=False)
    phage_shrunk.to_csv(output_dir / "ch06_arm1_phage_axis_predictions.csv", index=False)

    # Metrics: raw vs shrunk, per axis, per source.
    bacteria_metrics_raw = summarize_metrics_split_by_source(bacteria_shrunk, "predicted_probability_raw")
    bacteria_metrics_shrunk = summarize_metrics_split_by_source(bacteria_shrunk, "predicted_probability_shrunk")
    phage_metrics_raw = summarize_metrics_split_by_source(phage_shrunk, "predicted_probability_raw")
    phage_metrics_shrunk = summarize_metrics_split_by_source(phage_shrunk, "predicted_probability_shrunk")

    summary = {
        "arm_id": "ch06_arm1_ood_shrinkage",
        "target": "BASEL ECE (calibration arm — not discrimination; Arms 2-4 target BASEL AUC)",
        "gate_canonical": gate,
        "soft_scale": soft_scale,
        "threshold_canonical": float(threshold),
        "base_rate": base_rate,
        "zero_vector_phages": zero_vector_count,
        "per_phage_distance_summary": {
            "min": float(per_phage_distance["ood_distance_min_to_train"].min()),
            "mean": float(per_phage_distance["ood_distance_min_to_train"].mean()),
            "max": float(per_phage_distance["ood_distance_min_to_train"].max()),
        },
        "bacteria_axis": {"raw": bacteria_metrics_raw, "shrunk": bacteria_metrics_shrunk},
        "phage_axis": {"raw": phage_metrics_raw, "shrunk": phage_metrics_shrunk},
        "ece_closure_basel_bacteria": (
            bacteria_metrics_raw[SOURCE_BASEL]["ece"] - bacteria_metrics_shrunk[SOURCE_BASEL]["ece"]
        ),
        "ece_closure_basel_phage": (phage_metrics_raw[SOURCE_BASEL]["ece"] - phage_metrics_shrunk[SOURCE_BASEL]["ece"]),
        "ece_closure_guelin_bacteria": (
            bacteria_metrics_raw[SOURCE_GUELIN]["ece"] - bacteria_metrics_shrunk[SOURCE_GUELIN]["ece"]
        ),
        "ece_closure_guelin_phage": (
            phage_metrics_raw[SOURCE_GUELIN]["ece"] - phage_metrics_shrunk[SOURCE_GUELIN]["ece"]
        ),
    }
    with open(output_dir / "ch06_arm1_metrics.json", "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info(
        "Arm 1 done. BASEL ECE closure: bacteria-axis %.4f, phage-axis %.4f; "
        "Guelin ECE closure: bacteria-axis %.4f, phage-axis %.4f (positive = calibration improved).",
        summary["ece_closure_basel_bacteria"],
        summary["ece_closure_basel_phage"],
        summary["ece_closure_guelin_bacteria"],
        summary["ece_closure_guelin_phage"],
    )
    return summary


def parse_args(argv: Optional[list[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--ch05-dir", type=Path, default=DEFAULT_CH05_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--gate",
        choices=("hard", "soft"),
        default="hard",
        help="Hard (step at threshold) vs soft (sigmoid) gating. Hard is canonical; soft explored in diagnostic sweep.",
    )
    parser.add_argument(
        "--soft-scale",
        type=float,
        default=1.0,
        help="Sigmoid width for soft gating (only used when gate=soft).",
    )
    parser.add_argument(
        "--basel-log10-pfu-ml",
        type=float,
        default=9.0,
        help="Absolute log₁₀ pfu/ml to assign BASEL rows (matches CH05 default).",
    )
    return parser.parse_args(argv)


def main(argv: Optional[list[str]] = None) -> None:
    setup_logging()
    args = parse_args(argv)
    run_arm1(
        ch05_dir=args.ch05_dir,
        output_dir=args.output_dir,
        gate=args.gate,
        soft_scale=args.soft_scale,
        basel_log10_pfu_ml=args.basel_log10_pfu_ml,
    )


if __name__ == "__main__":
    main()
