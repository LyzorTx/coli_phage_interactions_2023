#!/usr/bin/env python3
"""CH09 Arm 3 analysis: compare CH04 baseline vs neat-only-positive-filtered retrain.

Runs once the Arm 3 CH04 retrain (`ch04_eval --drop-high-titer-only-positives`)
has landed its artifacts. Loads both sets of pair predictions, computes AUC +
Brier + ECE on the held-out pairs, reports the delta, and writes a small
summary JSON.

Acceptance criterion from plan.yml Arm 3:
> If excluding high-titer-only positives lowers raw ECE by > 5 pp, the label-
> permissiveness mechanism is confirmed at training-side (not just post-hoc
> calibration).

Reads:
  - lyzortx/generated_outputs/ch04_chisel_baseline/ch04_predictions.csv
  - lyzortx/generated_outputs/ch09_calibration_layer/arm3_filtered_retrain/ch04_predictions.csv

Writes:
  - lyzortx/generated_outputs/ch09_calibration_layer/ch09_arm3_delta.json
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch09_calibration_layer import expected_calibration_error

LOGGER = logging.getLogger(__name__)

BASELINE_PREDICTIONS = Path("lyzortx/generated_outputs/ch04_chisel_baseline/ch04_predictions.csv")
ARM3_PREDICTIONS = Path("lyzortx/generated_outputs/ch09_calibration_layer/arm3_filtered_retrain/ch04_predictions.csv")
OUTPUT_PATH = Path("lyzortx/generated_outputs/ch09_calibration_layer/ch09_arm3_delta.json")

ACCEPTANCE_ECE_DELTA = 0.05  # 5 pp


def _metrics(preds: pd.DataFrame) -> dict[str, float]:
    probs = preds["predicted_probability"].to_numpy(dtype=float)
    labels = preds["label_row_binary"].to_numpy(dtype=int)
    if len(np.unique(labels)) < 2:
        return {"auc": float("nan"), "brier": float("nan"), "ece": float("nan"), "n_pairs": int(len(preds))}
    return {
        "auc": float(roc_auc_score(labels, probs)),
        "brier": float(brier_score_loss(labels, probs)),
        "ece": float(expected_calibration_error(labels, probs)),
        "n_pairs": int(len(preds)),
    }


def run_arm3_analysis(
    *,
    baseline_path: Path = BASELINE_PREDICTIONS,
    arm3_path: Path = ARM3_PREDICTIONS,
    output_path: Path = OUTPUT_PATH,
) -> dict[str, object]:
    """Diff CH04 baseline vs Arm 3 retrain metrics. Writes summary JSON; returns dict."""
    if not baseline_path.exists():
        raise FileNotFoundError(f"CH04 baseline predictions missing: {baseline_path}")
    if not arm3_path.exists():
        raise FileNotFoundError(f"CH04 Arm 3 retrain predictions missing: {arm3_path}")

    baseline = pd.read_csv(baseline_path)
    arm3 = pd.read_csv(arm3_path)
    LOGGER.info("Baseline: %d pairs. Arm 3: %d pairs.", len(baseline), len(arm3))

    # Pairs retained in the Arm 3 training set may differ — the filter drops training
    # positives but keeps all pairs in evaluation (evaluation is at max-concentration per
    # pair). So baseline and arm3 pair sets should align at the pair level — check and
    # intersect on pair_id to be safe.
    common_ids = set(baseline["pair_id"]) & set(arm3["pair_id"])
    baseline_aligned = baseline[baseline["pair_id"].isin(common_ids)].sort_values("pair_id").reset_index(drop=True)
    arm3_aligned = arm3[arm3["pair_id"].isin(common_ids)].sort_values("pair_id").reset_index(drop=True)
    LOGGER.info("Aligned on %d common pairs (baseline %d, arm3 %d).", len(common_ids), len(baseline), len(arm3))

    baseline_metrics = _metrics(baseline_aligned)
    arm3_metrics = _metrics(arm3_aligned)

    delta = {
        "auc": arm3_metrics["auc"] - baseline_metrics["auc"],
        "brier": arm3_metrics["brier"] - baseline_metrics["brier"],
        "ece": arm3_metrics["ece"] - baseline_metrics["ece"],
    }

    ece_reduction = -delta["ece"]  # positive = ECE went down = calibration improved
    label_mechanism_confirmed = ece_reduction > ACCEPTANCE_ECE_DELTA

    summary = {
        "task_id": "CH09 Arm 3",
        "baseline_predictions": str(baseline_path),
        "arm3_predictions": str(arm3_path),
        "n_common_pairs": len(common_ids),
        "baseline_metrics": baseline_metrics,
        "arm3_metrics": arm3_metrics,
        "delta": delta,
        "ece_reduction_absolute": ece_reduction,
        "acceptance_threshold_ece_reduction": ACCEPTANCE_ECE_DELTA,
        "label_permissiveness_mechanism_confirmed": bool(label_mechanism_confirmed),
        "verdict": (
            "Label-permissiveness mechanism CONFIRMED at training side "
            f"(ECE dropped {ece_reduction * 100:.1f} pp, exceeds {ACCEPTANCE_ECE_DELTA * 100:.0f} pp threshold)."
            if label_mechanism_confirmed
            else f"Label-permissiveness mechanism NOT confirmed at training side "
            f"(ECE change {ece_reduction * 100:+.1f} pp, below {ACCEPTANCE_ECE_DELTA * 100:.0f} pp threshold). "
            "Post-hoc calibration (Arm 1) remains the primary remedy."
        ),
    }

    output_path.parent.mkdir(parents=True, exist_ok=True)
    with open(output_path, "w", encoding="utf-8") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info(
        "CH09 Arm 3 verdict: baseline AUC=%.4f Brier=%.4f ECE=%.4f; arm3 AUC=%.4f Brier=%.4f ECE=%.4f; "
        "Δ AUC=%+.4f, Δ Brier=%+.4f, Δ ECE=%+.4f. Mechanism confirmed: %s",
        baseline_metrics["auc"],
        baseline_metrics["brier"],
        baseline_metrics["ece"],
        arm3_metrics["auc"],
        arm3_metrics["brier"],
        arm3_metrics["ece"],
        delta["auc"],
        delta["brier"],
        delta["ece"],
        label_mechanism_confirmed,
    )
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-predictions", type=Path, default=BASELINE_PREDICTIONS)
    parser.add_argument("--arm3-predictions", type=Path, default=ARM3_PREDICTIONS)
    parser.add_argument("--output-path", type=Path, default=OUTPUT_PATH)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    run_arm3_analysis(
        baseline_path=args.baseline_predictions,
        arm3_path=args.arm3_predictions,
        output_path=args.output_path,
    )


if __name__ == "__main__":
    main()
