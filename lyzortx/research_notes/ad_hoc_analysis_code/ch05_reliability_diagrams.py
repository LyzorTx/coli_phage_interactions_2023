#!/usr/bin/env python3
"""CH05 post-hoc: per-decile reliability tables for Guelin + BASEL on both axes.

Reads the existing CH05 bacteria-axis and phage-axis predictions, bins predicted
probability into 10 deciles, and reports observed lysis rate vs mean predicted
probability per decile per (source × axis). This quantifies the calibration
divergence flagged by the second reviewer: BASEL mid-P reliability was 30-55 pp off
vs Guelin, and the Brier gap survives Straboviridae exclusion.

No model rerun — reads predictions CSVs. Outputs:
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_reliability_tables.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import DEFAULT_OUTPUT_DIR, SOURCE_BASEL, SOURCE_GUELIN

LOGGER = logging.getLogger(__name__)

N_BINS = 10
OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_reliability_tables.csv"


def _reliability_rows(axis: str, source: str, df: pd.DataFrame) -> list[dict[str, object]]:
    preds = df["predicted_probability"].astype(float).to_numpy()
    labels = df["label_row_binary"].astype(int).to_numpy()
    # Fixed decile edges on the [0, 1] predicted-probability range so Guelin and BASEL
    # bins are comparable (quantile bins would recenter by source and obscure the gap).
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
                "axis": axis,
                "source": source,
                "bin_lo": round(float(edges[b]), 3),
                "bin_hi": round(float(edges[b + 1]), 3),
                "n_pairs": n,
                "mean_predicted": round(mean_pred, 4),
                "observed_rate": round(obs_rate, 4),
                "reliability_gap": round(obs_rate - mean_pred, 4),
            }
        )
    return rows


def main() -> None:
    setup_logging()
    rows: list[dict[str, object]] = []
    for axis in ("bacteria_axis", "phage_axis"):
        preds = pd.read_csv(Path(DEFAULT_OUTPUT_DIR) / f"ch05_{axis}_predictions.csv")
        for source in (SOURCE_GUELIN, SOURCE_BASEL):
            subset = preds[preds["source"] == source]
            if subset.empty:
                continue
            rows.extend(_reliability_rows(axis, source, subset))

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    LOGGER.info("Wrote %s (%d rows)", OUTPUT_CSV, len(df))

    print()
    print("=== CH05 reliability by decile (|gap| = |observed − predicted|) ===")
    for axis in ("bacteria_axis", "phage_axis"):
        print(f"\n--- {axis} ---")
        subset = df[df["axis"] == axis]
        print(subset.to_string(index=False))

    print("\n--- Cross-source reliability delta per bin (BASEL − Guelin) ---")
    for axis in ("bacteria_axis", "phage_axis"):
        gue = df[(df["axis"] == axis) & (df["source"] == SOURCE_GUELIN)].set_index(["bin_lo", "bin_hi"])
        bas = df[(df["axis"] == axis) & (df["source"] == SOURCE_BASEL)].set_index(["bin_lo", "bin_hi"])
        joined = gue.join(bas, lsuffix="_g", rsuffix="_b", how="inner")
        if joined.empty:
            continue
        joined["delta_observed"] = joined["observed_rate_b"] - joined["observed_rate_g"]
        joined["delta_gap"] = joined["reliability_gap_b"] - joined["reliability_gap_g"]
        print(f"\n{axis}:")
        print(
            joined[
                [
                    "n_pairs_g",
                    "n_pairs_b",
                    "observed_rate_g",
                    "observed_rate_b",
                    "delta_observed",
                    "reliability_gap_g",
                    "reliability_gap_b",
                    "delta_gap",
                ]
            ]
            .round(4)
            .to_string()
        )


if __name__ == "__main__":
    main()
