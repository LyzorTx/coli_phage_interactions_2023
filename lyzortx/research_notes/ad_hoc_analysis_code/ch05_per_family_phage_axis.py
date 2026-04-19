#!/usr/bin/env python3
"""CH05 post-hoc: per-family phage-axis AUC + within-family vs family-novel contrast.

Motivated by self-review observation that CH05's phage-axis StratifiedKFold by ICTV
family holds out ~1 phage per fold for families with ≥10 members, keeping siblings in
training — i.e. the aggregate AUC is biased toward *within-family* generalization and is
an optimistic proxy for the "completely new family" deployment case.

Rare families (<10 phages) are collapsed into an "other" bucket at fold-assignment time,
so held-out phages from those families have NO same-family siblings in training — that's
the family-novel regime. Computing per-family AUC lets us quantify the contrast.

Outputs: prints a table to stdout and writes
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_per_family_breakdown.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import DEFAULT_OUTPUT_DIR
from lyzortx.pipeline.autoresearch.sx01_eval import N_FOLDS
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

MIN_PAIRS_FOR_REPORTING = 10
PHAGE_AXIS_PREDICTIONS = Path(DEFAULT_OUTPUT_DIR) / "ch05_phage_axis_predictions.csv"
OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_per_family_breakdown.csv"


def main() -> None:
    setup_logging()
    preds = pd.read_csv(PHAGE_AXIS_PREDICTIONS)
    family_map = load_unified_phage_family_map()
    preds["family"] = preds["phage"].map(family_map).fillna("UNKNOWN")

    # Global panel counts determine whether a family is within-family (≥N_FOLDS) or
    # family-novel (collapsed to "other" bucket in assign_phage_folds).
    from collections import Counter

    global_counts = Counter(family_map.values())
    rows = []
    for family in sorted(preds["family"].unique()):
        subset = preds[preds["family"] == family]
        n_pairs = len(subset)
        n_phages = int(subset["phage"].nunique())
        panel_n = int(global_counts.get(family, 0))
        if n_pairs < MIN_PAIRS_FOR_REPORTING or subset["label_row_binary"].nunique() < 2:
            continue
        regime = "within-family" if panel_n >= N_FOLDS else "family-novel (rare, collapsed to 'other')"
        auc = float(roc_auc_score(subset["label_row_binary"], subset["predicted_probability"]))
        brier = float(brier_score_loss(subset["label_row_binary"], subset["predicted_probability"]))
        rows.append(
            {
                "family": family,
                "generalization_regime": regime,
                "panel_phage_count": panel_n,
                "n_phages_scored": n_phages,
                "n_pairs": n_pairs,
                "pos_rate": round(float(subset["label_row_binary"].mean()), 4),
                "auc": round(auc, 4),
                "brier": round(brier, 4),
            }
        )
    df = pd.DataFrame(rows).sort_values(["generalization_regime", "panel_phage_count"], ascending=[True, False])
    df.to_csv(OUTPUT_CSV, index=False)
    LOGGER.info("Wrote %s (%d rows)", OUTPUT_CSV, len(df))

    within_mask = df["generalization_regime"] == "within-family"
    print()
    print("=== Per-family phage-axis AUC/Brier ===")
    print(df.to_string(index=False))
    print()
    if within_mask.any() and (~within_mask).any():
        within_auc_weighted = (df.loc[within_mask, "auc"] * df.loc[within_mask, "n_pairs"]).sum() / df.loc[
            within_mask, "n_pairs"
        ].sum()
        novel_auc_weighted = (df.loc[~within_mask, "auc"] * df.loc[~within_mask, "n_pairs"]).sum() / df.loc[
            ~within_mask, "n_pairs"
        ].sum()
        print(
            f"Pair-weighted AUC — within-family regime: {within_auc_weighted:.4f}; "
            f"family-novel regime: {novel_auc_weighted:.4f}; "
            f"gap: {within_auc_weighted - novel_auc_weighted:+.4f}"
        )


if __name__ == "__main__":
    main()
