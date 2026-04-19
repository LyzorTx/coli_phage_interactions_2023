#!/usr/bin/env python3
"""CH05 post-hoc: Straboviridae-exclusion diagnostic on existing predictions.

Motivated by reviewer observation that the CH05 cross-source BASEL AUC gap + BASEL
calibration divergence may trace to Straboviridae prior collapse (family-bias-straboviridae):
broad-host Straboviridae (62-71% lysis rate) dominate model rankings and may be over-predicted
on BASEL panels that happen to carry different Straboviridae proportions than Guelin.

Reads the existing CH05 phage-axis and bacteria-axis predictions + the unified ICTV family
map. Recomputes aggregate AUC/Brier with and without Straboviridae for each (source × axis)
combination. A material shrink in the BASEL-minus-Guelin delta after Straboviridae exclusion
would flag the family as the operative confounder; a stable delta rules it out.

No model rerun — reads predictions CSVs + family map. Outputs:
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_straboviridae_exclusion.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.ch05_eval import DEFAULT_OUTPUT_DIR, SOURCE_BASEL, SOURCE_GUELIN
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

STRABOVIRIDAE_FAMILY = "Straboviridae"
OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_straboviridae_exclusion.csv"


def _auc_brier(df: pd.DataFrame) -> tuple[float, float, int, int, int]:
    labels = df["label_row_binary"].astype(int).to_numpy()
    preds = df["predicted_probability"].astype(float).to_numpy()
    n_pairs = len(df)
    n_pos = int(labels.sum())
    n_phages = int(df["phage"].nunique())
    auc = float(roc_auc_score(labels, preds)) if len(set(labels.tolist())) > 1 else float("nan")
    brier = float(brier_score_loss(labels, preds))
    return auc, brier, n_pairs, n_pos, n_phages


def _row(axis: str, source: str | None, include_strabo: bool, df: pd.DataFrame) -> dict[str, object]:
    auc, brier, n_pairs, n_pos, n_phages = _auc_brier(df)
    return {
        "axis": axis,
        "source": source or "combined",
        "straboviridae": "included" if include_strabo else "excluded",
        "n_phages": n_phages,
        "n_pairs": n_pairs,
        "pos_rate": round(n_pos / n_pairs, 4) if n_pairs else float("nan"),
        "auc": round(auc, 4),
        "brier": round(brier, 4),
    }


def main() -> None:
    setup_logging()
    family_map = load_unified_phage_family_map()

    rows: list[dict[str, object]] = []
    for axis in ("bacteria_axis", "phage_axis"):
        preds_path = Path(DEFAULT_OUTPUT_DIR) / f"ch05_{axis}_predictions.csv"
        preds = pd.read_csv(preds_path)
        preds["family"] = preds["phage"].map(family_map).fillna("UNKNOWN")

        for source in (None, SOURCE_GUELIN, SOURCE_BASEL):
            subset_all = preds if source is None else preds[preds["source"] == source]
            if subset_all.empty:
                continue
            rows.append(_row(axis, source, True, subset_all))
            subset_no_strabo = subset_all[subset_all["family"] != STRABOVIRIDAE_FAMILY]
            if not subset_no_strabo.empty and subset_no_strabo["label_row_binary"].nunique() > 1:
                rows.append(_row(axis, source, False, subset_no_strabo))

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    LOGGER.info("Wrote %s (%d rows)", OUTPUT_CSV, len(df))

    print()
    print("=== CH05 Straboviridae-exclusion diagnostic (existing predictions, no rerun) ===")
    print(df.to_string(index=False))
    print()

    # Cross-source delta shift: does excluding Straboviridae close the BASEL-minus-Guelin AUC gap?
    for axis in ("bacteria_axis", "phage_axis"):
        pivot = df[df["axis"] == axis].pivot_table(
            index="straboviridae", columns="source", values="auc", aggfunc="first"
        )
        if {"guelin", "basel"}.issubset(pivot.columns):
            delta_incl = float(pivot.loc["included", "basel"] - pivot.loc["included", "guelin"])
            delta_excl = float(pivot.loc["excluded", "basel"] - pivot.loc["excluded", "guelin"])
            print(
                f"{axis}: BASEL − Guelin AUC | Straboviridae included: {delta_incl:+.4f}, "
                f"excluded: {delta_excl:+.4f}, shift: {delta_excl - delta_incl:+.4f}"
            )


if __name__ == "__main__":
    main()
