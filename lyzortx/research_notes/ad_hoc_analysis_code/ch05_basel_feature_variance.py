#!/usr/bin/env python3
"""CH05 post-hoc: Guelin vs BASEL per-phage feature-slot variance comparison.

Motivated by reviewer hypothesis: BASEL's mid-P over-prediction (reliability gap −21
to −27 pp worse than Guelin in 0.5-0.9 bins) may be driven by zero-fill-heavy BASEL
phage features. If BASEL phages carry systematically lower variance across
`phage_projection`, `phage_stats`, and any other phage-side slots, the all-pairs
model would fall back to a "default phage prior" — pushing BASEL predictions toward
mid-range values regardless of actual host compatibility.

For each phage slot, compute:
  - n_phages per source
  - mean per-feature coefficient of variation (CV) per source
  - fraction of features with variance < 1e-12 (i.e. effectively constant) per source
  - fraction of rows that are exact zero vectors per source

Reads the feature cache — no model rerun. Outputs:
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_basel_feature_variance.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import load_module_from_path
from lyzortx.pipeline.autoresearch.ch05_eval import (
    DEFAULT_CACHE_DIR,
    DEFAULT_CANDIDATE_DIR,
    DEFAULT_OUTPUT_DIR,
    SOURCE_BASEL,
    SOURCE_GUELIN,
    load_unified_row_frame,
)
from lyzortx.pipeline.autoresearch.sx03_eval import patch_context_with_extended_slots

LOGGER = logging.getLogger(__name__)

PHAGE_SLOTS = ["phage_projection", "phage_stats", "phage_rbp_struct"]
OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_basel_feature_variance.csv"
VAR_FLOOR = 1e-12


def _summarise_slot(slot_name: str, source_label: str, features: pd.DataFrame) -> dict[str, object]:
    n_phages, n_features = features.shape
    arr = features.to_numpy(dtype=float)
    per_feat_var = arr.var(axis=0)
    per_feat_mean = arr.mean(axis=0)
    safe_mean = np.where(np.abs(per_feat_mean) < 1e-12, np.nan, per_feat_mean)
    per_feat_cv = np.sqrt(per_feat_var) / np.abs(safe_mean)
    n_constant = int((per_feat_var < VAR_FLOOR).sum())
    n_zero_rows = int((np.abs(arr).sum(axis=1) < VAR_FLOOR).sum())
    return {
        "slot": slot_name,
        "source": source_label,
        "n_phages": n_phages,
        "n_features": n_features,
        "frac_constant_features": round(n_constant / max(n_features, 1), 4),
        "frac_zero_rows": round(n_zero_rows / max(n_phages, 1), 4),
        "mean_per_feature_cv": round(float(np.nanmean(per_feat_cv)), 4),
        "median_per_feature_cv": round(float(np.nanmedian(per_feat_cv)), 4),
        "mean_per_feature_variance": round(float(per_feat_var.mean()), 6),
    }


def main() -> None:
    setup_logging()
    unified = load_unified_row_frame()
    phage_source = unified.groupby("phage")["source"].first().to_dict()

    candidate_module = load_module_from_path("ch05_candidate", Path(DEFAULT_CANDIDATE_DIR) / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=Path(DEFAULT_CACHE_DIR), include_host_defense=True)
    patch_context_with_extended_slots(context)

    rows: list[dict[str, object]] = []
    for slot_name in PHAGE_SLOTS:
        try:
            feature_table = candidate_module.build_entity_feature_table(
                context.slot_artifacts, slot_names=[slot_name], entity_key="phage"
            )
        except Exception as exc:  # slot may not be present in cache
            LOGGER.warning("Skipping slot %s: %s", slot_name, exc)
            continue
        feature_table = feature_table.set_index("phage")
        feature_table = feature_table.select_dtypes(include=[np.number])
        guelin_phages = [p for p in feature_table.index if phage_source.get(p) == SOURCE_GUELIN]
        basel_phages = [p for p in feature_table.index if phage_source.get(p) == SOURCE_BASEL]
        if guelin_phages:
            rows.append(_summarise_slot(slot_name, SOURCE_GUELIN, feature_table.loc[guelin_phages]))
        if basel_phages:
            rows.append(_summarise_slot(slot_name, SOURCE_BASEL, feature_table.loc[basel_phages]))

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    LOGGER.info("Wrote %s (%d rows)", OUTPUT_CSV, len(df))

    print()
    print("=== CH05 Guelin vs BASEL per-phage feature-slot variance ===")
    print(df.to_string(index=False))
    print()

    # Per-slot delta: higher frac_constant / lower CV on BASEL → zero-fill hypothesis supported.
    print("=== Cross-source deltas (BASEL − Guelin) ===")
    pivot = df.pivot_table(
        index="slot",
        columns="source",
        values=["frac_constant_features", "frac_zero_rows", "mean_per_feature_cv"],
    )
    print(pivot.round(4).to_string())


if __name__ == "__main__":
    main()
