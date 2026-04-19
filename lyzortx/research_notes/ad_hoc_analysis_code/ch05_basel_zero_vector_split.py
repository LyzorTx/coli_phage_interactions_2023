#!/usr/bin/env python3
"""CH05 post-hoc: split BASEL phages by zero-vector vs non-zero on phage_projection.

Tests whether the BASEL mid-P over-prediction localizes to the 13 BASEL phages with
all-zero phage_projection vectors (no TL17 BLAST hits) or spreads across non-zero
BASEL phages too. If localizes → "fill in 13 zero-vector phages" is the fix; if
spreads → TL17 reference bank itself is panel-biased.

Also checks whether the 13 zero-vector BASEL phages overlap with the UNKNOWN
(no ICTV family) BASEL phages — i.e., is "taxonomically novel" the same set as
"no TL17 hit"?

Reads predictions CSVs + feature cache. Outputs:
  lyzortx/generated_outputs/ch05_unified_kfold/ch05_basel_zero_vector_split.csv
"""

from __future__ import annotations

import logging
from pathlib import Path

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

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
from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

LOGGER = logging.getLogger(__name__)

OUTPUT_CSV = Path(DEFAULT_OUTPUT_DIR) / "ch05_basel_zero_vector_split.csv"


def _summary(label: str, df: pd.DataFrame) -> dict[str, object]:
    if df.empty:
        return {"label": label, "n_pairs": 0, "n_phages": 0}
    labels = df["label_row_binary"].astype(int).to_numpy()
    preds = df["predicted_probability"].astype(float).to_numpy()
    pos_mask = labels == 1
    neg_mask = labels == 0
    auc = float(roc_auc_score(labels, preds)) if len(set(labels.tolist())) > 1 else float("nan")
    brier = float(brier_score_loss(labels, preds))
    return {
        "label": label,
        "n_phages": int(df["phage"].nunique()),
        "n_pairs": len(df),
        "pos_rate": round(float(pos_mask.mean()), 4),
        "mean_pred_overall": round(float(preds.mean()), 4),
        "mean_pred_pos": round(float(preds[pos_mask].mean()), 4) if pos_mask.any() else float("nan"),
        "mean_pred_neg": round(float(preds[neg_mask].mean()), 4) if neg_mask.any() else float("nan"),
        "auc": round(auc, 4),
        "brier": round(brier, 4),
    }


def main() -> None:
    setup_logging()
    unified = load_unified_row_frame()
    phage_source = unified.groupby("phage")["source"].first().to_dict()

    candidate_module = load_module_from_path("ch05_candidate", Path(DEFAULT_CANDIDATE_DIR) / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=Path(DEFAULT_CACHE_DIR), include_host_defense=True)
    patch_context_with_extended_slots(context)

    feature_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=["phage_projection"], entity_key="phage"
    ).set_index("phage")
    numeric = feature_table.select_dtypes(include=[np.number])
    row_abs_sum = numeric.abs().sum(axis=1)
    zero_vector_phages = set(row_abs_sum[row_abs_sum < 1e-12].index)

    basel_phages = {p for p, src in phage_source.items() if src == SOURCE_BASEL}
    guelin_phages = {p for p, src in phage_source.items() if src == SOURCE_GUELIN}
    basel_zero = basel_phages & zero_vector_phages
    basel_nonzero = basel_phages - zero_vector_phages
    guelin_zero = guelin_phages & zero_vector_phages
    guelin_nonzero = guelin_phages - zero_vector_phages

    family_map = load_unified_phage_family_map()
    basel_unknown = {p for p in basel_phages if p not in family_map or not family_map[p]}
    overlap = basel_zero & basel_unknown

    LOGGER.info(
        "BASEL phage split: %d zero-vector, %d non-zero; Guelin: %d zero-vector, %d non-zero",
        len(basel_zero),
        len(basel_nonzero),
        len(guelin_zero),
        len(guelin_nonzero),
    )
    LOGGER.info(
        "BASEL zero-vector ∩ BASEL UNKNOWN-family = %d of %d zero-vector / %d UNKNOWN",
        len(overlap),
        len(basel_zero),
        len(basel_unknown),
    )
    if basel_zero:
        LOGGER.info("BASEL zero-vector phages: %s", sorted(basel_zero))

    rows: list[dict[str, object]] = []
    for axis in ("bacteria_axis", "phage_axis"):
        preds_path = Path(DEFAULT_OUTPUT_DIR) / f"ch05_{axis}_predictions.csv"
        preds = pd.read_csv(preds_path)
        for label, subset_phages in [
            (f"{axis} | Guelin", guelin_phages),
            (f"{axis} | BASEL", basel_phages),
            (f"{axis} | BASEL zero-vector ({len(basel_zero)} phages)", basel_zero),
            (f"{axis} | BASEL non-zero ({len(basel_nonzero)} phages)", basel_nonzero),
        ]:
            df = preds[preds["phage"].isin(subset_phages)]
            rows.append(_summary(label, df))

    df = pd.DataFrame(rows)
    df.to_csv(OUTPUT_CSV, index=False)
    print()
    print("=== BASEL zero-vector vs non-zero reliability split ===")
    print(df.to_string(index=False))
    print()
    print(
        f"BASEL zero-vector ∩ UNKNOWN-family overlap: {len(overlap)}/{len(basel_zero)} zero-vector, "
        f"{len(overlap)}/{len(basel_unknown)} UNKNOWN-family"
    )
    print(f"Zero-vector BASEL phages: {sorted(basel_zero)}")
    print(f"UNKNOWN-family BASEL phages (no ICTV): {sorted(basel_unknown)}")


if __name__ == "__main__":
    main()
