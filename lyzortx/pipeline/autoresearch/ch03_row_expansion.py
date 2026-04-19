#!/usr/bin/env python3
"""CH03: Row-expanded training matrix loader and any_lysis rollup scaffolding.

Builds a (bacterium, phage, log_dilution, replicate, X, Y) row frame by joining
raw_interactions.csv against the ST02 pair table and ST03 split assignments. Each
raw observation (~318,816 rows across 369×96 = 35,424 pairs) carries the
pair-level metadata needed to drop into the existing SX10 pipeline once rolled
back up to pair level.

CH03 is plumbing-only: the training-time rollup collapses rows back to pair-level
`label_any_lysis` using the ST01B rule (hard_label=1 iff any row has score='1';
hard_label=0 iff no score='1' and score_0_count ≥ 5; otherwise unresolved/drop).
Rows with score='n' are treated as missing (not negative), matching SX10 label
semantics. CH04 removes the rollup and trains on the row-level binary score.

Artifacts:
  - lyzortx/generated_outputs/ch03_row_expansion/ch03_expanded_training_frame.csv
  - lyzortx/generated_outputs/ch03_row_expansion/ch03_regression_check.json
"""

from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd

from lyzortx.pipeline.autoresearch.candidate_replay import (
    DEFAULT_ST02_PAIR_TABLE_PATH,
    DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH,
)
from lyzortx.pipeline.autoresearch.sx01_eval import RAW_INTERACTIONS_PATH

LOGGER = logging.getLogger(__name__)

# ST01B label rule thresholds (reproduce lyzortx/pipeline/steel_thread_v0/steps/st01_label_policy.py).
ST01B_MIN_SCORE_0_COUNT_FOR_HARD_NEGATIVE = 5
RAW_SCORE_POSITIVE = "1"
RAW_SCORE_NEGATIVE = "0"
RAW_SCORE_UNINTERPRETABLE = "n"


def load_raw_observation_rows(raw_path: Path = RAW_INTERACTIONS_PATH) -> pd.DataFrame:
    """Load raw_interactions.csv as one row per observation.

    Columns: bacteria, bacteria_index, phage, image, replicate, plate, log_dilution,
    X, Y, score. 318,816 rows across 369 bacteria × 96 phages with log_dilution in
    {0, -1, -2, -4} (Guelin protocol: 3+2+3+1 = 9 replicates per pair).
    """
    if not raw_path.exists():
        raise FileNotFoundError(f"raw_interactions.csv not found: {raw_path}")
    df = pd.read_csv(raw_path, sep=";", dtype={"score": str})
    df["pair_id"] = df["bacteria"].astype(str) + "__" + df["phage"].astype(str)
    LOGGER.info(
        "Loaded raw observations: %d rows, %d pairs, score distribution=%s",
        len(df),
        df["pair_id"].nunique(),
        df["score"].value_counts(dropna=False).to_dict(),
    )
    return df


def load_row_expanded_frame(
    raw_path: Path = RAW_INTERACTIONS_PATH,
    st02_path: Path = DEFAULT_ST02_PAIR_TABLE_PATH,
    st03_path: Path = DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH,
) -> pd.DataFrame:
    """Merge raw observation rows with pair-level ST02/ST03 metadata.

    Every raw row gets the pair-level cv_group, split_holdout, is_hard_trainable,
    label_hard_any_lysis, training_weight_v3, and host/phage feature columns that
    downstream SX10 training reads. This is the canonical row-expanded training
    frame that CH04 will train on directly.
    """
    rows = load_raw_observation_rows(raw_path)
    st02 = pd.read_csv(st02_path, dtype=str, keep_default_na=False)
    st03 = pd.read_csv(st03_path, dtype=str, keep_default_na=False)
    pair_meta = st02.merge(
        st03[["pair_id", "split_holdout", "is_hard_trainable"]],
        on="pair_id",
        validate="one_to_one",
    )

    # Avoid duplicating bacteria/phage columns from ST02 — the raw rows already carry them.
    drop_from_pair = {"bacteria", "phage"}
    pair_columns = [c for c in pair_meta.columns if c not in drop_from_pair]
    merged = rows.merge(pair_meta[pair_columns], on="pair_id", how="left", validate="many_to_one")

    missing_pair_meta = merged[merged["cv_group"].isna() | (merged["cv_group"] == "")]
    if not missing_pair_meta.empty:
        raise ValueError(
            f"{len(missing_pair_meta)} raw rows have no ST02 pair metadata "
            f"(example pair_ids: {missing_pair_meta['pair_id'].head().tolist()})"
        )

    LOGGER.info(
        "Row-expanded frame: %d rows, %d pairs, %d columns",
        len(merged),
        merged["pair_id"].nunique(),
        len(merged.columns),
    )
    return merged


def rollup_row_scores_to_pair_any_lysis(row_frame: pd.DataFrame) -> pd.DataFrame:
    """Collapse a row-expanded frame to pair level using the ST01B any_lysis rule.

    Returns one row per pair_id with `rollup_label_hard_any_lysis` ∈ {"0", "1", ""},
    where "" marks unresolved pairs (no score='1' and fewer than 5 score='0' rows).
    Rows with score='n' are treated as missing — they count neither as lysis nor as
    clear negatives, mirroring the ST01B policy.
    """
    scores = row_frame.groupby("pair_id")["score"].agg(list)
    rollup_rows = []
    for pair_id, score_list in scores.items():
        s1 = sum(1 for v in score_list if v == RAW_SCORE_POSITIVE)
        s0 = sum(1 for v in score_list if v == RAW_SCORE_NEGATIVE)
        sn = sum(1 for v in score_list if v == RAW_SCORE_UNINTERPRETABLE)
        if s1 > 0:
            label = "1"
        elif s0 >= ST01B_MIN_SCORE_0_COUNT_FOR_HARD_NEGATIVE:
            label = "0"
        else:
            label = ""
        rollup_rows.append(
            {
                "pair_id": pair_id,
                "rollup_score_1_count": s1,
                "rollup_score_0_count": s0,
                "rollup_score_n_count": sn,
                "rollup_label_hard_any_lysis": label,
            }
        )
    return pd.DataFrame(rollup_rows)


def verify_rollup_matches_st02(
    row_frame: pd.DataFrame,
    st02_path: Path = DEFAULT_ST02_PAIR_TABLE_PATH,
) -> dict[str, int]:
    """Recompute any_lysis from row scores and confirm it matches ST02 exactly.

    This is the core CH03 regression check: the ST01B rule applied to the row-expanded
    frame must reproduce every pair's `label_hard_any_lysis` value in the canonical
    ST02 pair table. Any mismatch indicates the row expansion pipeline has drifted
    from the pair-level rollup and would contaminate CH04's per-row training.
    """
    rollup = rollup_row_scores_to_pair_any_lysis(row_frame)
    st02 = pd.read_csv(st02_path, dtype=str, keep_default_na=False)
    compared = rollup.merge(
        st02[["pair_id", "label_hard_any_lysis", "obs_score_1_count", "obs_score_0_count", "obs_score_n_count"]],
        on="pair_id",
        validate="one_to_one",
    )
    mismatches = compared[compared["rollup_label_hard_any_lysis"] != compared["label_hard_any_lysis"]]
    summary = {
        "n_pairs": len(compared),
        "n_mismatches": len(mismatches),
        "n_recomputed_positive": int((compared["rollup_label_hard_any_lysis"] == "1").sum()),
        "n_recomputed_negative": int((compared["rollup_label_hard_any_lysis"] == "0").sum()),
        "n_recomputed_unresolved": int((compared["rollup_label_hard_any_lysis"] == "").sum()),
    }
    if not mismatches.empty:
        LOGGER.error(
            "Rollup disagrees with ST02 on %d pairs; first 5: %s",
            len(mismatches),
            mismatches.head().to_dict(orient="records"),
        )
    return summary


def collapse_to_pair_level_for_training(row_frame: pd.DataFrame) -> pd.DataFrame:
    """Reduce the row-expanded frame back to the pair-level representation.

    Pair-level metadata must be constant within a pair_id. The function raises
    loudly if any pair-level column takes multiple distinct values within one
    pair, which would indicate a corrupted row-expansion join. Row-level columns
    (replicate, log_dilution, X, Y, score, image, plate) are discarded at this
    stage — CH03 keeps pair-level label semantics unchanged vs SX10. CH04
    replaces this function with a per-row trainer.
    """
    row_only = {"bacteria_index", "image", "replicate", "plate", "log_dilution", "X", "Y", "score"}
    keep = [c for c in row_frame.columns if c not in row_only]
    pair_level_columns = [c for c in keep if c != "pair_id"]
    inconsistent = row_frame.groupby("pair_id")[pair_level_columns].nunique()
    offending = inconsistent[(inconsistent > 1).any(axis=1)]
    if not offending.empty:
        first_pair = offending.index[0]
        offending_cols = offending.columns[(offending.loc[first_pair] > 1)].tolist()
        raise ValueError(f"pair-level metadata varies within pair_id {first_pair!r} for columns {offending_cols}")
    return row_frame[keep].drop_duplicates(subset=["pair_id"]).reset_index(drop=True)
