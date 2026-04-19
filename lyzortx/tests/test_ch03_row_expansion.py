"""Tests for CH03 row-expansion and any_lysis rollup scaffolding."""

from __future__ import annotations

import pandas as pd
import pytest

from lyzortx.pipeline.autoresearch.ch03_row_expansion import (
    ST01B_MIN_SCORE_0_COUNT_FOR_HARD_NEGATIVE,
    collapse_to_pair_level_for_training,
    rollup_row_scores_to_pair_any_lysis,
)


def _row(pair_id: str, score: str, replicate: int = 1, log_dilution: int = 0) -> dict:
    bacteria, phage = pair_id.split("__")
    return {
        "pair_id": pair_id,
        "bacteria": bacteria,
        "phage": phage,
        "replicate": replicate,
        "log_dilution": log_dilution,
        "score": score,
    }


def test_rollup_positive_any_lysis_when_any_score_one() -> None:
    rows = pd.DataFrame(
        [
            _row("bac_A__phage_X", "0"),
            _row("bac_A__phage_X", "0"),
            _row("bac_A__phage_X", "1"),
            _row("bac_A__phage_X", "0"),
        ]
    )
    out = rollup_row_scores_to_pair_any_lysis(rows)
    row = out.iloc[0]
    assert row["rollup_label_hard_any_lysis"] == "1"
    assert row["rollup_score_1_count"] == 1
    assert row["rollup_score_0_count"] == 3


def test_rollup_hard_negative_requires_min_score_zero_count() -> None:
    # Need ≥5 score='0' rows to call a hard 0.
    n = ST01B_MIN_SCORE_0_COUNT_FOR_HARD_NEGATIVE
    rows_hard = pd.DataFrame([_row("bac_A__phage_X", "0", replicate=i) for i in range(n)])
    assert rollup_row_scores_to_pair_any_lysis(rows_hard).iloc[0]["rollup_label_hard_any_lysis"] == "0"

    rows_soft = pd.DataFrame([_row("bac_A__phage_X", "0", replicate=i) for i in range(n - 1)])
    assert rollup_row_scores_to_pair_any_lysis(rows_soft).iloc[0]["rollup_label_hard_any_lysis"] == ""


def test_rollup_treats_n_as_missing_not_negative() -> None:
    # 4 score='n' rows — should NOT become hard 0 because 'n' is missing, not negative.
    rows = pd.DataFrame([_row("bac_A__phage_X", "n", replicate=i) for i in range(9)])
    out = rollup_row_scores_to_pair_any_lysis(rows)
    row = out.iloc[0]
    assert row["rollup_label_hard_any_lysis"] == ""
    assert row["rollup_score_n_count"] == 9
    assert row["rollup_score_0_count"] == 0


def test_rollup_score_one_overrides_insufficient_zero_count() -> None:
    # Single score='1' beats any number of 'n's and 0s.
    rows = pd.DataFrame(
        [
            _row("bac_A__phage_X", "n"),
            _row("bac_A__phage_X", "n"),
            _row("bac_A__phage_X", "1"),
        ]
    )
    assert rollup_row_scores_to_pair_any_lysis(rows).iloc[0]["rollup_label_hard_any_lysis"] == "1"


def test_collapse_deduplicates_rows_to_pair_level() -> None:
    rows = pd.DataFrame(
        [
            {
                "pair_id": "bac_A__phage_X",
                "bacteria": "bac_A",
                "phage": "phage_X",
                "cv_group": "42",
                "label_hard_any_lysis": "1",
                "training_weight_v3": "1.0",
                "replicate": 1,
                "log_dilution": 0,
                "score": "0",
            },
            {
                "pair_id": "bac_A__phage_X",
                "bacteria": "bac_A",
                "phage": "phage_X",
                "cv_group": "42",
                "label_hard_any_lysis": "1",
                "training_weight_v3": "1.0",
                "replicate": 2,
                "log_dilution": -1,
                "score": "1",
            },
            {
                "pair_id": "bac_B__phage_X",
                "bacteria": "bac_B",
                "phage": "phage_X",
                "cv_group": "42",
                "label_hard_any_lysis": "0",
                "training_weight_v3": "1.0",
                "replicate": 1,
                "log_dilution": 0,
                "score": "0",
            },
        ]
    )
    collapsed = collapse_to_pair_level_for_training(rows)
    assert len(collapsed) == 2
    assert set(collapsed["pair_id"]) == {"bac_A__phage_X", "bac_B__phage_X"}
    # Row-level columns must be stripped from the collapsed frame.
    assert {"replicate", "log_dilution", "score"}.isdisjoint(collapsed.columns)
    # Pair-level metadata survives.
    assert {"pair_id", "bacteria", "phage", "cv_group", "label_hard_any_lysis"}.issubset(collapsed.columns)


def test_collapse_rejects_pair_with_inconsistent_metadata() -> None:
    rows = pd.DataFrame(
        [
            {
                "pair_id": "bac_A__phage_X",
                "bacteria": "bac_A",
                "phage": "phage_X",
                "cv_group": "42",
                "replicate": 1,
                "score": "0",
            },
            {
                "pair_id": "bac_A__phage_X",
                "bacteria": "bac_A",
                "phage": "phage_X",
                "cv_group": "43",
                "replicate": 2,
                "score": "1",
            },
        ]
    )
    with pytest.raises(ValueError, match="pair-level metadata varies"):
        collapse_to_pair_level_for_training(rows)
