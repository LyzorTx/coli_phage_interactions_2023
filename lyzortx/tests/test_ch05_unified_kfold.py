"""Tests for CH05 unified Guelin+BASEL frame and phage-axis fold assignment."""

from __future__ import annotations

from unittest.mock import patch

import pandas as pd

from lyzortx.pipeline.autoresearch.ch05_eval import (
    BASEL_LOG_DILUTION,
    BASEL_REPLICATE_ID,
    SOURCE_BASEL,
    SOURCE_GUELIN,
    assign_phage_folds,
    load_basel_as_row_frame,
)


def test_assign_phage_folds_is_deterministic_and_stratifies() -> None:
    phages = [f"p{i:02d}" for i in range(20)]
    family = {p: ("famA" if int(p[1:]) < 10 else "famB") for p in phages}
    folds_a = assign_phage_folds(phages, family, n_splits=5)
    folds_b = assign_phage_folds(phages, family, n_splits=5)
    assert folds_a == folds_b
    # Each fold should carry both families roughly evenly.
    per_fold_family: dict[int, dict[str, int]] = {}
    for phage, fold in folds_a.items():
        per_fold_family.setdefault(fold, {"famA": 0, "famB": 0})[family[phage]] += 1
    for fold, counts in per_fold_family.items():
        assert counts["famA"] == 2, f"Fold {fold} famA count {counts['famA']}"
        assert counts["famB"] == 2, f"Fold {fold} famB count {counts['famB']}"


def test_assign_phage_folds_collapses_rare_families() -> None:
    # 10 phages in famA, 1 in famB, 1 in famC (both rare) — rares collapse to "other".
    phages = [f"p{i:02d}" for i in range(12)]
    family = {phages[i]: "famA" for i in range(10)}
    family[phages[10]] = "famB"
    family[phages[11]] = "famC"
    folds = assign_phage_folds(phages, family, n_splits=5)
    # Every phage receives a fold, assignment covers all 12.
    assert set(folds.keys()) == set(phages)
    assert set(folds.values()) == set(range(5))


def test_assign_phage_folds_handles_missing_family() -> None:
    phages = [f"p{i:02d}" for i in range(10)]
    # Only half the phages have a family entry.
    family = {phages[i]: "famA" for i in range(5)}
    folds = assign_phage_folds(phages, family, n_splits=5)
    assert set(folds.keys()) == set(phages)
    # All folds populated.
    assert set(folds.values()) == set(range(5))


def test_load_basel_as_row_frame_inherits_cv_group_from_guelin() -> None:
    basel_fake = pd.DataFrame(
        [
            {
                "pair_id": "ECOR-12__phage_X",
                "bacteria": "ECOR-12",
                "phage": "phage_X",
                "interaction": 1,
                "source": "basel",
            },
            {
                "pair_id": "ECOR-13__phage_Y",
                "bacteria": "ECOR-13",
                "phage": "phage_Y",
                "interaction": 0,
                "source": "basel",
            },
            {
                "pair_id": "UNKNOWN__phage_Z",
                "bacteria": "UNKNOWN",
                "phage": "phage_Z",
                "interaction": 1,
                "source": "basel",
            },
        ]
    )
    guelin_map = {"ECOR-12": "cv42", "ECOR-13": "cv42"}  # UNKNOWN absent → row dropped
    with patch("lyzortx.pipeline.autoresearch.ch05_eval.load_basel_interactions", return_value=basel_fake):
        result = load_basel_as_row_frame(guelin_map)
    assert len(result) == 2  # UNKNOWN dropped
    assert set(result["cv_group"]) == {"cv42"}
    assert set(result["log_dilution"]) == {BASEL_LOG_DILUTION}
    assert set(result["replicate"]) == {BASEL_REPLICATE_ID}
    assert set(result["source"]) == {SOURCE_BASEL}
    assert set(result["split_holdout"]) == {"train_non_holdout"}
    assert set(result["is_hard_trainable"]) == {"1"}
    positive = result[result["pair_id"] == "ECOR-12__phage_X"].iloc[0]
    assert positive["score"] == "1"
    negative = result[result["pair_id"] == "ECOR-13__phage_Y"].iloc[0]
    assert negative["score"] == "0"


def test_source_constants_are_distinct() -> None:
    assert SOURCE_GUELIN != SOURCE_BASEL
    assert SOURCE_GUELIN.lower() == "guelin"
    assert SOURCE_BASEL.lower() == "basel"
