"""Tests for cv_group-based fold assignment in sx01_eval."""

from __future__ import annotations

import pandas as pd
import pytest

from lyzortx.pipeline.autoresearch.sx01_eval import (
    FOLD_SALT,
    N_FOLDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)


def test_bacteria_sharing_cv_group_land_in_same_fold() -> None:
    mapping = {
        "bac_A": "group_1",
        "bac_B": "group_1",
        "bac_C": "group_2",
        "bac_D": "group_2",
        "bac_E": "group_3",
    }
    folds = assign_bacteria_folds(mapping)
    assert folds["bac_A"] == folds["bac_B"], "bacteria in group_1 must share a fold"
    assert folds["bac_C"] == folds["bac_D"], "bacteria in group_2 must share a fold"
    assert all(0 <= f < N_FOLDS for f in folds.values())


def test_assignment_is_deterministic() -> None:
    mapping = {"bac_A": "g1", "bac_B": "g2", "bac_C": "g3"}
    first = assign_bacteria_folds(mapping)
    second = assign_bacteria_folds(mapping)
    assert first == second


def test_different_salt_changes_assignment() -> None:
    mapping = {"bac_A": "g1", "bac_B": "g2", "bac_C": "g3", "bac_D": "g4"}
    default = assign_bacteria_folds(mapping, salt=FOLD_SALT)
    alt = assign_bacteria_folds(mapping, salt="different_salt")
    assert default != alt


def test_empty_mapping_raises() -> None:
    with pytest.raises(ValueError, match="empty"):
        assign_bacteria_folds({})


def test_missing_cv_group_raises() -> None:
    with pytest.raises(ValueError, match="no cv_group"):
        assign_bacteria_folds({"bac_A": "g1", "bac_B": ""})


def test_bacteria_to_cv_group_map_builds_from_frame() -> None:
    frame = pd.DataFrame(
        {
            "bacteria": ["bac_A", "bac_A", "bac_B", "bac_C"],
            "cv_group": ["g1", "g1", "g1", "g2"],
            "phage": ["p1", "p2", "p1", "p1"],
        }
    )
    mapping = bacteria_to_cv_group_map(frame)
    assert mapping == {"bac_A": "g1", "bac_B": "g1", "bac_C": "g2"}


def test_bacteria_to_cv_group_map_rejects_conflicts() -> None:
    frame = pd.DataFrame(
        {
            "bacteria": ["bac_A", "bac_A"],
            "cv_group": ["g1", "g2"],
        }
    )
    with pytest.raises(ValueError, match="multiple cv_group"):
        bacteria_to_cv_group_map(frame)


def test_bacteria_to_cv_group_map_missing_column() -> None:
    frame = pd.DataFrame({"bacteria": ["bac_A"]})
    with pytest.raises(KeyError, match="cv_group"):
        bacteria_to_cv_group_map(frame)
