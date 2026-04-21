"""Tests for CH05 unified Guelin+BASEL frame and phage-axis fold assignment."""

from __future__ import annotations

from pathlib import Path
from unittest.mock import patch

import pandas as pd
import pytest

from lyzortx.pipeline.autoresearch.candidate_replay import BootstrapMetricCI
from lyzortx.pipeline.autoresearch.ch05_eval import (
    BASEL_LOG10_PFU_ML,
    BASEL_LOG_DILUTION_SENTINEL,
    BASEL_REPLICATE_ID,
    SOURCE_BASEL,
    _bootstrap_by_unit,
    assign_phage_folds,
    load_basel_as_row_frame,
    run_ch05_eval,
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
        ]
    )
    guelin_map = {"ECOR-12": "cv42", "ECOR-13": "cv43"}
    with patch("lyzortx.pipeline.autoresearch.ch05_eval.load_basel_interactions", return_value=basel_fake):
        result = load_basel_as_row_frame(guelin_map)
    assert len(result) == 2
    assert set(result["cv_group"]) == {"cv42", "cv43"}
    assert set(result["log_dilution"]) == {BASEL_LOG_DILUTION_SENTINEL}
    assert set(result["log10_pfu_ml"]) == {BASEL_LOG10_PFU_ML}
    assert set(result["replicate"]) == {BASEL_REPLICATE_ID}
    assert set(result["source"]) == {SOURCE_BASEL}
    assert set(result["split_holdout"]) == {"train_non_holdout"}
    assert set(result["is_hard_trainable"]) == {"1"}
    positive = result[result["pair_id"] == "ECOR-12__phage_X"].iloc[0]
    assert positive["score"] == "1"
    negative = result[result["pair_id"] == "ECOR-13__phage_Y"].iloc[0]
    assert negative["score"] == "0"


def test_load_basel_as_row_frame_fails_fast_on_missing_cv_group() -> None:
    """AGENTS.md fail-fast: unexpected empty join must raise, not warn+drop."""
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
                "pair_id": "UNKNOWN_BACTERIUM__phage_Y",
                "bacteria": "UNKNOWN_BACTERIUM",
                "phage": "phage_Y",
                "interaction": 0,
                "source": "basel",
            },
        ]
    )
    guelin_map = {"ECOR-12": "cv42"}
    with patch("lyzortx.pipeline.autoresearch.ch05_eval.load_basel_interactions", return_value=basel_fake):
        with pytest.raises(ValueError, match="no Guelin cv_group mapping"):
            load_basel_as_row_frame(guelin_map)


def test_bootstrap_by_unit_resamples_at_specified_unit_level() -> None:
    """Cluster-bootstrap: resample bacteria (not pair rows) → CI width reflects
    the independent sample size at bacterium level, not pair level."""
    # Two bacteria, each with 50 perfectly separable pair predictions. Pair-level bootstrap
    # would compute AUC on n=100, but the effective independent sample size is 2 bacteria.
    rows = []
    for bac in ["bac_A", "bac_B"]:
        for i in range(50):
            rows.append(
                {
                    "bacteria": bac,
                    "phage": f"p{i}",
                    "pair_id": f"{bac}__p{i}",
                    "label_row_binary": 1 if i < 25 else 0,
                    "predicted_probability": 0.9 if i < 25 else 0.1,
                }
            )
    cis = _bootstrap_by_unit(rows, unit_key="bacteria", bootstrap_samples=100, bootstrap_random_state=1)
    assert "holdout_roc_auc" in cis
    assert "holdout_brier_score" in cis
    assert cis["holdout_roc_auc"].point_estimate == pytest.approx(1.0)
    # 100 resamples requested, though some may be skipped if they resample one bacterium twice
    # (degenerate single-class). bootstrap_samples_used ≤ requested.
    assert cis["holdout_roc_auc"].bootstrap_samples_requested == 100
    assert cis["holdout_roc_auc"].bootstrap_samples_used <= 100


def test_bootstrap_by_unit_skips_degenerate_single_class_resamples() -> None:
    """If a resample draws only one bacterium and all its labels are identical, AUC is
    undefined. The bootstrap should skip that resample, not crash or fabricate an AUC."""
    rows = [
        {
            "bacteria": "all_positive",
            "phage": f"p{i}",
            "pair_id": f"ap__p{i}",
            "label_row_binary": 1,
            "predicted_probability": 0.9,
        }
        for i in range(5)
    ] + [
        {
            "bacteria": "all_negative",
            "phage": f"p{i}",
            "pair_id": f"an__p{i}",
            "label_row_binary": 0,
            "predicted_probability": 0.1,
        }
        for i in range(5)
    ]
    # Bootstrap samples n_bacteria=2 with replacement. ~half of resamples pick the same
    # bacterium twice → all one class → AUC skipped. The function should complete.
    cis = _bootstrap_by_unit(rows, unit_key="bacteria", bootstrap_samples=200, bootstrap_random_state=42)
    assert isinstance(cis["holdout_roc_auc"], BootstrapMetricCI)
    # Some resamples should have been skipped.
    assert cis["holdout_roc_auc"].bootstrap_samples_used < 200


def test_run_ch05_eval_rejects_phage_slot_dir_with_missing_slots(tmp_path: Path) -> None:
    """A typoed --phage-slot-dir must hard-raise, not silently fall back to baseline slots.

    Without this guard, sx03_eval.patch_context_with_extended_slots logs a warning and
    keeps the original slot, which would produce a self-consistent but wrong result
    (training under a different phage feature set than the caller intended).
    """
    empty_slots = tmp_path / "nonexistent_slots"
    empty_slots.mkdir()
    # Create only one of the three expected slot subdirs to prove we check all three.
    (empty_slots / "phage_projection").mkdir()
    (empty_slots / "phage_projection" / "features.csv").write_text("phage\n")

    with pytest.raises(FileNotFoundError, match="phage_stats.*phage_rbp_struct"):
        run_ch05_eval(
            device_type="cpu",
            output_dir=tmp_path / "out",
            cache_dir=tmp_path / "cache",  # Never reached — guard fires before cache load.
            candidate_dir=tmp_path / "cand",
            phage_slot_dir=empty_slots,
        )
