"""Tests for CH04 per-row training frame and pair-max-concentration selection."""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lyzortx.pipeline.autoresearch.ch03_row_expansion import GUELIN_NEAT_LOG10_PFU_ML
from lyzortx.pipeline.autoresearch.ch04_eval import (
    CONCENTRATION_FEATURE_COLUMN,
    build_clean_row_training_frame,
    compute_aggregate_auc_brier,
    select_pair_max_concentration_rows,
)


def _raw_row(pair_id: str, score: str, log_dilution: int, replicate: int = 1) -> dict:
    bacteria, phage = pair_id.split("__")
    return {
        "bacteria": bacteria,
        "phage": phage,
        "pair_id": pair_id,
        "log_dilution": log_dilution,
        "log10_pfu_ml": GUELIN_NEAT_LOG10_PFU_ML + float(log_dilution),
        "replicate": replicate,
        "score": score,
        "split_holdout": "train_non_holdout",
        "is_hard_trainable": "1",
    }


def test_build_clean_row_training_frame_drops_n_rows() -> None:
    rows = pd.DataFrame(
        [
            _raw_row("bac_A__phage_X", "0", 0),
            _raw_row("bac_A__phage_X", "1", 0, replicate=2),
            _raw_row("bac_A__phage_X", "n", -1),
            _raw_row("bac_B__phage_X", "n", 0),
        ]
    )
    clean = build_clean_row_training_frame(rows)
    assert len(clean) == 2
    assert set(clean["score"]) == {"0", "1"}
    assert set(clean["label_row_binary"]) == {0, 1}
    assert CONCENTRATION_FEATURE_COLUMN in clean.columns
    assert (clean[CONCENTRATION_FEATURE_COLUMN] == clean["log10_pfu_ml"].astype(float)).all()
    # Guelin encoding: neat (log_dilution=0) → 8.7; -1 → 7.7 under GUELIN_NEAT_LOG10_PFU_ML + log_dilution.
    assert (clean[CONCENTRATION_FEATURE_COLUMN] >= 4.0).all()
    assert (clean[CONCENTRATION_FEATURE_COLUMN] <= 9.0).all()


def test_build_clean_row_training_frame_preserves_all_clean_rows() -> None:
    rows = pd.DataFrame([_raw_row("bac_A__phage_X", "0", d, replicate=r) for d in [0, -1, -2, -4] for r in [1, 2]])
    clean = build_clean_row_training_frame(rows)
    assert len(clean) == len(rows)


def _pred_row(pair_id: str, log_dilution: int, replicate: int, label: int, proba: float) -> dict:
    bacteria, phage = pair_id.split("__")
    return {
        "pair_id": pair_id,
        "bacteria": bacteria,
        "phage": phage,
        "log_dilution": log_dilution,
        "log10_pfu_ml": GUELIN_NEAT_LOG10_PFU_ML + float(log_dilution),
        "replicate": replicate,
        "label_row_binary": label,
        "predicted_probability": proba,
    }


def test_select_pair_max_concentration_keeps_highest_log_dilution_per_pair() -> None:
    # For a pair with observations at log_dilution ∈ {0, -1, -2}, log_dilution=0 (log10_pfu_ml=8.7)
    # is the highest actual concentration — the evaluator should pick that row.
    predictions = pd.DataFrame(
        [
            _pred_row("bac_A__phage_X", 0, 1, 1, 0.9),
            _pred_row("bac_A__phage_X", -1, 1, 0, 0.4),
            _pred_row("bac_A__phage_X", -2, 1, 0, 0.1),
            _pred_row("bac_B__phage_X", -2, 1, 0, 0.2),
            _pred_row("bac_B__phage_X", -4, 1, 0, 0.05),
        ]
    )
    pair_rows = select_pair_max_concentration_rows(predictions)
    assert len(pair_rows) == 2
    pair_a = pair_rows[pair_rows["pair_id"] == "bac_A__phage_X"].iloc[0]
    pair_b = pair_rows[pair_rows["pair_id"] == "bac_B__phage_X"].iloc[0]
    assert pair_a["log10_pfu_ml"] == pytest.approx(GUELIN_NEAT_LOG10_PFU_ML)
    assert pair_a["predicted_probability"] == pytest.approx(0.9)
    assert pair_a["label_row_binary"] == 1
    # bac_B's highest observed concentration is log_dilution=-2 (not every pair is tested at 0).
    assert pair_b["log10_pfu_ml"] == pytest.approx(GUELIN_NEAT_LOG10_PFU_ML - 2.0)
    assert pair_b["predicted_probability"] == pytest.approx(0.2)


def test_select_pair_max_concentration_averages_replicates_at_top() -> None:
    # Three replicates at log_dilution=0 — averaged; any replicate with label=1 wins the label.
    predictions = pd.DataFrame(
        [_pred_row("bac_A__phage_X", 0, r, label, p) for r, label, p in [(1, 0, 0.6), (2, 1, 0.8), (3, 0, 0.7)]]
    )
    pair_rows = select_pair_max_concentration_rows(predictions)
    assert len(pair_rows) == 1
    row = pair_rows.iloc[0]
    assert row["predicted_probability"] == pytest.approx((0.6 + 0.8 + 0.7) / 3)
    assert row["label_row_binary"] == 1
    assert row["n_replicates_at_max"] == 3


def test_compute_aggregate_auc_brier() -> None:
    rows = [
        {"label_row_binary": 1, "predicted_probability": 0.9, "bacteria": "A"},
        {"label_row_binary": 0, "predicted_probability": 0.2, "bacteria": "A"},
        {"label_row_binary": 1, "predicted_probability": 0.85, "bacteria": "B"},
        {"label_row_binary": 0, "predicted_probability": 0.3, "bacteria": "B"},
    ]
    result = compute_aggregate_auc_brier(rows)
    assert result["auc"] == pytest.approx(1.0)  # perfect separation
    assert result["brier"] == pytest.approx(np.mean([(0.9 - 1) ** 2, 0.2**2, (0.85 - 1) ** 2, 0.3**2]))
