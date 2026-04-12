"""Tests for SPANDEX evaluation metrics (nDCG, mAP)."""

import pytest

from lyzortx.pipeline.autoresearch.spandex_metrics import (
    compute_binary_metrics,
    compute_per_bacterium_map,
    compute_per_bacterium_ndcg,
    evaluate_holdout_rows,
)


def _make_rows(bacteria: str, phages_relevance_pred: list[tuple[str, float, float]]) -> list[dict]:
    """Helper: create rows for one bacterium with (phage, relevance, predicted_prob)."""
    return [
        {
            "bacteria": bacteria,
            "phage": phage,
            "mlc_score": rel,
            "label_binary": int(rel > 0) if rel is not None else None,
            "predicted_probability": pred,
        }
        for phage, rel, pred in phages_relevance_pred
    ]


class TestNDCG:
    def test_perfect_ranking(self):
        """Model ranks phages in exact MLC order -> nDCG = 1.0."""
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.95),
                ("P2", 3.0, 0.80),
                ("P3", 1.0, 0.60),
                ("P4", 0.0, 0.10),
            ],
        )
        ndcg = compute_per_bacterium_ndcg(rows)
        assert ndcg is not None
        assert ndcg == pytest.approx(1.0, abs=0.001)

    def test_inverted_ranking(self):
        """Model ranks lowest-MLC phage highest -> nDCG < 1.0."""
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.10),
                ("P2", 1.0, 0.95),
                ("P3", 0.0, 0.50),
            ],
        )
        ndcg = compute_per_bacterium_ndcg(rows)
        assert ndcg is not None
        assert ndcg < 1.0

    def test_all_negatives_returns_none(self):
        """No positives -> nDCG undefined."""
        rows = _make_rows("B1", [("P1", 0.0, 0.5), ("P2", 0.0, 0.3)])
        assert compute_per_bacterium_ndcg(rows) is None

    def test_partial_ground_truth(self):
        """Rows with None relevance are excluded."""
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.95),
                ("P2", 0.0, 0.10),
            ],
        ) + [{"bacteria": "B1", "phage": "P3", "mlc_score": None, "predicted_probability": 0.80}]
        ndcg = compute_per_bacterium_ndcg(rows)
        assert ndcg is not None
        assert ndcg == pytest.approx(1.0, abs=0.001)

    def test_multi_bacterium_average(self):
        """nDCG is averaged across bacteria."""
        rows = (
            _make_rows("B1", [("P1", 4.0, 0.9), ("P2", 0.0, 0.1)])  # perfect
            + _make_rows("B2", [("P1", 4.0, 0.1), ("P2", 0.0, 0.9)])  # inverted
        )
        ndcg = compute_per_bacterium_ndcg(rows)
        assert ndcg is not None
        assert 0.5 < ndcg < 1.0  # average of 1.0 and inverted (nDCG floor is ~0.63 for 2 items)


class TestMAP:
    def test_perfect_ranking(self):
        """Positives ranked above all negatives -> mAP = 1.0."""
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.95),
                ("P2", 2.0, 0.80),
                ("P3", 0.0, 0.10),
                ("P4", 0.0, 0.05),
            ],
        )
        mAP = compute_per_bacterium_map(rows)
        assert mAP is not None
        assert mAP == pytest.approx(1.0, abs=0.001)

    def test_all_same_class_returns_none(self):
        """All positives or all negatives -> mAP undefined."""
        rows = _make_rows("B1", [("P1", 4.0, 0.9), ("P2", 2.0, 0.5)])
        assert compute_per_bacterium_map(rows) is None

    def test_partial_ground_truth(self):
        """Rows with None label are excluded."""
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.95),
                ("P2", 0.0, 0.10),
            ],
        ) + [{"bacteria": "B1", "phage": "P3", "label_binary": None, "predicted_probability": 0.80}]
        mAP = compute_per_bacterium_map(rows)
        assert mAP is not None
        assert mAP == pytest.approx(1.0, abs=0.001)


class TestBinaryMetrics:
    def test_basic(self):
        rows = _make_rows(
            "B1",
            [
                ("P1", 4.0, 0.95),
                ("P2", 0.0, 0.10),
            ],
        )
        metrics = compute_binary_metrics(rows)
        assert metrics["roc_auc"] == pytest.approx(1.0)
        assert metrics["brier_score"] is not None
        assert metrics["brier_score"] < 0.1


class TestEvaluateHoldoutRows:
    def test_returns_all_metrics(self):
        rows = _make_rows("B1", [("P1", 4.0, 0.9), ("P2", 0.0, 0.1)]) + _make_rows(
            "B2", [("P1", 3.0, 0.8), ("P2", 0.0, 0.2)]
        )
        metrics = evaluate_holdout_rows(rows)
        assert "holdout_ndcg" in metrics
        assert "holdout_map" in metrics
        assert "holdout_roc_auc" in metrics
        assert "holdout_brier_score" in metrics
        assert all(v is not None for v in metrics.values())
