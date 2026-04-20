"""Tests for CH07 double cross-validation cell enumeration and pair-level bootstrap."""

from __future__ import annotations

import numpy as np

from lyzortx.pipeline.autoresearch.ch07_both_axis_holdout import (
    _bootstrap_pair_level,
    _cells_plan,
)


def test_cells_plan_enumerates_cartesian_product() -> None:
    plan = _cells_plan(3, 4)
    assert len(plan) == 12
    assert plan[0] == (0, 0)
    assert plan[-1] == (2, 3)
    # Row-major: bacteria-fold outer, phage-fold inner.
    assert plan[4] == (1, 0)


def test_bootstrap_pair_level_produces_finite_cis() -> None:
    rng = np.random.default_rng(42)
    n = 200
    labels = rng.integers(0, 2, size=n)
    preds = rng.uniform(0.1, 0.5, size=n)
    # Inject a signal so AUC > 0.5 (keep preds in [0, 1]).
    preds = np.clip(preds + 0.3 * labels, 0.0, 1.0)
    rows = [
        {
            "pair_id": f"b{i}__p{i}",
            "bacteria": f"b{i}",
            "phage": f"p{i % 10}",
            "label_row_binary": int(labels[i]),
            "predicted_probability": float(preds[i]),
        }
        for i in range(n)
    ]
    cis = _bootstrap_pair_level(rows, bootstrap_samples=100, bootstrap_random_state=1)
    assert cis["holdout_roc_auc"].point_estimate is not None
    assert cis["holdout_roc_auc"].ci_low is not None
    assert cis["holdout_roc_auc"].ci_high is not None
    assert cis["holdout_roc_auc"].ci_low <= cis["holdout_roc_auc"].point_estimate
    assert cis["holdout_roc_auc"].point_estimate <= cis["holdout_roc_auc"].ci_high
    # CI width is positive for a non-degenerate sample.
    assert cis["holdout_roc_auc"].ci_high > cis["holdout_roc_auc"].ci_low


def test_bootstrap_pair_level_handles_single_class_draws() -> None:
    # Rows with mostly-zero labels — some resamples will be all-zero and skipped.
    rng = np.random.default_rng(1)
    n = 60
    labels = np.zeros(n, dtype=int)
    labels[:5] = 1  # only 5 positives
    preds = rng.uniform(0, 1, size=n)
    rows = [
        {
            "pair_id": f"b{i}__p{i}",
            "bacteria": f"b{i}",
            "phage": f"p{i}",
            "label_row_binary": int(labels[i]),
            "predicted_probability": float(preds[i]),
        }
        for i in range(n)
    ]
    cis = _bootstrap_pair_level(rows, bootstrap_samples=50, bootstrap_random_state=2)
    # Some resamples may drop below 2 classes; we expect <= requested.
    assert 0 <= cis["holdout_roc_auc"].bootstrap_samples_used <= 50
    assert cis["holdout_brier_score"].point_estimate is not None
