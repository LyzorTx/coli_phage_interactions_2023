"""Tests for CH08 paired bootstrap and slot attachment plumbing."""

from __future__ import annotations

import numpy as np

from lyzortx.pipeline.autoresearch.ch08_wave2_reaudit import _bootstrap_paired_by_bacterium


def _make_pair_rows(n_pairs: int, seed: int, signal_boost: float = 0.0) -> list[dict]:
    rng = np.random.default_rng(seed)
    labels = rng.integers(0, 2, size=n_pairs)
    base = rng.uniform(0.1, 0.5, size=n_pairs)
    preds = np.clip(base + signal_boost * labels, 0.0, 1.0)
    return [
        {
            "pair_id": f"b{i // 5}__p{i}",
            "bacteria": f"b{i // 5}",
            "phage": f"p{i}",
            "label_row_binary": int(labels[i]),
            "predicted_probability": float(preds[i]),
        }
        for i in range(n_pairs)
    ]


def test_paired_bootstrap_returns_baseline_variant_and_delta() -> None:
    baseline = _make_pair_rows(200, seed=1, signal_boost=0.2)
    variant = _make_pair_rows(200, seed=1, signal_boost=0.35)  # stronger signal
    report = _bootstrap_paired_by_bacterium(baseline, variant, bootstrap_samples=50)
    assert report["n_common_pairs"] == 200
    assert report["baseline"]["auc_point"] is not None
    assert report["variant"]["auc_point"] is not None
    # Variant boost should give variant > baseline on point AUC.
    assert report["variant"]["auc_point"] >= report["baseline"]["auc_point"]
    # Delta has CI around the point estimate.
    assert report["delta"]["auc_ci"]["low"] is not None
    assert report["delta"]["auc_ci"]["high"] is not None
    assert report["delta"]["auc_ci"]["low"] <= report["delta"]["auc_point"]
    assert report["delta"]["auc_point"] <= report["delta"]["auc_ci"]["high"]


def test_paired_bootstrap_null_case_delta_centered_near_zero() -> None:
    baseline = _make_pair_rows(300, seed=7, signal_boost=0.25)
    # Identical variant → delta should center at 0 exactly.
    variant = [dict(r) for r in baseline]
    report = _bootstrap_paired_by_bacterium(baseline, variant, bootstrap_samples=100)
    assert abs(report["delta"]["auc_point"]) < 1e-9
    assert abs(report["delta"]["brier_point"]) < 1e-9
    assert abs(report["delta"]["auc_ci"]["low"] or 0.0) < 1e-6
    assert abs(report["delta"]["auc_ci"]["high"] or 0.0) < 1e-6


def test_paired_bootstrap_on_nonoverlapping_pairs_uses_intersection() -> None:
    baseline = _make_pair_rows(100, seed=2, signal_boost=0.2)
    variant = _make_pair_rows(100, seed=2, signal_boost=0.2)
    # Drop half of variant pairs; bootstrap must restrict to common pairs.
    variant = variant[:50]
    report = _bootstrap_paired_by_bacterium(baseline, variant, bootstrap_samples=25)
    assert report["n_common_pairs"] == 50
