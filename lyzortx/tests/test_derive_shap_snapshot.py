"""Unit tests for lyzortx/pipeline/autoresearch/derive_shap_snapshot.py.

The end-to-end SHAP pipeline is too slow for CI (30-60 min wallclock) and requires
gitignored CH05 artifacts, so this suite focuses on deterministic helpers: the SHAP
aggregator that collapses per-row contributions into per-pair max-concentration
contributions, and the feature-value aggregator. Both are reusable by EX04 follow-ups.
"""

from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from lyzortx.pipeline.autoresearch.derive_shap_snapshot import (
    _aggregate_fold_feature_values_to_pair_max_conc,
    _aggregate_fold_shap_to_pair_max_conc,
)


def _mock_pair_source() -> pd.Series:
    return pd.Series(
        {
            "001__phageA": "guelin",
            "001__phageB": "guelin",
            "002__phageA": "basel",
        },
        name="source",
    )


def test_shap_aggregator_picks_max_concentration_row_per_pair() -> None:
    """A pair observed at three concentrations should contribute its max-conc SHAP only."""
    holdout_design = pd.DataFrame(
        {
            "pair_id": ["001__phageA", "001__phageA", "001__phageA", "001__phageB"],
            "bacteria": ["001", "001", "001", "001"],
            "phage": ["phageA", "phageA", "phageA", "phageB"],
            "log10_pfu_ml": [4.7, 7.7, 8.7, 8.7],
            "feat_x": [0.1, 0.2, 0.3, 0.4],
            "feat_y": [1.0, 2.0, 3.0, 4.0],
        }
    )
    shap_matrix = np.array(
        [
            [0.11, 0.22],  # 4.7
            [0.12, 0.24],  # 7.7
            [0.13, 0.26],  # 8.7 -> selected for phageA
            [0.90, 0.80],  # 8.7 -> selected for phageB
        ]
    )
    rfe_features = ["feat_x", "feat_y"]
    agg = _aggregate_fold_shap_to_pair_max_conc(
        holdout_design=holdout_design,
        shap_matrix=shap_matrix,
        rfe_features=rfe_features,
        pair_source=_mock_pair_source(),
    )
    assert set(agg["pair_id"]) == {"001__phageA", "001__phageB"}
    phage_a = agg[agg["pair_id"] == "001__phageA"].iloc[0]
    assert phage_a["shap__feat_x"] == pytest.approx(0.13)
    assert phage_a["shap__feat_y"] == pytest.approx(0.26)
    assert phage_a["source"] == "guelin"
    phage_b = agg[agg["pair_id"] == "001__phageB"].iloc[0]
    assert phage_b["shap__feat_x"] == pytest.approx(0.90)
    assert phage_b["source"] == "guelin"


def test_shap_aggregator_averages_replicates_at_same_max_conc() -> None:
    """If a pair has multiple replicates at max_conc, mean SHAP across replicates."""
    holdout_design = pd.DataFrame(
        {
            "pair_id": ["001__phageA", "001__phageA", "001__phageA"],
            "bacteria": ["001"] * 3,
            "phage": ["phageA"] * 3,
            "log10_pfu_ml": [8.7, 8.7, 8.7],  # three replicates at same max
            "feat_x": [0.1, 0.2, 0.3],
        }
    )
    shap_matrix = np.array([[0.10], [0.20], [0.30]])
    agg = _aggregate_fold_shap_to_pair_max_conc(
        holdout_design=holdout_design,
        shap_matrix=shap_matrix,
        rfe_features=["feat_x"],
        pair_source=_mock_pair_source(),
    )
    assert len(agg) == 1
    assert agg.iloc[0]["shap__feat_x"] == pytest.approx(0.20)


def test_feature_value_aggregator_carries_raw_features_for_max_conc_pair() -> None:
    holdout_design = pd.DataFrame(
        {
            "pair_id": ["001__phageA", "001__phageA", "002__phageA"],
            "bacteria": ["001", "001", "002"],
            "phage": ["phageA", "phageA", "phageA"],
            "log10_pfu_ml": [4.7, 8.7, 8.7],
            "feat_x": [0.1, 0.2, 0.4],
            "feat_y": [1.0, 2.0, 4.0],
        }
    )
    agg = _aggregate_fold_feature_values_to_pair_max_conc(
        holdout_design=holdout_design,
        rfe_features=["feat_x", "feat_y"],
        pair_source=_mock_pair_source(),
    )
    assert set(agg["pair_id"]) == {"001__phageA", "002__phageA"}
    pa = agg[agg["pair_id"] == "001__phageA"].iloc[0]
    assert pa["feat_x"] == pytest.approx(0.2)  # max-conc row's raw value
    assert pa["feat_y"] == pytest.approx(2.0)
    assert pa["source"] == "guelin"
    pb = agg[agg["pair_id"] == "002__phageA"].iloc[0]
    assert pb["source"] == "basel"
