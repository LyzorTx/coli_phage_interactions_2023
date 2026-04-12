import sys
import types

import numpy as np
import pytest

from lyzortx.pipeline.track_g.steps.train_v1_binary_classifier import (
    FeatureSpace,
    build_feature_space,
    compute_top3_hit_rate,
    fit_final_estimator,
    merge_expanded_feature_rows,
    make_lightgbm_estimator,
    select_best_candidate,
)


def test_merge_expanded_feature_rows_adds_split_phage_and_pair_features() -> None:
    merged = merge_expanded_feature_rows(
        track_c_pair_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "label_hard_any_lysis": "1",
                "host_surface_lps_core_type": "R1",
            }
        ],
        split_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "cv_group": "G1",
                "split_holdout": "train_non_holdout",
                "split_cv5_fold": "0",
                "is_hard_trainable": "1",
            }
        ],
        phage_feature_blocks=[
            [{"phage": "P1", "phage_gc_content": "0.5"}],
            [{"phage": "P1", "phage_viridic_mds_00": "0.2"}],
        ],
        pair_feature_blocks=[
            [{"pair_id": "B1__P1", "bacteria": "B1", "phage": "P1", "lookup_available": "1"}],
            [
                {
                    "pair_id": "B1__P1",
                    "bacteria": "B1",
                    "phage": "P1",
                    "isolation_host_distance": "0.3",
                }
            ],
        ],
    )

    assert merged[0]["cv_group"] == "G1"
    assert merged[0]["phage_gc_content"] == "0.5"
    assert merged[0]["lookup_available"] == "1"
    assert merged[0]["isolation_host_distance"] == "0.3"


def test_merge_expanded_feature_rows_can_zero_fill_missing_pair_features() -> None:
    merged = merge_expanded_feature_rows(
        track_c_pair_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "label_hard_any_lysis": "1",
            },
            {
                "pair_id": "B2__P1",
                "bacteria": "B2",
                "phage": "P1",
                "label_hard_any_lysis": "0",
            },
        ],
        split_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "cv_group": "G1",
                "split_holdout": "train_non_holdout",
                "split_cv5_fold": "0",
                "is_hard_trainable": "1",
            },
            {
                "pair_id": "B2__P1",
                "bacteria": "B2",
                "phage": "P1",
                "cv_group": "G2",
                "split_holdout": "holdout_test",
                "split_cv5_fold": "-1",
                "is_hard_trainable": "1",
            },
        ],
        phage_feature_blocks=[[{"phage": "P1", "phage_gc_content": "0.5"}]],
        pair_feature_blocks=[[{"pair_id": "B1__P1", "bacteria": "B1", "phage": "P1", "lookup_available": "1"}]],
        allow_missing_pair_features=True,
    )

    assert merged[1]["lookup_available"] == 0.0


def test_merge_expanded_feature_rows_still_raises_on_non_holdout_miss_with_zero_fill_enabled() -> None:
    with pytest.raises(KeyError, match="Missing pair-level feature row for pair_id B1__P1"):
        merge_expanded_feature_rows(
            track_c_pair_rows=[
                {
                    "pair_id": "B1__P1",
                    "bacteria": "B1",
                    "phage": "P1",
                    "label_hard_any_lysis": "1",
                }
            ],
            split_rows=[
                {
                    "pair_id": "B1__P1",
                    "bacteria": "B1",
                    "phage": "P1",
                    "cv_group": "G1",
                    "split_holdout": "train_non_holdout",
                    "split_cv5_fold": "0",
                    "is_hard_trainable": "1",
                }
            ],
            phage_feature_blocks=[[{"phage": "P1", "phage_gc_content": "0.5"}]],
            pair_feature_blocks=[[{"pair_id": "B2__P1", "bacteria": "B2", "phage": "P1", "lookup_available": "1"}]],
            allow_missing_pair_features=True,
        )


def test_build_feature_space_keeps_v0_columns_and_adds_track_specific_blocks() -> None:
    feature_space = build_feature_space(
        st02_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "host_pathotype": "pt",
                "host_mouse_killed_10": "0",
            }
        ],
        track_c_pair_rows=[
            {
                "pair_id": "B1__P1",
                "bacteria": "B1",
                "phage": "P1",
                "host_pathotype": "pt",
                "host_mouse_killed_10": "0",
                "host_surface_lps_core_type": "R1",
                "host_phylogeny_umap_00": "0.1",
            }
        ],
        track_d_feature_columns=["phage_gc_content", "phage_viridic_mds_00"],
        track_e_feature_columns=["lookup_available", "isolation_host_distance"],
    )

    assert "host_pathotype" in feature_space.categorical_columns
    assert "host_surface_lps_core_type" in feature_space.categorical_columns
    assert "host_phylogeny_umap_00" in feature_space.numeric_columns
    assert "phage_gc_content" in feature_space.numeric_columns
    assert "lookup_available" in feature_space.numeric_columns


def test_compute_top3_hit_rate_reports_all_and_susceptible_denominators() -> None:
    metrics = compute_top3_hit_rate(
        [
            {"bacteria": "B1", "phage": "P1", "label_hard_any_lysis": "1", "predicted_probability": 0.9},
            {"bacteria": "B1", "phage": "P2", "label_hard_any_lysis": "0", "predicted_probability": 0.8},
            {"bacteria": "B1", "phage": "P3", "label_hard_any_lysis": "0", "predicted_probability": 0.7},
            {"bacteria": "B2", "phage": "P1", "label_hard_any_lysis": "0", "predicted_probability": 0.9},
            {"bacteria": "B2", "phage": "P2", "label_hard_any_lysis": "0", "predicted_probability": 0.8},
            {"bacteria": "B2", "phage": "P3", "label_hard_any_lysis": "0", "predicted_probability": 0.7},
        ],
        probability_key="predicted_probability",
    )

    assert metrics["strain_count"] == 2
    assert metrics["hit_count"] == 1
    assert metrics["top3_hit_rate_all_strains"] == 0.5
    assert metrics["susceptible_strain_count"] == 1
    assert metrics["top3_hit_rate_susceptible_only"] == 1.0


def test_fit_final_estimator_forwards_sample_weights_when_requested() -> None:
    captured: dict[str, object] = {}

    class FakeEstimator:
        def fit(self, X, y, sample_weight=None):
            captured["sample_weight"] = list(sample_weight) if sample_weight is not None else None
            captured["train_rows"] = len(y)
            return self

        def predict_proba(self, X):
            return np.array([[0.25, 0.75]] * X.shape[0])

    rows = [
        {
            "pair_id": "B1__P1",
            "bacteria": "B1",
            "phage": "P1",
            "split_holdout": "train_non_holdout",
            "is_hard_trainable": "1",
            "label_hard_any_lysis": "1",
            "feature_a": "1.0",
            "effective_training_weight": "0.2",
        },
        {
            "pair_id": "B1__P2",
            "bacteria": "B1",
            "phage": "P2",
            "split_holdout": "holdout_test",
            "is_hard_trainable": "1",
            "label_hard_any_lysis": "0",
            "feature_a": "2.0",
        },
    ]
    feature_space = FeatureSpace(
        categorical_columns=(),
        numeric_columns=("feature_a",),
        track_c_additional_columns=(),
        track_d_columns=(),
        track_e_columns=(),
    )

    fit_final_estimator(
        rows,
        feature_space,
        estimator_factory=lambda params, seed_offset: FakeEstimator(),
        params={},
        sample_weight_key="effective_training_weight",
    )

    assert captured["train_rows"] == 1
    assert captured["sample_weight"] == [0.2]


def test_select_best_candidate_prefers_auc_then_top3_then_brier() -> None:
    best = select_best_candidate(
        [
            {
                "params": {"name": "a"},
                "summary": {
                    "mean_roc_auc": 0.81,
                    "mean_top3_hit_rate_all_strains": 0.80,
                    "mean_brier_score": 0.20,
                },
            },
            {
                "params": {"name": "b"},
                "summary": {
                    "mean_roc_auc": 0.81,
                    "mean_top3_hit_rate_all_strains": 0.85,
                    "mean_brier_score": 0.25,
                },
            },
        ]
    )

    assert best["params"]["name"] == "b"


def test_make_lightgbm_estimator_enables_determinism_without_forcing_single_thread(monkeypatch) -> None:
    captured: dict[str, object] = {}

    class FakeLGBMClassifier:
        def __init__(self, **kwargs):
            captured.update(kwargs)

    fake_lightgbm = types.SimpleNamespace(LGBMClassifier=FakeLGBMClassifier)
    monkeypatch.setitem(sys.modules, "lightgbm", fake_lightgbm)

    estimator = make_lightgbm_estimator({"num_leaves": 31}, 2, base_random_state=17)

    assert estimator.__class__.__name__ == "FakeLGBMClassifier"
    assert captured["objective"] == "binary"
    assert captured["class_weight"] == "balanced"
    assert captured["random_state"] == 19
    assert captured["deterministic"] is True
    assert captured["force_col_wise"] is True
    assert "n_jobs" not in captured
