"""Tests for lyzortx/explainability_ui/build_snapshot.py.

Unit-tests the snapshot extractor's parsers and assemblers against real CH05/CH09
artifact shapes (either mock fixtures or, when available, the live artifacts under
`lyzortx/generated_outputs/`). The end-to-end `write_snapshot` roundtrip is covered
by a gated integration test that skips if the artifacts are missing (they're
gitignored and expensive to regenerate).
"""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any

import pandas as pd
import pytest

from lyzortx.explainability_ui.build_snapshot import (
    _extract_ece_from_report,
    _feature_to_slot,
    _parse_variant_name,
    build_cross_source_json,
    build_feature_importance_json,
    build_reliability_json,
    build_slot_manifest_json,
    build_summary_json,
    write_snapshot,
)

REPO_ROOT = Path(__file__).resolve().parents[2]
CH05_DIR = REPO_ROOT / "lyzortx/generated_outputs/ch05_unified_kfold"
CH09_DIR = REPO_ROOT / "lyzortx/generated_outputs/ch09_calibration_layer"


# ---- Parser unit tests (deterministic, no fixture dependency) ----


def test_feature_to_slot_handles_all_known_prefixes() -> None:
    cases = {
        "host_defense__AbiD": "host_defense",
        "host_surface__host_capsule_profile_cluster_19_score": "host_surface",
        "host_typing__host_serotype": "host_typing",
        "host_stats__host_gc_content": "host_stats",
        "phage_projection__recep_frac_GluI": "phage_projection",
        "phage_stats__phage_gc_content": "phage_stats",
        "phage_kmer__AAAA": "phage_kmer",
        "pair_depo_capsule__cluster_1_x_capsule_type_K1": "pair_depo_capsule",
        "pair_receptor_omp__OmpC_x_OmpC_score": "pair_receptor_omp",
        "pair_concentration__log10_pfu_ml": "pair_concentration",
        "mystery_feature__x": "other",
    }
    for feature, expected_slot in cases.items():
        assert _feature_to_slot(feature) == expected_slot, feature


def test_parse_variant_name_covers_all_ch09_shapes() -> None:
    """CH09 emits 8 variants: {bacteria,phage}_axis × {guelin,basel}_{raw,calibrated-variant}."""
    assert _parse_variant_name("bacteria_axis_guelin_raw") == ("bacteria", "guelin", "raw")
    assert _parse_variant_name("bacteria_axis_guelin_loof_calibrated") == (
        "bacteria",
        "guelin",
        "loof_calibrated",
    )
    assert _parse_variant_name("bacteria_axis_basel_raw") == ("bacteria", "basel", "raw")
    assert _parse_variant_name("bacteria_axis_basel_guelin_calibrator_applied") == (
        "bacteria",
        "basel",
        "guelin_calibrator_applied",
    )
    assert _parse_variant_name("phage_axis_guelin_raw") == ("phage", "guelin", "raw")
    assert _parse_variant_name("phage_axis_basel_guelin_calibrator_applied") == (
        "phage",
        "basel",
        "guelin_calibrator_applied",
    )


def test_parse_variant_name_rejects_unknown_axis() -> None:
    with pytest.raises(ValueError, match="Unrecognized CH09 variant name"):
        _parse_variant_name("some_other_axis_guelin_raw")


def test_parse_variant_name_rejects_unknown_source() -> None:
    with pytest.raises(ValueError, match="Unrecognized CH09 variant source"):
        _parse_variant_name("bacteria_axis_mysterious_raw")


def test_extract_ece_from_report_flattens_axis_source_to_variant_name() -> None:
    report = {
        "axis_source_metrics": {
            "bacteria_axis": {
                "guelin_raw": {"auc": 0.81, "brier": 0.17, "ece": 0.12, "n_pairs": 35403},
                "basel_raw": {"auc": 0.73, "brier": 0.21, "ece": 0.20, "n_pairs": 1240},
            },
            "phage_axis": {
                "guelin_raw": {"auc": 0.88, "brier": 0.13, "ece": 0.11, "n_pairs": 35403},
            },
        }
    }
    result = _extract_ece_from_report(report)
    assert result["bacteria_axis_guelin_raw"] == pytest.approx(0.12)
    assert result["bacteria_axis_basel_raw"] == pytest.approx(0.20)
    assert result["phage_axis_guelin_raw"] == pytest.approx(0.11)


# ---- Assembler unit tests with inline mock data ----


def _mock_combined_summary() -> dict[str, Any]:
    return {
        "task_id": "CH05",
        "scorecard": "AUC + Brier (nDCG/mAP/top-k retired)",
        "per_phage_blending": "retired (see per-phage-retired-under-chisel)",
        "bacteria_axis": {
            "n_pairs": 36643,
            "n_bacteria": 369,
            "n_phages": 148,
            "aggregate": {
                "holdout_roc_auc": {
                    "point_estimate": 0.807921,
                    "ci_low": 0.793432,
                    "ci_high": 0.822287,
                },
                "holdout_brier_score": {
                    "point_estimate": 0.176269,
                    "ci_low": 0.168760,
                    "ci_high": 0.183997,
                },
            },
        },
        "phage_axis": {
            "n_pairs": 36643,
            "n_bacteria": 369,
            "n_phages": 148,
            "aggregate": {
                "holdout_roc_auc": {
                    "point_estimate": 0.887042,
                    "ci_low": 0.865808,
                    "ci_high": 0.905509,
                },
                "holdout_brier_score": {
                    "point_estimate": 0.135156,
                    "ci_low": 0.122669,
                    "ci_high": 0.148942,
                },
            },
            "cross_source": [
                {
                    "source": "guelin",
                    "n_pairs": 35403,
                    "n_phages": 96,
                    "auc_point": 0.887091,
                    "auc_low": 0.867058,
                    "auc_high": 0.907535,
                    "brier_point": 0.134624,
                    "brier_low": 0.121341,
                    "brier_high": 0.148332,
                },
                {
                    "source": "basel",
                    "n_pairs": 1240,
                    "n_phages": 52,
                    "auc_point": 0.895159,
                    "auc_low": 0.838062,
                    "auc_high": 0.942456,
                    "brier_point": 0.150337,
                    "brier_low": 0.124741,
                    "brier_high": 0.177961,
                },
            ],
        },
        "phage_axis_generalization_gap_auc": -0.079121,
        "cross_source_auc_delta": 0.008068,
        "ch04_baseline_auc": 0.808276,
        "elapsed_seconds": 3224.4,
    }


def _mock_ch09_report() -> dict[str, Any]:
    return {
        "task_id": "CH09",
        "axis_source_metrics": {
            "bacteria_axis": {
                "guelin_raw": {"auc": 0.81036, "brier": 0.17484, "ece": 0.12195, "n_pairs": 35403},
                "basel_raw": {"auc": 0.73921, "brier": 0.21708, "ece": 0.20466, "n_pairs": 1240},
            },
            "phage_axis": {
                "guelin_raw": {"auc": 0.88709, "brier": 0.13462, "ece": 0.11139, "n_pairs": 35403},
                "basel_raw": {"auc": 0.89516, "brier": 0.15034, "ece": 0.18794, "n_pairs": 1240},
            },
        },
    }


def test_build_summary_exposes_headline_numbers_at_4dp() -> None:
    summary = build_summary_json(_mock_combined_summary())
    assert summary["task_id"] == "CH05"
    assert summary["axes"]["bacteria"]["auc"]["point"] == pytest.approx(0.807921)
    assert summary["axes"]["bacteria"]["auc"]["lo"] == pytest.approx(0.793432)
    assert summary["axes"]["bacteria"]["auc"]["hi"] == pytest.approx(0.822287)
    assert summary["axes"]["phage"]["auc"]["point"] == pytest.approx(0.887042)
    assert summary["axes"]["phage"]["brier"]["point"] == pytest.approx(0.135156)
    assert summary["axes"]["bacteria"]["n_pairs"] == 36643
    assert summary["axes"]["phage"]["n_phages"] == 148


def test_build_cross_source_carries_phage_axis_cis_and_bact_axis_point_only() -> None:
    result = build_cross_source_json(_mock_combined_summary(), _mock_ch09_report())
    assert len(result["phage"]) == 2
    basel_phage = next(r for r in result["phage"] if r["source"] == "basel")
    assert basel_phage["auc"]["point"] == pytest.approx(0.895159)
    assert basel_phage["auc"]["lo"] == pytest.approx(0.838062)
    # Bacteria-axis has CH09 point estimates only; CIs stay null (CH05 doesn't bootstrap
    # bact-axis cross-source, so presenting CIs would be fabrication).
    assert len(result["bacteria"]) == 2
    basel_bact = next(r for r in result["bacteria"] if r["source"] == "basel")
    assert basel_bact["auc"]["point"] == pytest.approx(0.73921)
    assert basel_bact["auc"]["lo"] is None
    assert basel_bact["auc"]["hi"] is None
    assert basel_bact["ece"] == pytest.approx(0.20466)


def test_build_feature_importance_attaches_slot_tag(tmp_path: Path) -> None:
    """FI CSVs get a `slot` column for UI coloring; prefix mapping must match."""
    bact_csv = tmp_path / "ch05_bacteria_axis_feature_importance.csv"
    phage_csv = tmp_path / "ch05_phage_axis_feature_importance.csv"
    pd.DataFrame(
        {
            "feature": [
                "host_typing__host_serotype",
                "phage_stats__phage_gc_content",
                "pair_concentration__log10_pfu_ml",
            ],
            "is_concentration_feature": [False, False, True],
            "mean_importance": [2500.0, 570.0, 320.0],
            "n_folds_selected": [30, 30, 30],
        }
    ).to_csv(bact_csv, index=False)
    pd.DataFrame(
        {
            "feature": ["host_typing__host_o_type"],
            "is_concentration_feature": [False],
            "mean_importance": [300.0],
            "n_folds_selected": [30],
        }
    ).to_csv(phage_csv, index=False)
    fi = build_feature_importance_json(tmp_path)
    assert fi["bacteria_axis"][0]["slot"] == "host_typing"
    assert fi["bacteria_axis"][1]["slot"] == "phage_stats"
    assert fi["bacteria_axis"][2]["slot"] == "pair_concentration"
    assert fi["phage_axis"][0]["slot"] == "host_typing"


def test_build_reliability_groups_by_variant_and_attaches_ece(tmp_path: Path) -> None:
    tables = pd.DataFrame(
        {
            "variant": [
                "bacteria_axis_guelin_raw",
                "bacteria_axis_guelin_raw",
                "bacteria_axis_guelin_loof_calibrated",
                "phage_axis_basel_raw",
            ],
            "bin_lo": [0.0, 0.1, 0.0, 0.0],
            "bin_hi": [0.1, 0.2, 0.1, 0.1],
            "n_pairs": [9084, 5042, 9000, 200],
            "mean_predicted": [0.04, 0.15, 0.05, 0.06],
            "observed_rate": [0.05, 0.12, 0.05, 0.03],
            "reliability_gap": [0.01, -0.03, 0.0, -0.03],
        }
    )
    tables.to_csv(tmp_path / "ch09_reliability_tables.csv", index=False)
    report = {
        "axis_source_metrics": {
            "bacteria_axis": {
                "guelin_raw": {"auc": 0.81, "brier": 0.17, "ece": 0.122, "n_pairs": 35403},
                "guelin_loof_calibrated": {
                    "auc": 0.805,
                    "brier": 0.15,
                    "ece": 0.006,
                    "n_pairs": 35403,
                },
            },
            "phage_axis": {
                "basel_raw": {"auc": 0.90, "brier": 0.15, "ece": 0.188, "n_pairs": 1240},
            },
        }
    }
    (tmp_path / "ch09_calibration_report.json").write_text(json.dumps(report))
    result = build_reliability_json(tmp_path)
    assert set(result["variants"].keys()) == {
        "bacteria_axis_guelin_raw",
        "bacteria_axis_guelin_loof_calibrated",
        "phage_axis_basel_raw",
    }
    raw = result["variants"]["bacteria_axis_guelin_raw"]
    assert raw["axis"] == "bacteria"
    assert raw["source"] == "guelin"
    assert raw["variant"] == "raw"
    assert len(raw["bins"]) == 2
    assert raw["ece"] == pytest.approx(0.122)
    assert result["variants"]["bacteria_axis_guelin_loof_calibrated"]["ece"] == pytest.approx(0.006)


def test_build_slot_manifest_counts_features_per_axis() -> None:
    fi = {
        "bacteria_axis": [
            {
                "feature": "host_typing__a",
                "slot": "host_typing",
                "is_concentration_feature": False,
                "mean_importance": 1000,
                "n_folds_selected": 30,
            },
            {
                "feature": "host_typing__b",
                "slot": "host_typing",
                "is_concentration_feature": False,
                "mean_importance": 500,
                "n_folds_selected": 30,
            },
            {
                "feature": "phage_stats__gc",
                "slot": "phage_stats",
                "is_concentration_feature": False,
                "mean_importance": 300,
                "n_folds_selected": 30,
            },
        ],
        "phage_axis": [
            {
                "feature": "host_typing__a",
                "slot": "host_typing",
                "is_concentration_feature": False,
                "mean_importance": 1100,
                "n_folds_selected": 30,
            },
        ],
    }
    manifest = build_slot_manifest_json(fi)
    slots = {s["slot_name"]: s for s in manifest["slots"]}
    host_typing = slots["host_typing"]
    assert host_typing["axes"]["bacteria_axis"]["feature_count"] == 2
    assert host_typing["axes"]["bacteria_axis"]["cumulative_importance"] == pytest.approx(1500.0)
    assert host_typing["axes"]["phage_axis"]["feature_count"] == 1
    # Slots with no features on an axis still appear (with count=0) so UI can show an empty card.
    phage_stats = slots["phage_stats"]
    assert phage_stats["axes"]["phage_axis"]["feature_count"] == 0


# ---- End-to-end integration test (skip if artifacts missing) ----


_SKIP_REASON = (
    "CH05/CH09 artifacts are gitignored; regenerate with "
    "`python -m lyzortx.pipeline.autoresearch.ch05_eval` to enable this test."
)


@pytest.mark.skipif(
    not (CH05_DIR / "ch05_combined_summary.json").exists() or not (CH09_DIR / "ch09_calibration_report.json").exists(),
    reason=_SKIP_REASON,
)
def test_write_snapshot_emits_all_seven_jsons_with_expected_headline_numbers(tmp_path: Path) -> None:
    """End-to-end: extractor → 7 JSONs, bact-axis AUC 0.807921, phage-axis AUC 0.887042."""
    out_dir = tmp_path / "snapshot"
    written = write_snapshot(out_dir=out_dir, ch05_dir=CH05_DIR, ch09_dir=CH09_DIR)
    expected_names = {
        "ch05_summary.json",
        "ch05_cross_source.json",
        "ch05_feature_importance.json",
        "ch05_predictions_bact.json",
        "ch05_predictions_phage.json",
        "ch05_reliability.json",
        "slot_manifest.json",
    }
    assert set(written.keys()) == expected_names
    with (out_dir / "ch05_summary.json").open() as f:
        summary = json.load(f)
    assert summary["axes"]["bacteria"]["auc"]["point"] == pytest.approx(0.807921, abs=1e-5)
    assert summary["axes"]["phage"]["auc"]["point"] == pytest.approx(0.887042, abs=1e-5)
    assert summary["axes"]["bacteria"]["n_pairs"] == 36643
    with (out_dir / "ch05_feature_importance.json").open() as f:
        fi = json.load(f)
    assert len(fi["bacteria_axis"]) > 200
    assert len(fi["phage_axis"]) > 200
    assert fi["bacteria_axis"][0]["feature"] == "host_typing__host_serotype"
    assert fi["bacteria_axis"][0]["slot"] == "host_typing"
