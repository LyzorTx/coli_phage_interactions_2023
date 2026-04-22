#!/usr/bin/env python3
"""Extract CH05 + CH09 artifacts into normalized JSONs for the explainability UI.

The static `index.html` fetches these JSONs (either from a local server during
development or from GitHub Release assets in production, per EX03) and renders six
views: Overview, Cross-source, Calibration, Feature importance, Per-slot breakdown,
Predictions table.

This script is deterministic and idempotent: re-running it with the same source
artifacts produces byte-identical output (apart from the `snapshot_generated_at` UTC
timestamp, which is intentionally the only non-deterministic field).

Usage:
    python -m lyzortx.explainability_ui.build_snapshot \
        --out .scratch/explainability_ui/data
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.runtime_contract import SLOT_SPECS

LOGGER = logging.getLogger(__name__)

DEFAULT_CH05_DIR = Path("lyzortx/generated_outputs/ch05_unified_kfold")
DEFAULT_CH09_DIR = Path("lyzortx/generated_outputs/ch09_calibration_layer")
DEFAULT_OUT_DIR = Path(".scratch/explainability_ui/data")

CONCENTRATION_FEATURE_PREFIX = "pair_concentration__"


def _read_combined_summary(ch05_dir: Path) -> dict[str, Any]:
    path = ch05_dir / "ch05_combined_summary.json"
    if not path.exists():
        raise FileNotFoundError(
            f"Missing CH05 combined summary at {path}. Run "
            "`python -m lyzortx.pipeline.autoresearch.ch05_eval --device-type cpu` first."
        )
    with path.open(encoding="utf-8") as f:
        return json.load(f)


def _ci_block(ci_dict: dict[str, Any]) -> dict[str, float | None]:
    """Normalize a CH05 {point_estimate, ci_low, ci_high} dict into UI-friendly shape."""
    return {
        "point": ci_dict.get("point_estimate"),
        "lo": ci_dict.get("ci_low"),
        "hi": ci_dict.get("ci_high"),
    }


def build_summary_json(combined: dict[str, Any]) -> dict[str, Any]:
    """Shape the headline numbers for the Overview tab."""
    bact = combined["bacteria_axis"]
    phage = combined["phage_axis"]
    return {
        "task_id": combined["task_id"],
        "scorecard": combined["scorecard"],
        "per_phage_blending": combined["per_phage_blending"],
        "frame": (
            "CHISEL unified 148-phage × 369-bacteria panel, per-row binary training, "
            "all-pairs only, Arm 3 Moriniere per-receptor-class k-mer fractions as the "
            "canonical phage_projection slot."
        ),
        "elapsed_seconds": combined["elapsed_seconds"],
        "ch04_baseline_auc": combined.get("ch04_baseline_auc"),
        "cross_source_auc_delta": combined.get("cross_source_auc_delta"),
        "phage_axis_generalization_gap_auc": combined.get("phage_axis_generalization_gap_auc"),
        "axes": {
            "bacteria": {
                "fold_strategy": "10-fold via CH02 cv_group hash",
                "n_pairs": bact["n_pairs"],
                "n_bacteria": bact["n_bacteria"],
                "n_phages": bact["n_phages"],
                "auc": _ci_block(bact["aggregate"]["holdout_roc_auc"]),
                "brier": _ci_block(bact["aggregate"]["holdout_brier_score"]),
            },
            "phage": {
                "fold_strategy": "10-fold StratifiedKFold by ICTV family",
                "n_pairs": phage["n_pairs"],
                "n_bacteria": phage["n_bacteria"],
                "n_phages": phage["n_phages"],
                "auc": _ci_block(phage["aggregate"]["holdout_roc_auc"]),
                "brier": _ci_block(phage["aggregate"]["holdout_brier_score"]),
            },
        },
    }


def build_cross_source_json(combined: dict[str, Any], ch09_report: dict[str, Any]) -> dict[str, Any]:
    """Per-axis Guelin vs BASEL breakdown.

    Phage-axis entries carry bootstrap CIs (from `ch05_phage_axis.cross_source`); bact-axis
    entries are point estimates only (CH05 doesn't bootstrap bact-axis cross-source — the
    paired bootstrap is phage-axis scoped by design). Bact-axis numbers come from CH09's
    `axis_source_metrics.bacteria_axis.{guelin_raw, basel_raw}` which reports point AUC/
    Brier/ECE without resampling.
    """
    phage_rows = []
    for row in combined["phage_axis"]["cross_source"]:
        phage_rows.append(
            {
                "source": row["source"],
                "n_pairs": row["n_pairs"],
                "n_phages": row["n_phages"],
                "auc": {"point": row["auc_point"], "lo": row["auc_low"], "hi": row["auc_high"]},
                "brier": {
                    "point": row["brier_point"],
                    "lo": row["brier_low"],
                    "hi": row["brier_high"],
                },
            }
        )
    bacteria_rows = []
    bact_metrics = ch09_report["axis_source_metrics"]["bacteria_axis"]
    for source_label, key in (("guelin", "guelin_raw"), ("basel", "basel_raw")):
        entry = bact_metrics[key]
        bacteria_rows.append(
            {
                "source": source_label,
                "n_pairs": entry["n_pairs"],
                "auc": {"point": entry["auc"], "lo": None, "hi": None},
                "brier": {"point": entry["brier"], "lo": None, "hi": None},
                "ece": entry["ece"],
            }
        )
    return {"bacteria": bacteria_rows, "phage": phage_rows}


def _feature_to_slot(feature: str) -> str:
    """Map a namespaced feature name (e.g. `host_typing__host_serotype`) to its slot key."""
    if feature.startswith(CONCENTRATION_FEATURE_PREFIX):
        return "pair_concentration"
    for spec in SLOT_SPECS:
        if feature.startswith(spec.column_prefix):
            return spec.slot_name
    if feature.startswith("pair_depo_capsule__"):
        return "pair_depo_capsule"
    if feature.startswith("pair_receptor_omp__"):
        return "pair_receptor_omp"
    return "other"


def build_feature_importance_json(ch05_dir: Path) -> dict[str, Any]:
    """Load both per-axis FI CSVs and attach slot tags for coloring in the UI."""
    result: dict[str, Any] = {}
    for axis_name, filename in (
        ("bacteria_axis", "ch05_bacteria_axis_feature_importance.csv"),
        ("phage_axis", "ch05_phage_axis_feature_importance.csv"),
    ):
        path = ch05_dir / filename
        if not path.exists():
            raise FileNotFoundError(
                f"Missing feature importance CSV at {path}. Ensure CH05 ran under the "
                "EX01 patch (ch05_eval.py emits these alongside per-axis metrics)."
            )
        df = pd.read_csv(path)
        df["slot"] = df["feature"].map(_feature_to_slot)
        result[axis_name] = df.to_dict(orient="records")
    return result


def build_predictions_json(ch05_dir: Path, axis_key: str) -> list[dict[str, Any]]:
    """Flatten a per-axis predictions CSV into a compact JSON list for the table tab."""
    csv_name = f"ch05_{axis_key}_predictions.csv"
    path = ch05_dir / csv_name
    if not path.exists():
        raise FileNotFoundError(f"Missing {csv_name} at {path}.")
    df = pd.read_csv(path)
    out = df.loc[
        :,
        [
            "bacteria",
            "phage",
            "source",
            "label_row_binary",
            "predicted_probability",
            "log10_pfu_ml",
        ],
    ].rename(
        columns={
            "label_row_binary": "label",
            "predicted_probability": "predicted",
        }
    )
    return out.to_dict(orient="records")


def build_reliability_json(ch09_dir: Path) -> dict[str, Any]:
    """Group CH09 reliability tables by (axis, source, variant) and attach per-variant ECE."""
    tables_path = ch09_dir / "ch09_reliability_tables.csv"
    report_path = ch09_dir / "ch09_calibration_report.json"
    for path in (tables_path, report_path):
        if not path.exists():
            raise FileNotFoundError(f"Missing CH09 artifact at {path}.")
    tables = pd.read_csv(tables_path)
    with report_path.open(encoding="utf-8") as f:
        report = json.load(f)

    variants: dict[str, dict[str, Any]] = {}
    for variant_name, group in tables.groupby("variant"):
        axis, source, variant_kind = _parse_variant_name(str(variant_name))
        bins = group[["bin_lo", "bin_hi", "n_pairs", "mean_predicted", "observed_rate", "reliability_gap"]].to_dict(
            orient="records"
        )
        variants[str(variant_name)] = {
            "axis": axis,
            "source": source,
            "variant": variant_kind,
            "bins": bins,
        }
    ece_map = _extract_ece_from_report(report)
    for variant_name, body in variants.items():
        body["ece"] = ece_map.get(variant_name)
    return {"variants": variants}


def _parse_variant_name(name: str) -> tuple[str, str, str]:
    """Parse CH09 variant names into (axis, source, variant_kind).

    Examples:
        bacteria_axis_guelin_raw            -> (bacteria, guelin, raw)
        phage_axis_basel_guelin_calibrator_applied -> (phage, basel, guelin_calibrator_applied)
        bacteria_axis_guelin_loof_calibrated -> (bacteria, guelin, loof_calibrated)
    """
    if name.startswith("bacteria_axis_"):
        axis = "bacteria"
        rest = name[len("bacteria_axis_") :]
    elif name.startswith("phage_axis_"):
        axis = "phage"
        rest = name[len("phage_axis_") :]
    else:
        raise ValueError(f"Unrecognized CH09 variant name: {name!r}")
    if rest.startswith("guelin_"):
        source = "guelin"
        variant_kind = rest[len("guelin_") :]
    elif rest.startswith("basel_"):
        source = "basel"
        variant_kind = rest[len("basel_") :]
    else:
        raise ValueError(f"Unrecognized CH09 variant source in: {name!r}")
    return axis, source, variant_kind


def _extract_ece_from_report(report: dict[str, Any]) -> dict[str, float]:
    """Map each reliability variant name to its ECE from the CH09 report."""
    result: dict[str, float] = {}
    for axis_key, axis_block in report["axis_source_metrics"].items():
        for entry_key, entry in axis_block.items():
            result[f"{axis_key}_{entry_key}"] = entry["ece"]
    return result


def build_slot_manifest_json(feature_importance: dict[str, Any]) -> dict[str, Any]:
    """Per-slot metadata: canonical SlotSpec info + per-axis feature count and cumulative importance."""
    slots: dict[str, dict[str, Any]] = {}
    for spec in SLOT_SPECS:
        slots[spec.slot_name] = {
            "slot_name": spec.slot_name,
            "block_role": spec.block_role,
            "entity_key": spec.entity_key,
            "column_prefix": spec.column_prefix,
            "description": spec.description,
        }
    # Pair-level and concentration slots are not in SlotSpec but show up in feature importance.
    for pseudo_slot in ("pair_depo_capsule", "pair_receptor_omp", "pair_concentration", "other"):
        slots.setdefault(
            pseudo_slot,
            {
                "slot_name": pseudo_slot,
                "block_role": "pair" if pseudo_slot.startswith("pair_") else "meta",
                "entity_key": "pair",
                "column_prefix": f"{pseudo_slot}__",
                "description": (
                    "Pairwise cross-term block (host × phage)."
                    if pseudo_slot.startswith("pair_") and pseudo_slot != "pair_concentration"
                    else (
                        "Concentration feature encoded as absolute log₁₀ pfu/ml."
                        if pseudo_slot == "pair_concentration"
                        else "Uncategorized feature source."
                    )
                ),
            },
        )
    # Attach per-axis counts and cumulative importance from FI.
    for axis_key in ("bacteria_axis", "phage_axis"):
        for slot in slots.values():
            slot.setdefault("axes", {})
            slot["axes"][axis_key] = {"feature_count": 0, "cumulative_importance": 0.0}
        for feature in feature_importance[axis_key]:
            slot = slots.get(feature["slot"])
            if slot is None:
                continue
            slot["axes"][axis_key]["feature_count"] += 1
            slot["axes"][axis_key]["cumulative_importance"] += float(feature["mean_importance"])
    return {"slots": list(slots.values())}


def write_snapshot(out_dir: Path, ch05_dir: Path, ch09_dir: Path) -> dict[str, Path]:
    """Write all seven snapshot JSONs to `out_dir`. Returns a name → path map."""
    out_dir.mkdir(parents=True, exist_ok=True)
    combined = _read_combined_summary(ch05_dir)
    report_path = ch09_dir / "ch09_calibration_report.json"
    if not report_path.exists():
        raise FileNotFoundError(f"Missing CH09 report at {report_path}.")
    with report_path.open(encoding="utf-8") as f:
        ch09_report = json.load(f)

    LOGGER.info("Reading CH05 combined summary + CH09 calibration report")
    summary = build_summary_json(combined)
    LOGGER.info("Building cross-source JSON (phage-axis w/ CIs, bacteria-axis point-only)")
    cross_source = build_cross_source_json(combined, ch09_report)
    LOGGER.info("Loading per-axis feature importance")
    feature_importance = build_feature_importance_json(ch05_dir)
    LOGGER.info("Loading bacteria-axis predictions")
    predictions_bact = build_predictions_json(ch05_dir, "bacteria_axis")
    LOGGER.info("Loading phage-axis predictions")
    predictions_phage = build_predictions_json(ch05_dir, "phage_axis")
    LOGGER.info("Building reliability JSON from CH09 tables")
    reliability = build_reliability_json(ch09_dir)
    LOGGER.info("Building slot manifest with per-axis counts + cumulative importance")
    slot_manifest = build_slot_manifest_json(feature_importance)

    snapshot_generated_at = datetime.now(timezone.utc).isoformat()
    summary["snapshot_generated_at"] = snapshot_generated_at

    outputs = {
        "ch05_summary.json": summary,
        "ch05_cross_source.json": cross_source,
        "ch05_feature_importance.json": feature_importance,
        "ch05_predictions_bact.json": predictions_bact,
        "ch05_predictions_phage.json": predictions_phage,
        "ch05_reliability.json": reliability,
        "slot_manifest.json": slot_manifest,
    }
    written: dict[str, Path] = {}
    for name, body in outputs.items():
        path = out_dir / name
        with path.open("w", encoding="utf-8") as f:
            json.dump(body, f, indent=2)
        size_kb = path.stat().st_size / 1024
        LOGGER.info("Wrote %s (%.1f KB)", path, size_kb)
        written[name] = path
    LOGGER.info("Snapshot complete at %s", snapshot_generated_at)
    return written


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--out",
        type=Path,
        default=DEFAULT_OUT_DIR,
        help=f"Output directory for snapshot JSONs (default: {DEFAULT_OUT_DIR})",
    )
    parser.add_argument(
        "--ch05-dir",
        type=Path,
        default=DEFAULT_CH05_DIR,
        help=f"CH05 artifact directory (default: {DEFAULT_CH05_DIR})",
    )
    parser.add_argument(
        "--ch09-dir",
        type=Path,
        default=DEFAULT_CH09_DIR,
        help=f"CH09 artifact directory (default: {DEFAULT_CH09_DIR})",
    )
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    start = datetime.now(timezone.utc)
    LOGGER.info("Snapshot build starting (out=%s, ch05=%s, ch09=%s)", args.out, args.ch05_dir, args.ch09_dir)
    write_snapshot(out_dir=args.out, ch05_dir=args.ch05_dir, ch09_dir=args.ch09_dir)
    elapsed = (datetime.now(timezone.utc) - start).total_seconds()
    LOGGER.info("Snapshot build complete in %.1f s", elapsed)


if __name__ == "__main__":
    main()
