#!/usr/bin/env python3
"""User-facing entry point for short AUTORESEARCH experiment runs."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import time
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
from lightgbm import LGBMClassifier
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.autoresearch.calibration import (
    calibrate_predictions,
    fit_isotonic_calibrator_cv,
)
from lyzortx.autoresearch.per_phage_model import (
    DEFAULT_BLEND_ALPHA,
    PerPhageResult,
    fit_per_phage_models,
    predict_per_phage,
)
from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.runtime_contract import (
    CACHE_CONTRACT_ID,
    CACHE_MANIFEST_FILENAME,
    DEFAULT_CACHE_DIR,
    DEFAULT_OUTPUT_ROOT,
    DISALLOWED_SEARCH_SPLITS,
    INNER_VAL_SPLIT,
    PROVENANCE_MANIFEST_FILENAME,
    SCHEMA_MANIFEST_FILENAME,
    SCHEMA_MANIFEST_ID,
    SEARCH_PAIR_TABLE_ID,
    SLOT_FEATURE_TABLE_FILENAME,
    SLOT_FEATURES_FILENAME,
    SLOT_SCHEMA_FILENAME,
    SLOT_SPEC_BY_NAME,
    SUPPORTED_SEARCH_SPLITS,
    TRAIN_SPLIT,
    sha256_strings,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_TRAIN_OUTPUT_DIR = DEFAULT_OUTPUT_ROOT / "train_runs" / "ar07_baseline"

BASELINE_ID = "ar07_raw_features_lightgbm_pair_scorer_v2"
PRIMARY_SEARCH_METRIC = "roc_auc"
SECONDARY_REPORT_ONLY_METRICS = ("top3_hit_rate", "brier_score")
FIXED_SINGLE_GPU_WALL_CLOCK_BUDGET_SECONDS = 1800
PAIR_SCORER_RANDOM_STATE = 7
PAIR_SCORER_PARAMS = {
    "n_estimators": 300,
    "learning_rate": 0.05,
    "num_leaves": 31,
    "min_child_samples": 10,
    "subsample": 0.8,
    "colsample_bytree": 0.8,
    "reg_alpha": 0.0,
    "reg_lambda": 0.0,
}
REQUIRED_BASELINE_SLOTS = (
    "host_surface",
    "host_typing",
    "host_stats",
    "phage_projection",
    "phage_stats",
)
OPTIONAL_ADDITIVE_ABLATION_SLOTS = ("host_defense",)

SLOT_PREFIXES = tuple(f"{slot}__" for slot in (*REQUIRED_BASELINE_SLOTS, *OPTIONAL_ADDITIVE_ABLATION_SLOTS))


@dataclass(frozen=True)
class SlotArtifact:
    slot_name: str
    entity_key: str
    feature_columns: tuple[str, ...]
    frame: pd.DataFrame


@dataclass(frozen=True)
class CacheContext:
    cache_dir: Path
    cache_manifest: dict[str, Any]
    schema_manifest: dict[str, Any]
    provenance_manifest: dict[str, Any]
    contract_manifest: dict[str, Any]
    split_frames: dict[str, pd.DataFrame]
    slot_artifacts: dict[str, SlotArtifact]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Prepared AUTORESEARCH search cache directory from prepare.py.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_TRAIN_OUTPUT_DIR,
        help="Directory for baseline metrics and prediction outputs.",
    )
    parser.add_argument(
        "--device-type",
        choices=("cpu", "gpu", "cuda"),
        default="cuda",
        help="LightGBM device type for the fixed-budget pair scorer. RunPod search should use one GPU.",
    )
    parser.add_argument(
        "--include-host-defense",
        action="store_true",
        help="Add the reserved host_defense block as an additive ablation on top of the adsorption-first baseline.",
    )
    parser.add_argument(
        "--variant",
        choices=("all-pairs", "per-phage-blend"),
        default="all-pairs",
        help=(
            "Model variant. 'all-pairs' = standard LightGBM on all (phage, host) pairs. "
            "'per-phage-blend' = fit per-phage sub-models on host-only features and blend "
            "with all-pairs predictions (AX02)."
        ),
    )
    parser.add_argument(
        "--blend-alpha",
        type=float,
        default=0.5,
        help="Blending weight for per-phage variant: alpha * per_phage + (1-alpha) * all_pairs. Default 0.5.",
    )
    parser.add_argument(
        "--calibrate",
        choices=("none", "isotonic"),
        default="none",
        help=(
            "Post-hoc calibration method. 'isotonic' = fit isotonic regression on "
            "3-fold cross-validated out-of-fold training predictions (AX05)."
        ),
    )
    return parser.parse_args(argv)


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv_frame(path: Path) -> pd.DataFrame:
    return pd.read_csv(path, dtype=str, keep_default_na=False)


def _json_default(value: object) -> object:
    if isinstance(value, Path):
        return str(value)
    if isinstance(value, np.floating):
        return float(value)
    if isinstance(value, np.integer):
        return int(value)
    return value


def write_json(path: Path, payload: dict[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2, default=_json_default), encoding="utf-8")


def write_rows(path: Path, rows: list[dict[str, object]]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    if not rows:
        raise ValueError(f"Cannot write an empty CSV artifact: {path}")
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def enforce_budget(start_time: float) -> None:
    elapsed = time.monotonic() - start_time
    if elapsed > FIXED_SINGLE_GPU_WALL_CLOCK_BUDGET_SECONDS:
        raise TimeoutError(
            "AUTORESEARCH train.py exceeded the fixed single-GPU wall-clock budget of "
            f"{FIXED_SINGLE_GPU_WALL_CLOCK_BUDGET_SECONDS} seconds."
        )


def validate_cache_manifest(cache_dir: Path) -> tuple[dict[str, Any], dict[str, Any], dict[str, Any]]:
    cache_manifest_path = cache_dir / CACHE_MANIFEST_FILENAME
    schema_manifest_path = cache_dir / SCHEMA_MANIFEST_FILENAME
    provenance_manifest_path = cache_dir / PROVENANCE_MANIFEST_FILENAME
    for path in (cache_manifest_path, schema_manifest_path, provenance_manifest_path):
        if not path.exists():
            raise FileNotFoundError(
                f"Prepared search cache is incomplete at {cache_dir}. Missing required manifest: {path.name}"
            )

    cache_manifest = read_json(cache_manifest_path)
    schema_manifest = read_json(schema_manifest_path)
    provenance_manifest = read_json(provenance_manifest_path)

    if cache_manifest.get("cache_contract_id") != CACHE_CONTRACT_ID:
        raise ValueError(
            "AUTORESEARCH cache contract mismatch: "
            f"expected {CACHE_CONTRACT_ID}, got {cache_manifest.get('cache_contract_id')!r}"
        )
    if schema_manifest.get("schema_manifest_id") != SCHEMA_MANIFEST_ID:
        raise ValueError(
            "AUTORESEARCH schema manifest mismatch: "
            f"expected {SCHEMA_MANIFEST_ID}, got {schema_manifest.get('schema_manifest_id')!r}"
        )
    if schema_manifest.get("pair_table_id") != SEARCH_PAIR_TABLE_ID:
        raise ValueError(
            "AUTORESEARCH pair-table contract mismatch: "
            f"expected {SEARCH_PAIR_TABLE_ID}, got {schema_manifest.get('pair_table_id')!r}"
        )
    return cache_manifest, schema_manifest, provenance_manifest


def resolve_slot_features_path(*, slot_dir: Path, slot_name: str) -> Path:
    for candidate in (SLOT_FEATURES_FILENAME, SLOT_FEATURE_TABLE_FILENAME):
        path = slot_dir / candidate
        if path.exists():
            return path
    return slot_dir / SLOT_FEATURES_FILENAME


def load_slot_artifact(
    *,
    cache_dir: Path,
    schema_manifest: dict[str, Any],
    cache_manifest: dict[str, Any],
    slot_name: str,
    require_materialized_features: bool,
) -> SlotArtifact:
    slot_spec = SLOT_SPEC_BY_NAME.get(slot_name)
    if slot_spec is None:
        raise ValueError(f"Unrecognized AUTORESEARCH slot: {slot_name}")
    top_level_slot = schema_manifest.get("feature_slots", {}).get(slot_name)
    if top_level_slot is None:
        raise ValueError(f"AUTORESEARCH schema manifest does not define slot {slot_name}.")
    slot_summary = (cache_manifest.get("feature_slots") or {}).get(slot_name)
    if slot_summary is None:
        raise ValueError(f"AUTORESEARCH cache manifest is missing slot {slot_name}.")

    schema_path = Path(str(slot_summary["schema_manifest_path"]))
    if not schema_path.exists():
        fallback_schema_path = cache_dir / "feature_slots" / slot_name / SLOT_SCHEMA_FILENAME
        if fallback_schema_path.exists():
            schema_path = fallback_schema_path
    if not schema_path.exists():
        raise FileNotFoundError(f"AUTORESEARCH slot schema manifest not found for {slot_name}: {schema_path}")
    slot_schema = read_json(schema_path)
    if slot_schema.get("schema_manifest_id") != SCHEMA_MANIFEST_ID:
        raise ValueError(f"AUTORESEARCH slot {slot_name} no longer points at the frozen schema manifest.")
    if slot_schema.get("slot_name") != slot_name:
        raise ValueError(
            f"AUTORESEARCH slot schema mismatch: expected {slot_name}, got {slot_schema.get('slot_name')}."
        )
    if slot_schema.get("join_keys") != slot_spec.join_keys:
        raise ValueError(f"AUTORESEARCH slot {slot_name} schema no longer matches the frozen join key.")

    feature_columns = tuple(str(column) for column in slot_schema["reserved_feature_columns"])
    if tuple(str(column) for column in top_level_slot["reserved_feature_columns"]) != feature_columns:
        raise ValueError(f"AUTORESEARCH slot {slot_name} bypassed the frozen top-level cache schema.")

    slot_dir = cache_dir / "feature_slots" / slot_name
    features_path = resolve_slot_features_path(slot_dir=slot_dir, slot_name=slot_name)
    if require_materialized_features and not feature_columns:
        raise ValueError(f"AUTORESEARCH baseline requires materialized columns for slot {slot_name}.")
    if not feature_columns:
        return SlotArtifact(
            slot_name=slot_name, entity_key=slot_spec.entity_key, feature_columns=(), frame=pd.DataFrame()
        )
    if not features_path.exists():
        raise FileNotFoundError(
            f"AUTORESEARCH slot {slot_name} is missing its materialized feature table at {features_path}."
        )

    frame = load_csv_frame(features_path)
    expected_header = [slot_spec.entity_key, *feature_columns]
    actual_header = list(frame.columns)
    if actual_header != expected_header:
        raise ValueError(
            f"AUTORESEARCH slot {slot_name} features.csv header mismatch. "
            f"Expected {len(expected_header)} columns, got {len(actual_header)}."
        )
    return SlotArtifact(
        slot_name=slot_name,
        entity_key=slot_spec.entity_key,
        feature_columns=feature_columns,
        frame=frame,
    )


def load_and_validate_cache(
    *,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    include_host_defense: bool = False,
) -> CacheContext:
    cache_manifest, schema_manifest, provenance_manifest = validate_cache_manifest(cache_dir)

    contract_manifest_path = Path(
        str(cache_manifest.get("provenance_manifest_path", "")).replace(
            PROVENANCE_MANIFEST_FILENAME, "ar01_split_benchmark_manifest_v1.json"
        )
    )
    if not contract_manifest_path.exists():
        candidate = cache_dir.parent / "ar01_split_benchmark_manifest_v1.json"
        if candidate.exists():
            contract_manifest_path = candidate
    contract_manifest = read_json(contract_manifest_path) if contract_manifest_path.exists() else {}

    split_frames: dict[str, pd.DataFrame] = {}
    pair_tables = cache_manifest.get("pair_tables", {})
    for split_name in SUPPORTED_SEARCH_SPLITS:
        split_info = pair_tables.get(split_name)
        if split_info is None:
            raise ValueError(f"Prepared cache is missing the required pair table for split '{split_name}'.")
        split_path = Path(str(split_info["path"]))
        if not split_path.exists():
            raise FileNotFoundError(f"Pair table for split '{split_name}' not found: {split_path}")
        split_frames[split_name] = load_csv_frame(split_path)

    selected_slots = list(REQUIRED_BASELINE_SLOTS)
    if include_host_defense:
        selected_slots.append("host_defense")
    slot_artifacts = {
        slot_name: load_slot_artifact(
            cache_dir=cache_dir,
            schema_manifest=schema_manifest,
            cache_manifest=cache_manifest,
            slot_name=slot_name,
            require_materialized_features=True,
        )
        for slot_name in selected_slots
    }

    return CacheContext(
        cache_dir=cache_dir,
        cache_manifest=cache_manifest,
        schema_manifest=schema_manifest,
        provenance_manifest=provenance_manifest,
        contract_manifest=contract_manifest,
        split_frames=split_frames,
        slot_artifacts=slot_artifacts,
    )


def detect_feature_types(
    frame: pd.DataFrame, feature_columns: Sequence[str]
) -> tuple[list[str], list[str], pd.DataFrame]:
    typed = frame.copy()
    numeric_columns: list[str] = []
    categorical_columns: list[str] = []
    for column in feature_columns:
        numeric_candidate = pd.to_numeric(typed[column], errors="coerce")
        non_empty_mask = typed[column].astype(str) != ""
        if bool(non_empty_mask.any()) and bool(numeric_candidate[non_empty_mask].notna().all()):
            typed[column] = numeric_candidate.astype(float)
            numeric_columns.append(column)
        else:
            typed[column] = typed[column].astype("category")
            categorical_columns.append(column)
    if not numeric_columns and not categorical_columns:
        raise ValueError("AUTORESEARCH baseline requires at least one feature column.")
    return numeric_columns, categorical_columns, typed


def build_entity_feature_table(
    slot_artifacts: dict[str, SlotArtifact],
    *,
    slot_names: Sequence[str],
    entity_key: str,
) -> pd.DataFrame:
    merged: pd.DataFrame | None = None
    for slot_name in slot_names:
        artifact = slot_artifacts[slot_name]
        if artifact.entity_key != entity_key:
            raise ValueError(f"AUTORESEARCH slot {slot_name} does not join on {entity_key}.")
        frame = artifact.frame.copy()
        if merged is None:
            merged = frame
        else:
            merged = merged.merge(frame, on=entity_key, how="outer", validate="one_to_one")
    if merged is None or merged.empty:
        raise ValueError(f"AUTORESEARCH entity table for {entity_key} is empty.")
    return merged


def type_entity_features(entity_table: pd.DataFrame, entity_key: str) -> tuple[pd.DataFrame, list[str], list[str]]:
    """Detect feature types and convert columns in place. Returns (typed_table, numeric_cols, categorical_cols)."""
    feature_columns = [col for col in entity_table.columns if col != entity_key]
    numeric_columns, categorical_columns, typed = detect_feature_types(entity_table, feature_columns)
    return typed, numeric_columns, categorical_columns


def build_raw_pair_design_matrix(
    pair_frame: pd.DataFrame,
    *,
    host_features: pd.DataFrame,
    phage_features: pd.DataFrame,
) -> pd.DataFrame:
    """Merge raw host and phage features onto pair rows. No SVD, no interaction terms."""
    merged = pair_frame.merge(host_features, on="bacteria", how="left", validate="many_to_one")
    merged = merged.merge(phage_features, on="phage", how="left", validate="many_to_one")

    feature_cols = [col for col in merged.columns if col.startswith(SLOT_PREFIXES)]
    if not feature_cols:
        raise ValueError("AUTORESEARCH pair design matrix has zero feature columns after merge.")
    # Check entity key join completeness (every pair must join to both a host and phage row).
    # NaN in individual feature columns is acceptable — LightGBM handles missing values natively.
    entity_key_cols = [col for col in ["bacteria", "phage"] if col in merged.columns]
    for key_col in entity_key_cols:
        prefix = "host_" if key_col == "bacteria" else "phage_"
        key_features = [col for col in feature_cols if col.startswith(prefix)]
        if key_features and merged[key_features].isna().all(axis=1).any():
            raise ValueError(f"AUTORESEARCH pair table has rows with no {key_col} features — join failure.")
    return merged


def build_pair_scorer(device_type: str) -> LGBMClassifier:
    estimator_params: dict[str, Any] = {
        **PAIR_SCORER_PARAMS,
        "objective": "binary",
        "class_weight": "balanced",
        "random_state": PAIR_SCORER_RANDOM_STATE,
        "n_jobs": 1,
        "verbosity": -1,
        "device_type": device_type,
    }
    if device_type == "cpu":
        estimator_params["deterministic"] = True
        estimator_params["force_col_wise"] = True
    return LGBMClassifier(**estimator_params)


def compute_top3_hit_rate(scored_rows: pd.DataFrame) -> float:
    if scored_rows.empty:
        raise ValueError("Cannot compute top-3 hit rate on an empty score table.")

    hit_flags: list[float] = []
    for _, bacteria_rows in scored_rows.sort_values(
        ["bacteria", "prediction", "phage"], ascending=[True, False, True]
    ).groupby(
        "bacteria",
        sort=True,
    ):
        top_rows = bacteria_rows.head(3)
        hit_flags.append(float((top_rows["label_any_lysis"].astype(int) == 1).any()))
    return float(sum(hit_flags) / len(hit_flags))


def safe_float(value: float) -> float:
    return float(f"{value:.6f}")


def validate_no_holdout_leakage(cache_dir: Path) -> None:
    """Ensure no holdout labels leaked into the search cache."""
    search_pairs_dir = cache_dir / "search_pairs"
    if not search_pairs_dir.exists():
        return
    for csv_path in search_pairs_dir.glob("*.csv"):
        frame = load_csv_frame(csv_path)
        if "split" in frame.columns:
            leaked = set(frame["split"].unique()) & set(DISALLOWED_SEARCH_SPLITS)
            if leaked:
                raise ValueError(
                    f"AUTORESEARCH search cache contains holdout labels at {csv_path.name} "
                    f"(disallowed splits: {leaked})"
                )


def validate_split_membership(context: CacheContext) -> None:
    """Verify pair IDs match the AR01 contract to detect split membership drift."""
    contract = context.contract_manifest
    if not contract or "split_contract" not in contract:
        return
    split_hashes = contract["split_contract"].get("split_hashes", {})
    for split_name, split_frame in context.split_frames.items():
        retained = split_frame.loc[split_frame["retained_for_autoresearch"] == "1"]
        expected_hash = split_hashes.get(split_name, {}).get("retained_pair_ids_sha256")
        if expected_hash is not None:
            actual_hash = sha256_strings(sorted(retained["pair_id"].tolist()))
            if actual_hash != expected_hash:
                raise ValueError(
                    f"AUTORESEARCH split membership drift detected for '{split_name}': "
                    f"expected {expected_hash[:12]}..., got {actual_hash[:12]}..."
                )


def run_baseline(
    *,
    context: CacheContext,
    device_type: str,
    include_host_defense: bool,
    variant: str = "all-pairs",
    blend_alpha: float = DEFAULT_BLEND_ALPHA,
    calibrate: str = "none",
    start_time: float,
) -> tuple[dict[str, Any], list[dict[str, object]]]:
    enforce_budget(start_time)

    host_slots = ["host_surface", "host_typing", "host_stats"]
    if include_host_defense:
        host_slots.append("host_defense")
    phage_slots = ["phage_projection", "phage_stats"]

    host_table = build_entity_feature_table(context.slot_artifacts, slot_names=host_slots, entity_key="bacteria")
    phage_table = build_entity_feature_table(context.slot_artifacts, slot_names=phage_slots, entity_key="phage")

    host_typed, host_numeric, host_categorical = type_entity_features(host_table, "bacteria")
    phage_typed, phage_numeric, phage_categorical = type_entity_features(phage_table, "phage")

    train_pairs = context.split_frames[TRAIN_SPLIT].loc[lambda frame: frame["retained_for_autoresearch"] == "1"].copy()
    inner_val_pairs = (
        context.split_frames[INNER_VAL_SPLIT].loc[lambda frame: frame["retained_for_autoresearch"] == "1"].copy()
    )

    train_design = build_raw_pair_design_matrix(train_pairs, host_features=host_typed, phage_features=phage_typed)
    inner_val_design = build_raw_pair_design_matrix(
        inner_val_pairs, host_features=host_typed, phage_features=phage_typed
    )

    feature_columns = [col for col in train_design.columns if col.startswith(SLOT_PREFIXES)]
    if not feature_columns:
        raise ValueError("AUTORESEARCH baseline constructed zero pair features.")

    numeric_feature_columns = host_numeric + phage_numeric
    categorical_feature_columns = host_categorical + phage_categorical
    categorical_in_features = [col for col in categorical_feature_columns if col in feature_columns]

    y_train = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    y_inner = inner_val_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)

    enforce_budget(start_time)
    estimator = build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[feature_columns],
        y_train,
        sample_weight=sample_weight,
        categorical_feature=categorical_in_features,
    )
    enforce_budget(start_time)

    all_pairs_predictions = estimator.predict_proba(inner_val_design[feature_columns])[:, 1]

    # Per-phage blending (AX02): fit per-phage sub-models on host-only features and blend.
    per_phage_result: PerPhageResult | None = None
    if variant == "per-phage-blend":
        host_feature_columns = [
            col
            for col in feature_columns
            if col.startswith(("host_surface__", "host_typing__", "host_stats__", "host_defense__"))
        ]
        if not host_feature_columns:
            raise ValueError("Per-phage variant requires host feature columns but none found in design matrix.")
        LOGGER.info(
            "Fitting per-phage sub-models on %d host features for %d training phages",
            len(host_feature_columns),
            train_design["phage"].nunique(),
        )
        per_phage_models = fit_per_phage_models(
            train_design,
            host_feature_columns,
            device_type=device_type,
        )
        predictions, per_phage_result = predict_per_phage(
            per_phage_models,
            inner_val_design,
            host_feature_columns,
            all_pairs_predictions,
            blend_alpha=blend_alpha,
        )
        LOGGER.info(
            "Per-phage blending: %d/%d phages fitted, %d fallback (alpha=%.2f)",
            per_phage_result.n_phages_fitted,
            per_phage_result.n_phages_total,
            per_phage_result.n_phages_fallback,
            blend_alpha,
        )
        enforce_budget(start_time)
    else:
        predictions = all_pairs_predictions

    # Post-hoc calibration (AX05): fit on OOF training predictions, apply to inner_val.
    calibration_applied = False
    if calibrate == "isotonic":
        LOGGER.info("Fitting isotonic calibrator via 3-fold CV on training predictions")
        calibrator = fit_isotonic_calibrator_cv(
            train_design,
            feature_columns,
            y_train,
            sample_weight,
            categorical_features=categorical_in_features,
            device_type=device_type,
            model_params=PAIR_SCORER_PARAMS,
        )
        predictions = calibrate_predictions(calibrator, predictions)
        calibration_applied = True
        LOGGER.info("Applied isotonic calibration to inner_val predictions")
        enforce_budget(start_time)

    scored_inner_val = inner_val_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].copy()
    scored_inner_val["prediction"] = predictions
    scored_inner_val["rank_within_bacteria"] = (
        scored_inner_val.sort_values(["bacteria", "prediction", "phage"], ascending=[True, False, True])
        .groupby("bacteria")
        .cumcount()
        .add(1)
        .sort_index()
    )

    if len(np.unique(y_inner)) < 2:
        raise ValueError("AUTORESEARCH inner_val split must contain both classes to produce ROC-AUC.")

    metrics = {
        "roc_auc": safe_float(float(roc_auc_score(y_inner, predictions))),
        "top3_hit_rate": safe_float(compute_top3_hit_rate(scored_inner_val)),
        "brier_score": safe_float(float(brier_score_loss(y_inner, predictions))),
    }
    summary: dict[str, Any] = {
        "baseline_id": BASELINE_ID,
        "variant": variant,
        "feature_space": {
            "type": "raw_slot_features",
            "host_slots": host_slots,
            "phage_slots": phage_slots,
            "numeric_columns": numeric_feature_columns,
            "categorical_columns": categorical_feature_columns,
            "total_feature_count": len(feature_columns),
        },
        "calibration": {
            "method": calibrate,
            "applied": calibration_applied,
        },
        "pair_scorer": {
            "type": "lightgbm_binary_classifier",
            "device_type": device_type,
            "params": dict(PAIR_SCORER_PARAMS),
        },
        "feature_columns": feature_columns,
        "metrics": metrics,
    }
    if per_phage_result is not None:
        summary["per_phage"] = {
            "blend_alpha": blend_alpha,
            "n_phages_fitted": per_phage_result.n_phages_fitted,
            "n_phages_fallback": per_phage_result.n_phages_fallback,
            "fitted_phages": list(per_phage_result.fitted_phages),
            "fallback_phages": list(per_phage_result.fallback_phages),
        }
    prediction_rows = [
        {
            "pair_id": str(row["pair_id"]),
            "bacteria": str(row["bacteria"]),
            "phage": str(row["phage"]),
            "label_any_lysis": int(row["label_any_lysis"]),
            "prediction": safe_float(float(row["prediction"])),
            "rank_within_bacteria": int(row["rank_within_bacteria"]),
        }
        for row in scored_inner_val.sort_values(["bacteria", "rank_within_bacteria", "phage"]).to_dict(orient="records")
    ]
    return summary, prediction_rows


def main(argv: Optional[Sequence[str]] = None) -> int:
    setup_logging()
    args = parse_args(argv)
    start_time = time.monotonic()

    context = load_and_validate_cache(
        cache_dir=args.cache_dir,
        include_host_defense=args.include_host_defense,
    )
    validate_no_holdout_leakage(context.cache_dir)
    validate_split_membership(context)
    LOGGER.info("AUTORESEARCH train sandbox validated cache at %s", args.cache_dir)
    LOGGER.info("Exported splits: %s", ", ".join(SUPPORTED_SEARCH_SPLITS))
    active_optional = []
    if args.include_host_defense:
        active_optional.append("host_defense")
    LOGGER.info(
        "Baseline feature slots: %s",
        ", ".join(REQUIRED_BASELINE_SLOTS + tuple(active_optional)),
    )
    LOGGER.info("Reserved-but-ignored baseline ablations: %s", ", ".join(OPTIONAL_ADDITIVE_ABLATION_SLOTS))
    LOGGER.info("Fixed single-GPU wall-clock budget: %d seconds", FIXED_SINGLE_GPU_WALL_CLOCK_BUDGET_SECONDS)
    LOGGER.info("train.py is the short-loop experiment surface; cache rebuilding belongs in prepare.py.")

    LOGGER.info("Model variant: %s", args.variant)
    if args.calibrate != "none":
        LOGGER.info("Calibration: %s", args.calibrate)
    baseline_summary, prediction_rows = run_baseline(
        context=context,
        device_type=args.device_type,
        include_host_defense=args.include_host_defense,
        variant=args.variant,
        blend_alpha=args.blend_alpha,
        calibrate=args.calibrate,
        start_time=start_time,
    )
    elapsed_seconds = safe_float(time.monotonic() - start_time)

    args.output_dir.mkdir(parents=True, exist_ok=True)
    prediction_path = args.output_dir / "ar07_inner_val_predictions.csv"
    summary_path = args.output_dir / "ar07_baseline_summary.json"
    write_rows(prediction_path, prediction_rows)
    output_summary = {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_dir": str(args.cache_dir),
        "output_dir": str(args.output_dir),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "schema_manifest_id": SCHEMA_MANIFEST_ID,
        "search_mutation_boundary": {
            "mutable_file": "lyzortx/autoresearch/train.py",
            "labels_out_of_bounds": True,
            "split_membership_out_of_bounds": True,
            "feature_extraction_out_of_bounds": True,
            "evaluation_contract_out_of_bounds": True,
        },
        "search_runtime_contract": {
            "device_type": args.device_type,
            "fixed_wall_clock_budget_seconds": FIXED_SINGLE_GPU_WALL_CLOCK_BUDGET_SECONDS,
            "cache_build_outside_budget": True,
            "cache_rebuild_allowed_in_train": False,
        },
        "baseline_contract": {
            "minimum_slots": list(REQUIRED_BASELINE_SLOTS),
            "host_defense_active": args.include_host_defense,
            "host_defense_role": "reserved_additive_ablation",
            "variant": args.variant,
            "blend_alpha": args.blend_alpha,
            "calibrate": args.calibrate,
            "primary_metric": PRIMARY_SEARCH_METRIC,
            "secondary_report_only_metrics": list(SECONDARY_REPORT_ONLY_METRICS),
        },
        "artifacts": {
            "inner_val_predictions_path": str(prediction_path),
        },
        "inner_val_metrics": dict(baseline_summary["metrics"]),
        "search_metric": {
            "name": PRIMARY_SEARCH_METRIC,
            "value": baseline_summary["metrics"][PRIMARY_SEARCH_METRIC],
        },
        "elapsed_seconds": elapsed_seconds,
        **baseline_summary,
    }
    write_json(summary_path, output_summary)

    LOGGER.info(
        "Primary search metric: inner_val_%s=%.6f", PRIMARY_SEARCH_METRIC, output_summary["search_metric"]["value"]
    )
    LOGGER.info(
        "Secondary report-only metrics: top3_hit_rate=%.6f, brier_score=%.6f",
        *[output_summary["inner_val_metrics"][metric_name] for metric_name in SECONDARY_REPORT_ONLY_METRICS],
    )
    LOGGER.info("Wrote AUTORESEARCH baseline summary to %s", summary_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
