#!/usr/bin/env python3
"""AR09: import RunPod candidates and replay them on the sealed AUTORESEARCH holdout."""

from __future__ import annotations

import argparse
import csv
import hashlib
import importlib.util
import json
import logging
import shutil
import tarfile
import tempfile
from collections import defaultdict
from contextlib import contextmanager
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Iterator, Mapping, Optional, Sequence

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch import build_contract, prepare_cache, runtime_contract
from lyzortx.pipeline.steel_thread_v0.steps._io_helpers import read_csv_rows, safe_round
from lyzortx.pipeline.track_c.steps.build_v1_host_feature_pair_table import DEFENSE_DERIVED_COLUMNS
from lyzortx.pipeline.track_g.steps import train_v1_binary_classifier

LOGGER = logging.getLogger(__name__)

DEFAULT_ST02_PAIR_TABLE_PATH = Path("lyzortx/generated_outputs/steel_thread_v0/intermediate/st02_pair_table.csv")
DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH = Path(
    "lyzortx/generated_outputs/steel_thread_v0/intermediate/st03_split_assignments.csv"
)
DEFAULT_CANDIDATES_DIR = runtime_contract.DEFAULT_OUTPUT_ROOT / "candidates"
DEFAULT_DECISION_BUNDLES_DIR = runtime_contract.DEFAULT_OUTPUT_ROOT / "decision_bundles"
IMPORT_MANIFEST_FILENAME = "ar09_import_manifest.json"
DECISION_BUNDLE_FILENAME = "ar09_decision_bundle.json"
SEED_METRICS_FILENAME = "ar09_seed_metrics.csv"
AGGREGATED_PREDICTIONS_FILENAME = "ar09_aggregated_holdout_predictions.csv"

DEFAULT_REPLICATION_SEEDS = (7, 17, 29)
DEFAULT_BOOTSTRAP_SAMPLES = 1000
PRIMARY_METRIC = "holdout_roc_auc"
BASELINE_ARM_ID = "locked_production_intent_comparator"
CANDIDATE_ARM_ID = "imported_candidate"

REQUIRED_CANDIDATE_FILES = (
    "train.py",
    "local_run_metadata.json",
    "autoresearch_runpod_bundle_manifest.json",
    "runpod_experiment.log",
    "runpod_execution_metadata.json",
)
OPTIONAL_CANDIDATE_FILES = (
    "ar07_baseline_summary.json",
    "ar07_inner_val_predictions.csv",
)
DEFAULT_FEATURE_LOCK_PATH = Path(
    "lyzortx/generated_outputs/track_g/tg05_feature_subset_sweep/tg05_locked_v1_feature_config.json"
)
DEFAULT_TG01_MODEL_SUMMARY_PATH = Path(build_contract.COMPARATOR_BENCHMARK["model_summary_path"])
FALLBACK_FEATURE_LOCK_PATH = Path("lyzortx/pipeline/track_g/v1_feature_configuration.json")
FALLBACK_TG01_BEST_PARAMS = {
    "learning_rate": 0.05,
    "min_child_samples": 25,
    "n_estimators": 300,
    "num_leaves": 31,
}


BOOTSTRAP_METRIC_NAMES = ("holdout_roc_auc", "holdout_brier_score")


@dataclass(frozen=True)
class BootstrapMetricCI:
    point_estimate: Optional[float]
    ci_low: Optional[float]
    ci_high: Optional[float]
    bootstrap_samples_requested: int
    bootstrap_samples_used: int


def load_v1_lock(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    locked_blocks = list(payload.get("winner_subset_blocks", []))
    if locked_blocks != ["defense", "phage_genomic"]:
        raise ValueError(
            "TL05 expects the current locked v1 baseline to be defense + phage_genomic; "
            f"found {locked_blocks!r} in {path}"
        )
    return dict(payload)


def load_tg01_lock(path: Path) -> dict[str, object]:
    with path.open("r", encoding="utf-8") as handle:
        payload = json.load(handle)
    return {
        "best_params": dict(payload["lightgbm"]["best_params"]),
        "holdout_binary_metrics": dict(payload["lightgbm"]["holdout_binary_metrics"]),
        "holdout_top3_metrics": dict(payload["lightgbm"]["holdout_top3_metrics"]),
    }


def partition_track_c_defense_columns(track_c_columns: Sequence[str]) -> tuple[str, ...]:
    partitioned = train_v1_binary_classifier.deduplicate_preserving_order(track_c_columns)
    if not partitioned:
        raise ValueError("Track C feature table has no columns.")
    defense_columns = [
        col for col in partitioned if col.startswith("host_defense_subtype_") or col in DEFENSE_DERIVED_COLUMNS
    ]
    if not defense_columns:
        raise ValueError("No defense subtype columns were found in Track C features.")
    return tuple(defense_columns)


def _evaluate_holdout_rows(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    y_true = [int(str(row["label_hard_any_lysis"])) for row in rows]
    y_prob = [float(row["predicted_probability"]) for row in rows]
    return {
        "binary": train_v1_binary_classifier.compute_binary_metrics(y_true, y_prob),
    }


def bootstrap_holdout_metric_cis(
    holdout_rows_by_arm: Mapping[str, Sequence[Mapping[str, object]]],
    *,
    bootstrap_samples: int,
    bootstrap_random_state: int,
    baseline_arm_id: str,
) -> dict[str, dict[str, BootstrapMetricCI]]:
    if bootstrap_samples < 1:
        raise ValueError("bootstrap_samples must be >= 1")

    if baseline_arm_id not in holdout_rows_by_arm:
        raise ValueError("Missing baseline arm for bootstrap evaluation.")

    holdout_by_bacteria: dict[str, list[Mapping[str, object]]] = defaultdict(list)
    for row in holdout_rows_by_arm[baseline_arm_id]:
        holdout_by_bacteria[str(row["bacteria"])].append(row)

    bacteria_ids = tuple(sorted(holdout_by_bacteria.keys()))
    if not bacteria_ids:
        raise ValueError("No holdout bacteria available for bootstrap evaluation.")

    arm_bacteria_sets: dict[str, dict[str, list[Mapping[str, object]]]] = {
        arm_id: {bacteria: [] for bacteria in bacteria_ids} for arm_id in holdout_rows_by_arm
    }
    for arm_id, rows in holdout_rows_by_arm.items():
        for row in rows:
            bacteria = str(row["bacteria"])
            if bacteria in arm_bacteria_sets[arm_id]:
                arm_bacteria_sets[arm_id][bacteria].append(row)
    for arm_id, bacteria_map in arm_bacteria_sets.items():
        missing = [bacteria for bacteria in bacteria_ids if not bacteria_map.get(bacteria)]
        if missing:
            raise ValueError(f"Missing holdout rows for arm {arm_id}: {', '.join(missing)}")

    rng = np.random.default_rng(bootstrap_random_state)
    metric_samples: dict[str, dict[str, list[float]]] = {
        arm_id: {metric_name: [] for metric_name in BOOTSTRAP_METRIC_NAMES} for arm_id in holdout_rows_by_arm
    }
    delta_samples: dict[str, dict[str, list[float]]] = {
        f"{arm_id}__delta_vs_{baseline_arm_id}": {metric_name: [] for metric_name in BOOTSTRAP_METRIC_NAMES}
        for arm_id in holdout_rows_by_arm
        if arm_id != baseline_arm_id
    }

    bacteria_count = len(bacteria_ids)
    progress_interval = max(1, bootstrap_samples // 5)
    for sample_index in range(bootstrap_samples):
        if sample_index == 0 or (sample_index + 1) % progress_interval == 0 or sample_index + 1 == bootstrap_samples:
            LOGGER.info(
                "Bootstrap progress: %d/%d paired holdout-strain resamples",
                sample_index + 1,
                bootstrap_samples,
            )
        sampled_bacteria_indices = rng.integers(0, bacteria_count, size=bacteria_count)
        sampled_rows_by_arm: dict[str, list[Mapping[str, object]]] = {}
        for arm_id, bacteria_map in arm_bacteria_sets.items():
            sampled_rows: list[Mapping[str, object]] = []
            for bacteria_index in sampled_bacteria_indices.tolist():
                sampled_rows.extend(bacteria_map[bacteria_ids[bacteria_index]])
            sampled_rows_by_arm[arm_id] = sampled_rows

        metrics_by_arm = {arm_id: _evaluate_holdout_rows(rows) for arm_id, rows in sampled_rows_by_arm.items()}
        baseline_metrics = metrics_by_arm[baseline_arm_id]
        for arm_id, metrics in metrics_by_arm.items():
            metric_samples[arm_id]["holdout_brier_score"].append(float(metrics["binary"]["brier_score"]))
            if metrics["binary"]["roc_auc"] is not None:
                metric_samples[arm_id]["holdout_roc_auc"].append(float(metrics["binary"]["roc_auc"]))

        for arm_id, metrics in metrics_by_arm.items():
            if arm_id == baseline_arm_id:
                continue
            delta_key = f"{arm_id}__delta_vs_{baseline_arm_id}"
            if baseline_metrics["binary"]["roc_auc"] is not None and metrics["binary"]["roc_auc"] is not None:
                delta_samples[delta_key]["holdout_roc_auc"].append(
                    float(metrics["binary"]["roc_auc"]) - float(baseline_metrics["binary"]["roc_auc"])
                )
            delta_samples[delta_key]["holdout_brier_score"].append(
                float(baseline_metrics["binary"]["brier_score"]) - float(metrics["binary"]["brier_score"])
            )

    def _ci(values: Sequence[float]) -> tuple[Optional[float], Optional[float], int]:
        if not values:
            return None, None, 0
        low, high = np.quantile(np.asarray(values, dtype=float), [0.025, 0.975])
        return safe_round(float(low)), safe_round(float(high)), len(values)

    actual_metrics_by_arm = {arm_id: _evaluate_holdout_rows(rows) for arm_id, rows in holdout_rows_by_arm.items()}

    ci_summary: dict[str, dict[str, BootstrapMetricCI]] = {}
    for arm_id, samples in metric_samples.items():
        actual_metrics = actual_metrics_by_arm[arm_id]
        ci_summary[arm_id] = {}
        for metric_name, sample_values in samples.items():
            ci_low, ci_high, used = _ci(sample_values)
            ci_summary[arm_id][metric_name] = BootstrapMetricCI(
                point_estimate=(
                    float(actual_metrics["binary"]["roc_auc"])
                    if metric_name == "holdout_roc_auc"
                    else float(actual_metrics["binary"]["brier_score"])
                ),
                ci_low=ci_low,
                ci_high=ci_high,
                bootstrap_samples_requested=bootstrap_samples,
                bootstrap_samples_used=used,
            )

    for delta_key, samples in delta_samples.items():
        ci_summary[delta_key] = {}
        for metric_name, sample_values in samples.items():
            ci_low, ci_high, used = _ci(sample_values)
            ci_summary[delta_key][metric_name] = BootstrapMetricCI(
                point_estimate=None,
                ci_low=ci_low,
                ci_high=ci_high,
                bootstrap_samples_requested=bootstrap_samples,
                bootstrap_samples_used=used,
            )

    return ci_summary


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    subparsers = parser.add_subparsers(dest="command", required=True)

    import_parser = subparsers.add_parser(
        "import-runpod-candidate",
        help="Import an AR08 RunPod candidate bundle into lyzortx/generated_outputs/autoresearch/candidates/.",
    )
    import_parser.add_argument(
        "--candidate-bundle", type=Path, required=True, help="Candidate bundle directory or tgz."
    )
    import_parser.add_argument("--candidate-id", help="Stable identifier for the imported candidate.")
    import_parser.add_argument(
        "--candidates-dir",
        type=Path,
        default=DEFAULT_CANDIDATES_DIR,
        help="Destination root for imported candidate bundles.",
    )

    replicate_parser = subparsers.add_parser(
        "replicate",
        help="Replay an imported candidate and the locked comparator on the sealed holdout.",
    )
    replicate_parser.add_argument("--candidate-dir", type=Path, required=True, help="Imported candidate directory.")
    replicate_parser.add_argument(
        "--cache-dir",
        type=Path,
        default=runtime_contract.DEFAULT_CACHE_DIR,
        help="Prepared AUTORESEARCH cache directory. The default path will be regenerated if missing.",
    )
    replicate_parser.add_argument(
        "--output-dir",
        type=Path,
        help="Destination directory for the auditable AR09 decision bundle.",
    )
    replicate_parser.add_argument(
        "--device-type",
        choices=("cpu", "gpu", "cuda"),
        default="cpu",
        help="LightGBM device type for repeated-seed replay.",
    )
    replicate_parser.add_argument(
        "--replication-seeds",
        type=int,
        nargs="+",
        default=list(DEFAULT_REPLICATION_SEEDS),
        help="Repeated random seeds used for both candidate and comparator replay.",
    )
    replicate_parser.add_argument(
        "--bootstrap-samples",
        type=int,
        default=DEFAULT_BOOTSTRAP_SAMPLES,
        help="Paired holdout-strain bootstrap samples for arm deltas.",
    )
    replicate_parser.add_argument(
        "--bootstrap-random-state",
        type=int,
        default=42,
        help="Random seed for paired holdout-strain bootstrap confidence intervals.",
    )
    replicate_parser.add_argument(
        "--use-st03-split",
        action="store_true",
        help="Evaluate on ST03 holdout (same 65-bacteria split as TL18) instead of AR01 holdout. Skips comparator.",
    )
    replicate_parser.add_argument(
        "--include-host-defense",
        action="store_true",
        help="Replay the optional host_defense block on top of the imported train.py baseline.",
    )
    replicate_parser.add_argument(
        "--include-pairwise-depo-capsule",
        action="store_true",
        help="Add pairwise depolymerase × capsule cross-terms (GT01).",
    )
    replicate_parser.add_argument(
        "--include-pairwise-receptor-omp",
        action="store_true",
        help="Add directed receptor × OMP cross-terms (GT02).",
    )
    replicate_parser.add_argument(
        "--variant",
        choices=("all-pairs", "per-phage-blend"),
        default="all-pairs",
        help="Model variant: all-pairs or per-phage-blend (AX02).",
    )
    replicate_parser.add_argument(
        "--blend-alpha",
        type=float,
        default=0.5,
        help="Per-phage blend weight (AX02).",
    )
    replicate_parser.add_argument(
        "--calibrate",
        choices=("none", "isotonic"),
        default="none",
        help="Post-hoc calibration method (AX05).",
    )
    replicate_parser.add_argument(
        "--skip-track-g-prerequisites",
        action="store_true",
        help="Require all Track G lock artifacts and prerequisites to exist instead of regenerating them.",
    )
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_json(path: Path) -> dict[str, Any]:
    with path.open("r", encoding="utf-8") as handle:
        return json.load(handle)


def write_json(path: Path, payload: Mapping[str, Any]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_rows(path: Path, rows: Sequence[Mapping[str, object]]) -> None:
    if not rows:
        raise ValueError(f"Cannot write empty CSV artifact: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = list(rows[0].keys())
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


@contextmanager
def candidate_source_root(candidate_bundle: Path) -> Iterator[Path]:
    if not candidate_bundle.exists():
        raise FileNotFoundError(f"Candidate bundle does not exist: {candidate_bundle}")
    if candidate_bundle.is_dir():
        yield candidate_bundle
        return
    suffixes = candidate_bundle.suffixes
    if suffixes[-2:] != [".tar", ".gz"] and suffixes[-1:] != [".tgz"]:
        raise ValueError(f"Candidate bundle must be a directory or .tgz/.tar.gz archive: {candidate_bundle}")
    with tempfile.TemporaryDirectory(prefix="ar09_candidate_import_") as temp_dir:
        temp_root = Path(temp_dir)
        with tarfile.open(candidate_bundle, "r:gz") as archive:
            archive.extractall(temp_root, filter="data")
        yield temp_root


def find_candidate_bundle_root(root: Path) -> Path:
    if all((root / name).exists() for name in REQUIRED_CANDIDATE_FILES):
        return root
    matches = [path for path in root.rglob("train.py") if path.parent.is_dir()]
    for path in matches:
        candidate_root = path.parent
        if all((candidate_root / name).exists() for name in REQUIRED_CANDIDATE_FILES):
            return candidate_root
    raise FileNotFoundError(
        "Candidate bundle is missing one or more required files: " + ", ".join(REQUIRED_CANDIDATE_FILES)
    )


def derive_candidate_id(candidate_root: Path) -> str:
    local_metadata_path = candidate_root / "local_run_metadata.json"
    if local_metadata_path.exists():
        metadata = read_json(local_metadata_path)
        run_id = str(metadata.get("github_run_id", "")).strip()
        run_attempt = str(metadata.get("github_run_attempt", "")).strip()
        if run_id and run_attempt:
            return f"runpod_run_{run_id}_attempt_{run_attempt}"
        if run_id:
            return f"runpod_run_{run_id}"
    return candidate_root.name.replace(".", "_")


def copy_candidate_files(*, source_root: Path, destination_root: Path) -> dict[str, dict[str, object]]:
    copied: dict[str, dict[str, object]] = {}
    for name in (*REQUIRED_CANDIDATE_FILES, *OPTIONAL_CANDIDATE_FILES):
        source_path = source_root / name
        if not source_path.exists():
            if name in OPTIONAL_CANDIDATE_FILES:
                continue
            raise FileNotFoundError(f"Required candidate file missing from RunPod bundle: {source_path}")
        destination_path = destination_root / name
        destination_path.parent.mkdir(parents=True, exist_ok=True)
        shutil.copy2(source_path, destination_path)
        copied[name] = {
            "path": str(destination_path),
            "sha256": sha256_file(destination_path),
        }
    return copied


def import_runpod_candidate(
    *,
    candidate_bundle: Path,
    candidates_dir: Path,
    candidate_id: Optional[str],
) -> Path:
    with candidate_source_root(candidate_bundle) as extracted_root:
        source_root = find_candidate_bundle_root(extracted_root)
        resolved_candidate_id = candidate_id or derive_candidate_id(source_root)
        destination_root = candidates_dir / resolved_candidate_id
        if destination_root.exists():
            raise FileExistsError(f"Imported candidate directory already exists: {destination_root}")
        copied_files = copy_candidate_files(source_root=source_root, destination_root=destination_root)

    import_manifest = {
        "task_id": "AR09",
        "candidate_id": resolved_candidate_id,
        "imported_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "source_bundle": {
            "path": str(candidate_bundle),
            "sha256": sha256_file(candidate_bundle) if candidate_bundle.is_file() else None,
            "is_archive": candidate_bundle.is_file(),
        },
        "imported_files": copied_files,
    }
    summary_path = destination_root / "ar07_baseline_summary.json"
    if summary_path.exists():
        import_manifest["runpod_inner_val_summary"] = read_json(summary_path)

    write_json(destination_root / IMPORT_MANIFEST_FILENAME, import_manifest)
    LOGGER.info("Imported RunPod candidate into %s", destination_root)
    return destination_root


def ensure_autoresearch_cache(cache_dir: Path) -> Path:
    if cache_dir.exists():
        return cache_dir
    if cache_dir != runtime_contract.DEFAULT_CACHE_DIR:
        raise FileNotFoundError(f"Prepared AUTORESEARCH cache not found: {cache_dir}")
    LOGGER.info("AUTORESEARCH cache missing at %s; rebuilding with prepare.py", cache_dir)
    exit_code = prepare_cache.main([])
    if exit_code != 0:
        raise RuntimeError(f"prepare.py failed while rebuilding AUTORESEARCH cache with exit code {exit_code}")
    if not cache_dir.exists():
        raise FileNotFoundError(f"prepare.py did not create the expected cache directory: {cache_dir}")
    return cache_dir


def load_module_from_path(module_name: str, path: Path) -> ModuleType:
    import sys

    spec = importlib.util.spec_from_file_location(module_name, path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not import module from {path}")
    module = importlib.util.module_from_spec(spec)
    sys.modules[module_name] = module
    spec.loader.exec_module(module)
    return module


def load_holdout_frame(contract_manifest: Mapping[str, object]) -> pd.DataFrame:
    pair_table_path = Path(str(contract_manifest["pair_table"]["path"]))
    if not pair_table_path.exists():
        raise FileNotFoundError(f"AR01 pair table not found for sealed holdout replication: {pair_table_path}")
    pair_table = pd.read_csv(pair_table_path, dtype=str, keep_default_na=False)
    holdout_frame = pair_table.loc[
        (pair_table["split"] == build_contract.HOLDOUT_SPLIT) & (pair_table["retained_for_autoresearch"] == "1")
    ].copy()
    if holdout_frame.empty:
        raise ValueError("AR01 holdout retained set is empty.")
    if any(label not in {"0", "1"} for label in holdout_frame["label_any_lysis"]):
        raise ValueError("AR01 holdout contains non-binary retained labels.")
    return holdout_frame


def _load_st03_joined_frame(
    st02_path: Path = DEFAULT_ST02_PAIR_TABLE_PATH,
    st03_path: Path = DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH,
) -> pd.DataFrame:
    """Join ST02 pair table with ST03 split assignments on pair_id."""
    if not st02_path.exists():
        raise FileNotFoundError(f"ST02 pair table not found: {st02_path}")
    if not st03_path.exists():
        raise FileNotFoundError(f"ST03 split assignments not found: {st03_path}")
    st02 = pd.read_csv(st02_path, dtype=str, keep_default_na=False)
    st03 = pd.read_csv(st03_path, dtype=str, keep_default_na=False)
    joined = st02.merge(st03[["pair_id", "split_holdout", "is_hard_trainable"]], on="pair_id", validate="one_to_one")
    if len(joined) != len(st02):
        raise ValueError(f"ST02/ST03 join lost rows: {len(st02)} → {len(joined)}")
    return joined


def load_st03_holdout_frame(
    st02_path: Path = DEFAULT_ST02_PAIR_TABLE_PATH,
    st03_path: Path = DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH,
) -> pd.DataFrame:
    """Load the ST03 holdout test set: 65 bacteria, ~6235 trainable pairs."""
    joined = _load_st03_joined_frame(st02_path, st03_path)
    holdout = joined.loc[(joined["split_holdout"] == "holdout_test") & (joined["is_hard_trainable"] == "1")].copy()
    if holdout.empty:
        raise ValueError("ST03 holdout retained set is empty.")
    holdout = holdout.rename(columns={"label_hard_any_lysis": "label_any_lysis"})
    LOGGER.info("ST03 holdout: %d pairs, %d bacteria", len(holdout), holdout["bacteria"].nunique())
    return holdout


def build_st03_training_frame(
    st02_path: Path = DEFAULT_ST02_PAIR_TABLE_PATH,
    st03_path: Path = DEFAULT_ST03_SPLIT_ASSIGNMENTS_PATH,
) -> pd.DataFrame:
    """Load the ST03 training set: ~29031 trainable pairs."""
    joined = _load_st03_joined_frame(st02_path, st03_path)
    training = joined.loc[
        (joined["split_holdout"] == "train_non_holdout") & (joined["is_hard_trainable"] == "1")
    ].copy()
    if training.empty:
        raise ValueError("ST03 training retained set is empty.")
    training = training.rename(columns={"label_hard_any_lysis": "label_any_lysis"})
    LOGGER.info("ST03 training: %d pairs, %d bacteria", len(training), training["bacteria"].nunique())
    return training


def build_candidate_training_frame(context: Any) -> pd.DataFrame:
    frames = []
    for split_name in (runtime_contract.TRAIN_SPLIT, runtime_contract.INNER_VAL_SPLIT):
        split_frame = context.split_frames.get(split_name)
        if split_frame is None:
            raise ValueError(f"Candidate cache context is missing split {split_name}")
        frames.append(split_frame.loc[split_frame["retained_for_autoresearch"] == "1"].copy())
    training_frame = pd.concat(frames, ignore_index=True)
    if training_frame.empty:
        raise ValueError("Candidate replay training frame is empty after combining train + inner_val.")
    return training_frame


@contextmanager
def temporary_module_attribute(module: ModuleType, name: str, value: object) -> Iterator[None]:
    had_value = hasattr(module, name)
    previous = getattr(module, name, None)
    setattr(module, name, value)
    try:
        yield
    finally:
        if had_value:
            setattr(module, name, previous)
        else:
            delattr(module, name)


def build_candidate_holdout_rows(
    *,
    candidate_module: ModuleType,
    context: Any,
    holdout_frame: pd.DataFrame,
    seed: int,
    device_type: str,
    include_host_defense: bool,
    include_pairwise_depo_capsule: bool = False,
    include_pairwise_receptor_omp: bool = False,
    variant: str = "all-pairs",
    blend_alpha: float = 0.5,
    calibrate: str = "none",
    training_frame_override: pd.DataFrame | None = None,
) -> list[dict[str, object]]:
    host_slots = ["host_surface", "host_typing", "host_stats"]
    if include_host_defense:
        host_slots.append("host_defense")
    phage_slots = ["phage_projection", "phage_stats"]

    host_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts,
        slot_names=host_slots,
        entity_key="bacteria",
    )
    phage_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts,
        slot_names=phage_slots,
        entity_key="phage",
    )

    host_typed, host_numeric, host_categorical = candidate_module.type_entity_features(host_table, "bacteria")
    phage_typed, phage_numeric, phage_categorical = candidate_module.type_entity_features(phage_table, "phage")

    training_frame = (
        training_frame_override if training_frame_override is not None else build_candidate_training_frame(context)
    )
    train_design = candidate_module.build_raw_pair_design_matrix(
        training_frame, host_features=host_typed, phage_features=phage_typed
    )
    holdout_design = candidate_module.build_raw_pair_design_matrix(
        holdout_frame, host_features=host_typed, phage_features=phage_typed
    )

    # Pairwise depolymerase × capsule cross-terms (GT01).
    if include_pairwise_depo_capsule:
        from lyzortx.pipeline.autoresearch.derive_pairwise_depo_capsule_features import (
            compute_pairwise_depo_capsule_features,
        )

        compute_pairwise_depo_capsule_features(train_design)
        compute_pairwise_depo_capsule_features(holdout_design)

    # Directed receptor × OMP cross-terms (GT02).
    if include_pairwise_receptor_omp:
        from lyzortx.pipeline.autoresearch.derive_pairwise_receptor_omp_features import (
            compute_pairwise_receptor_omp_features,
        )

        compute_pairwise_receptor_omp_features(train_design)
        compute_pairwise_receptor_omp_features(holdout_design)

    all_slot_names = host_slots + phage_slots
    replay_prefixes = tuple(f"{s}__" for s in all_slot_names)
    if include_pairwise_depo_capsule:
        replay_prefixes = (*replay_prefixes, "pair_depo_capsule__")
    if include_pairwise_receptor_omp:
        replay_prefixes = (*replay_prefixes, "pair_receptor_omp__")
    feature_columns = [col for col in train_design.columns if col.startswith(replay_prefixes)]
    if not feature_columns:
        raise ValueError("Candidate replay constructed zero pair features for holdout replication.")
    LOGGER.info(
        "AR09 replay: %d feature columns from %d slots (prefixes: %s)",
        len(feature_columns),
        len(all_slot_names),
        ", ".join(replay_prefixes),
    )

    categorical_in_features = [col for col in (host_categorical + phage_categorical) if col in feature_columns]
    y_train = train_design["label_any_lysis"].astype(int).to_numpy(dtype=int)
    sample_weight = train_design["training_weight_v3"].astype(float).to_numpy(dtype=float)
    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[feature_columns],
        y_train,
        sample_weight=sample_weight,
        categorical_feature=categorical_in_features,
    )
    all_pairs_predictions = estimator.predict_proba(holdout_design[feature_columns])[:, 1]

    # Log feature importance by slot.
    imp = estimator.feature_importances_
    slot_imp: dict[str, float] = {}
    for col, val in zip(feature_columns, imp):
        slot = col.split("__")[0]
        slot_imp[slot] = slot_imp.get(slot, 0) + val
    total_imp = sum(slot_imp.values()) or 1
    parts = [f"{s}={v / total_imp * 100:.1f}%" for s, v in sorted(slot_imp.items(), key=lambda x: -x[1])]
    LOGGER.info("AR09 feature importance by slot: %s", ", ".join(parts))

    # Per-phage blending (AX02).
    if variant == "per-phage-blend":
        from lyzortx.autoresearch.per_phage_model import fit_per_phage_models, predict_per_phage

        host_feature_columns = [
            col
            for col in feature_columns
            if col.startswith(("host_surface__", "host_typing__", "host_stats__", "host_defense__"))
        ]
        per_phage_models = fit_per_phage_models(
            train_design, host_feature_columns, device_type=device_type, random_state=seed
        )
        predictions, _ = predict_per_phage(
            per_phage_models, holdout_design, host_feature_columns, all_pairs_predictions, blend_alpha=blend_alpha
        )
    else:
        predictions = all_pairs_predictions

    # Isotonic calibration (AX05).
    if calibrate == "isotonic":
        from lyzortx.autoresearch.calibration import calibrate_predictions, fit_isotonic_calibrator_cv

        calibrator = fit_isotonic_calibrator_cv(
            train_design,
            feature_columns,
            y_train,
            sample_weight,
            categorical_features=categorical_in_features,
            device_type=device_type,
            model_params=dict(candidate_module.PAIR_SCORER_PARAMS),
        )
        predictions = calibrate_predictions(calibrator, predictions)

    rows = []
    for row, probability in zip(
        holdout_design.loc[:, ["pair_id", "bacteria", "phage", "label_any_lysis"]].to_dict(orient="records"),
        predictions,
    ):
        rows.append(
            {
                "arm_id": CANDIDATE_ARM_ID,
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(row["label_any_lysis"]),
                "predicted_probability": safe_round(float(probability)),
            }
        )
    return rows


def resolve_feature_lock_path(path: Path, expected_locked_blocks: Sequence[object] = ()) -> Path:
    if path.exists():
        # Check whether the generated lock matches the contract before using it.
        if expected_locked_blocks:
            try:
                payload = json.loads(path.read_text(encoding="utf-8"))
                actual = [str(b) for b in payload.get("winner_subset_blocks", [])]
                expected = [str(b) for b in expected_locked_blocks]
                if actual != expected and FALLBACK_FEATURE_LOCK_PATH.exists():
                    LOGGER.warning(
                        "Generated feature lock at %s has blocks %r (expected %r); using fallback %s",
                        path,
                        actual,
                        expected,
                        FALLBACK_FEATURE_LOCK_PATH,
                    )
                    return FALLBACK_FEATURE_LOCK_PATH
            except (json.JSONDecodeError, KeyError):
                pass
        return path
    if path == DEFAULT_FEATURE_LOCK_PATH and FALLBACK_FEATURE_LOCK_PATH.exists():
        return FALLBACK_FEATURE_LOCK_PATH
    raise FileNotFoundError(f"Locked comparator feature-lock artifact not found: {path}")


def load_comparator_params(path: Path) -> dict[str, object]:
    if path.exists():
        return load_tg01_lock(path)["best_params"]
    if path == DEFAULT_TG01_MODEL_SUMMARY_PATH:
        return dict(FALLBACK_TG01_BEST_PARAMS)
    raise FileNotFoundError(f"Locked comparator model summary not found: {path}")


def validate_comparator_feature_lock(path: Path, expected_locked_blocks: Sequence[object]) -> dict[str, object]:
    payload = load_v1_lock(path)
    expected_blocks = [str(block) for block in expected_locked_blocks]
    actual_blocks = [str(block) for block in payload.get("winner_subset_blocks", [])]
    if expected_blocks and actual_blocks != expected_blocks:
        raise ValueError(
            "Comparator feature lock does not match the AR01 contract: "
            f"expected {expected_blocks!r}, found {actual_blocks!r} in {path}"
        )
    return payload


def ensure_track_g_prerequisites(skip_prerequisites: bool) -> argparse.Namespace:
    args = train_v1_binary_classifier.parse_args([])
    required_paths = (
        args.st02_pair_table_path,
        args.track_c_pair_table_path,
        args.track_d_genome_kmer_path,
        args.track_d_distance_path,
    )
    if all(path.exists() for path in required_paths):
        return args
    if skip_prerequisites:
        missing = ", ".join(str(path) for path in required_paths if not path.exists())
        raise FileNotFoundError(f"Missing Track G prerequisite artifact(s): {missing}")
    LOGGER.info("Track G prerequisite artifacts missing; regenerating prerequisite outputs")
    train_v1_binary_classifier.ensure_prerequisite_outputs(args)
    return args


def build_comparator_feature_space(
    *,
    st02_rows: Sequence[Mapping[str, str]],
    track_c_rows: Sequence[Mapping[str, str]],
    track_d_genome_rows: Sequence[Mapping[str, str]],
    track_d_distance_rows: Sequence[Mapping[str, str]],
) -> train_v1_binary_classifier.FeatureSpace:
    st02_columns = tuple(st02_rows[0].keys())
    track_c_new_columns = [column for column in track_c_rows[0].keys() if column not in st02_columns]
    defense_columns = partition_track_c_defense_columns(track_c_new_columns)
    track_d_feature_columns = train_v1_binary_classifier.deduplicate_preserving_order(
        [column for column in track_d_genome_rows[0].keys() if column != "phage"]
        + [column for column in track_d_distance_rows[0].keys() if column != "phage"]
    )
    numeric_columns = train_v1_binary_classifier.deduplicate_preserving_order(
        list(train_v1_binary_classifier.V0_NUMERIC_FEATURE_COLUMNS)
        + list(defense_columns)
        + list(track_d_feature_columns)
    )
    return train_v1_binary_classifier.FeatureSpace(
        categorical_columns=tuple(train_v1_binary_classifier.V0_CATEGORICAL_FEATURE_COLUMNS),
        numeric_columns=numeric_columns,
        track_c_additional_columns=tuple(track_c_new_columns),
        track_d_columns=tuple(track_d_feature_columns),
        track_e_columns=(),
    )


def build_ar01_comparator_rows(
    *,
    ar01_rows: Sequence[Mapping[str, str]],
    st02_rows: Sequence[Mapping[str, str]],
    track_c_rows: Sequence[Mapping[str, str]],
    track_d_genome_rows: Sequence[Mapping[str, str]],
    track_d_distance_rows: Sequence[Mapping[str, str]],
) -> list[dict[str, object]]:
    st02_by_pair = {row["pair_id"]: dict(row) for row in st02_rows}
    track_c_by_pair = {row["pair_id"]: dict(row) for row in track_c_rows}
    track_d_genome_by_phage = {row["phage"]: dict(row) for row in track_d_genome_rows}
    track_d_distance_by_phage = {row["phage"]: dict(row) for row in track_d_distance_rows}
    st02_columns = set(st02_rows[0].keys())

    merged_rows = []
    for ar01_row in ar01_rows:
        pair_id = str(ar01_row["pair_id"])
        st02_row = st02_by_pair.get(pair_id)
        if st02_row is None:
            raise KeyError(f"ST02 is missing pair_id {pair_id} needed for comparator replay.")
        track_c_row = track_c_by_pair.get(pair_id)
        if track_c_row is None:
            raise KeyError(f"Track C pair table is missing pair_id {pair_id} needed for comparator replay.")
        phage = str(ar01_row["phage"])
        track_d_genome_row = track_d_genome_by_phage.get(phage)
        track_d_distance_row = track_d_distance_by_phage.get(phage)
        if track_d_genome_row is None or track_d_distance_row is None:
            raise KeyError(f"Track D feature blocks are missing phage {phage} needed for comparator replay.")

        merged = dict(st02_row)
        for column, value in track_c_row.items():
            if column not in st02_columns:
                merged[column] = value
        for feature_row in (track_d_genome_row, track_d_distance_row):
            for column, value in feature_row.items():
                if column != "phage":
                    merged[column] = value

        split_name = str(ar01_row["split"])
        merged.update(
            {
                "pair_id": pair_id,
                "bacteria": str(ar01_row["bacteria"]),
                "phage": phage,
                "label_hard_any_lysis": str(ar01_row["label_any_lysis"]),
                "training_weight_v3": str(ar01_row["training_weight_v3"]),
                "is_hard_trainable": "1",
                "split_holdout": "holdout_test" if split_name == build_contract.HOLDOUT_SPLIT else "train_non_holdout",
                "split_cv5_fold": "-1",
            }
        )
        merged_rows.append(merged)
    return merged_rows


def build_comparator_holdout_rows(
    *,
    contract_manifest: Mapping[str, object],
    skip_prerequisites: bool,
    seed: int,
) -> list[dict[str, object]]:
    comparator_manifest = dict(contract_manifest["current_locked_comparator_benchmark"])
    expected_blocks = comparator_manifest.get("locked_feature_blocks", ())
    feature_lock_path = resolve_feature_lock_path(
        Path(str(comparator_manifest["feature_lock_path"])), expected_locked_blocks=expected_blocks
    )
    # AR09 replays the locked comparator at the block level, not via a hand-maintained per-column allowlist.
    # The v1 lock records the winning subset blocks ("defense" + "phage_genomic"), so we validate that contract
    # here and then rebuild the full defense/phage-genomic feature space from the raw Track C/Track D columns.
    validate_comparator_feature_lock(feature_lock_path, expected_blocks)
    comparator_params = load_comparator_params(Path(str(comparator_manifest["model_summary_path"])))
    track_g_args = ensure_track_g_prerequisites(skip_prerequisites)

    st02_rows = read_csv_rows(track_g_args.st02_pair_table_path)
    track_c_rows = read_csv_rows(track_g_args.track_c_pair_table_path)
    track_d_genome_rows = read_csv_rows(track_g_args.track_d_genome_kmer_path)
    track_d_distance_rows = read_csv_rows(track_g_args.track_d_distance_path)
    ar01_pair_table_path = Path(str(contract_manifest["pair_table"]["path"]))
    ar01_rows = [
        row
        for row in read_csv_rows(ar01_pair_table_path)
        if row["retained_for_autoresearch"] == "1" and row["label_any_lysis"] in {"0", "1"}
    ]
    comparator_rows = build_ar01_comparator_rows(
        ar01_rows=ar01_rows,
        st02_rows=st02_rows,
        track_c_rows=track_c_rows,
        track_d_genome_rows=track_d_genome_rows,
        track_d_distance_rows=track_d_distance_rows,
    )
    feature_space = build_comparator_feature_space(
        st02_rows=st02_rows,
        track_c_rows=track_c_rows,
        track_d_genome_rows=track_d_genome_rows,
        track_d_distance_rows=track_d_distance_rows,
    )

    estimator_factory = lambda params, seed_offset: train_v1_binary_classifier.make_lightgbm_estimator(  # noqa: E731
        params,
        seed_offset,
        base_random_state=seed,
    )
    _, _, _, eval_rows, probabilities = train_v1_binary_classifier.fit_final_estimator(
        comparator_rows,
        feature_space,
        estimator_factory=estimator_factory,
        params=comparator_params,
        sample_weight_key="training_weight_v3",
    )

    rows = []
    for row, probability in zip(eval_rows, probabilities):
        rows.append(
            {
                "arm_id": BASELINE_ARM_ID,
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "label_hard_any_lysis": int(str(row["label_hard_any_lysis"])),
                "predicted_probability": safe_round(float(probability)),
            }
        )
    return rows


def summarize_seed_metrics(rows: Sequence[Mapping[str, object]]) -> dict[str, Optional[float]]:
    y_true = [int(row["label_hard_any_lysis"]) for row in rows]
    y_prob = [float(row["predicted_probability"]) for row in rows]
    binary = train_v1_binary_classifier.compute_binary_metrics(y_true, y_prob)
    return {
        "holdout_roc_auc": binary["roc_auc"],
        "holdout_brier_score": binary["brier_score"],
    }


def aggregate_seed_rows(rows: Sequence[Mapping[str, object]]) -> list[dict[str, object]]:
    frame = pd.DataFrame(rows)
    if frame.empty:
        raise ValueError("Cannot aggregate empty repeated-seed holdout rows.")
    grouped = (
        frame.groupby(["arm_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
            "predicted_probability"
        ]
        .mean()
        .sort_values(["arm_id", "bacteria", "phage"])
    )
    return [
        {
            "arm_id": str(row["arm_id"]),
            "pair_id": str(row["pair_id"]),
            "bacteria": str(row["bacteria"]),
            "phage": str(row["phage"]),
            "label_hard_any_lysis": int(row["label_hard_any_lysis"]),
            "predicted_probability": safe_round(float(row["predicted_probability"])),
        }
        for row in grouped.to_dict(orient="records")
    ]


def summarize_arm(rows: Sequence[Mapping[str, object]]) -> dict[str, object]:
    metrics = summarize_seed_metrics(rows)
    return {
        "holdout_roc_auc": metrics["holdout_roc_auc"],
        "holdout_brier_score": metrics["holdout_brier_score"],
    }


def bootstrap_metric_block(metric_ci: Any) -> dict[str, object]:
    return {
        "point_estimate": metric_ci.point_estimate,
        "ci_low": metric_ci.ci_low,
        "ci_high": metric_ci.ci_high,
        "bootstrap_samples_requested": metric_ci.bootstrap_samples_requested,
        "bootstrap_samples_used": metric_ci.bootstrap_samples_used,
    }


def build_decision_summary(
    *,
    candidate_rows: Sequence[Mapping[str, object]],
    comparator_rows: Sequence[Mapping[str, object]],
    bootstrap_summary: Mapping[str, Mapping[str, Any]],
) -> dict[str, object]:
    candidate_metrics = summarize_arm(candidate_rows)
    comparator_metrics = summarize_arm(comparator_rows)
    delta_key = f"{CANDIDATE_ARM_ID}__delta_vs_{BASELINE_ARM_ID}"
    delta_summary = bootstrap_summary[delta_key]

    promotion_rule = {
        "primary_metric": PRIMARY_METRIC,
        "promotion_requires": {
            "auc_delta_ci_low_gt_zero": True,
            "brier_improvement_ci_high_gte_zero": True,
        },
    }
    auc_clear = delta_summary["holdout_roc_auc"].ci_low is not None and delta_summary["holdout_roc_auc"].ci_low > 0.0
    brier_clear = (
        delta_summary["holdout_brier_score"].ci_high is not None and delta_summary["holdout_brier_score"].ci_high >= 0.0
    )
    if auc_clear and brier_clear:
        decision = "promote"
        rationale = "candidate clears the predeclared AUC lift rule without material Brier regression"
    else:
        decision = "no_honest_lift"
        reasons = []
        if not auc_clear:
            reasons.append("AUC delta stays within bootstrap noise")
        if not brier_clear:
            reasons.append("Brier score materially degrades")
        rationale = "; ".join(reasons)

    return {
        "promotion_rule": promotion_rule,
        "decision": decision,
        "decision_rationale": rationale,
        "arm_metrics": {
            CANDIDATE_ARM_ID: candidate_metrics,
            BASELINE_ARM_ID: comparator_metrics,
        },
        "bootstrap_summary": {
            arm_id: {metric_name: bootstrap_metric_block(metric_ci) for metric_name, metric_ci in metric_map.items()}
            for arm_id, metric_map in bootstrap_summary.items()
        },
    }


def replicate_candidate(args: argparse.Namespace) -> Path:
    cache_dir = ensure_autoresearch_cache(args.cache_dir)
    candidate_manifest_path = args.candidate_dir / IMPORT_MANIFEST_FILENAME
    if not candidate_manifest_path.exists():
        raise FileNotFoundError(f"Imported candidate manifest not found: {candidate_manifest_path}")
    candidate_manifest = read_json(candidate_manifest_path)
    candidate_id = str(candidate_manifest["candidate_id"])
    output_dir = args.output_dir or (DEFAULT_DECISION_BUNDLES_DIR / candidate_id)

    candidate_train_path = args.candidate_dir / "train.py"
    candidate_module = load_module_from_path(f"ar09_candidate_{candidate_id}", candidate_train_path)
    context = candidate_module.load_and_validate_cache(
        cache_dir=cache_dir,
        include_host_defense=args.include_host_defense,
    )

    use_st03 = getattr(args, "use_st03_split", False)
    if use_st03:
        holdout_frame = load_st03_holdout_frame()
        training_frame_override = build_st03_training_frame()
        LOGGER.info("AR09 using ST03 split (65-bacteria holdout, same as TL18)")
    else:
        holdout_frame = load_holdout_frame(context.contract_manifest)
        training_frame_override = None

    all_seed_rows: list[dict[str, object]] = []
    seed_metric_rows: list[dict[str, object]] = []
    for seed in args.replication_seeds:
        LOGGER.info("AR09 replay seed %d: candidate", seed)
        candidate_seed_rows = build_candidate_holdout_rows(
            candidate_module=candidate_module,
            context=context,
            holdout_frame=holdout_frame,
            seed=seed,
            device_type=args.device_type,
            include_host_defense=args.include_host_defense,
            include_pairwise_depo_capsule=getattr(args, "include_pairwise_depo_capsule", False),
            include_pairwise_receptor_omp=getattr(args, "include_pairwise_receptor_omp", False),
            variant=getattr(args, "variant", "all-pairs"),
            blend_alpha=getattr(args, "blend_alpha", 0.5),
            calibrate=getattr(args, "calibrate", "none"),
            training_frame_override=training_frame_override,
        )
        all_seed_rows.extend(candidate_seed_rows)
        candidate_seed_metrics = summarize_seed_metrics(candidate_seed_rows)
        seed_metric_rows.append({"arm_id": CANDIDATE_ARM_ID, "seed": seed, **candidate_seed_metrics})

        if not use_st03:
            LOGGER.info("AR09 replay seed %d: locked comparator", seed)
            comparator_seed_rows = build_comparator_holdout_rows(
                contract_manifest=context.contract_manifest,
                skip_prerequisites=args.skip_track_g_prerequisites,
                seed=seed,
            )
            all_seed_rows.extend(comparator_seed_rows)
            comparator_seed_metrics = summarize_seed_metrics(comparator_seed_rows)
            seed_metric_rows.append({"arm_id": BASELINE_ARM_ID, "seed": seed, **comparator_seed_metrics})

    aggregated_rows = aggregate_seed_rows(all_seed_rows)
    candidate_rows = [row for row in aggregated_rows if row["arm_id"] == CANDIDATE_ARM_ID]

    if use_st03:
        # Candidate-only bootstrap CIs (no comparator arm)
        bootstrap_summary = bootstrap_holdout_metric_cis(
            {CANDIDATE_ARM_ID: candidate_rows},
            bootstrap_samples=args.bootstrap_samples,
            bootstrap_random_state=args.bootstrap_random_state,
            baseline_arm_id=CANDIDATE_ARM_ID,
        )
        candidate_metrics = summarize_arm(candidate_rows)
        decision_summary = {
            "evaluation_mode": "st03_candidate_only",
            "arm_metrics": {CANDIDATE_ARM_ID: candidate_metrics},
            "bootstrap_summary": {
                arm_id: {
                    metric_name: bootstrap_metric_block(metric_ci) for metric_name, metric_ci in metric_map.items()
                }
                for arm_id, metric_map in bootstrap_summary.items()
            },
        }
    else:
        comparator_rows = [row for row in aggregated_rows if row["arm_id"] == BASELINE_ARM_ID]
        bootstrap_summary = bootstrap_holdout_metric_cis(
            {
                BASELINE_ARM_ID: comparator_rows,
                CANDIDATE_ARM_ID: candidate_rows,
            },
            bootstrap_samples=args.bootstrap_samples,
            bootstrap_random_state=args.bootstrap_random_state,
            baseline_arm_id=BASELINE_ARM_ID,
        )
        decision_summary = build_decision_summary(
            candidate_rows=candidate_rows,
            comparator_rows=comparator_rows,
            bootstrap_summary=bootstrap_summary,
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    write_rows(output_dir / AGGREGATED_PREDICTIONS_FILENAME, aggregated_rows)
    write_rows(output_dir / SEED_METRICS_FILENAME, seed_metric_rows)
    bundle = {
        "task_id": "AR09",
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "candidate_id": candidate_id,
        "candidate_provenance": candidate_manifest,
        "replication_inputs": {
            "cache_dir": str(cache_dir),
            "contract_manifest_path": str(
                context.provenance_manifest["source_contract"]["pair_contract_manifest_path"]
            ),
            "replication_seeds": list(args.replication_seeds),
            "device_type": args.device_type,
            "include_host_defense": args.include_host_defense,
            "include_pairwise_depo_capsule": getattr(args, "include_pairwise_depo_capsule", False),
            "include_pairwise_receptor_omp": getattr(args, "include_pairwise_receptor_omp", False),
            "variant": getattr(args, "variant", "all-pairs"),
            "blend_alpha": getattr(args, "blend_alpha", 0.5),
            "calibrate": getattr(args, "calibrate", "none"),
            "use_st03_split": use_st03,
        },
        "artifacts": {
            "aggregated_holdout_predictions_path": str(output_dir / AGGREGATED_PREDICTIONS_FILENAME),
            "seed_metrics_path": str(output_dir / SEED_METRICS_FILENAME),
        },
        **decision_summary,
    }
    write_json(output_dir / DECISION_BUNDLE_FILENAME, bundle)
    LOGGER.info("AR09 wrote decision bundle to %s", output_dir)
    return output_dir


def main(argv: Optional[Sequence[str]] = None) -> int:
    setup_logging()
    args = parse_args(argv)
    if args.command == "import-runpod-candidate":
        import_runpod_candidate(
            candidate_bundle=args.candidate_bundle,
            candidates_dir=args.candidates_dir,
            candidate_id=args.candidate_id,
        )
        return 0
    if args.command == "replicate":
        replicate_candidate(args)
        return 0
    raise ValueError(f"Unsupported command: {args.command}")


if __name__ == "__main__":
    raise SystemExit(main())
