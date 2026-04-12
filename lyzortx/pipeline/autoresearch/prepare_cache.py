#!/usr/bin/env python3
"""AR02: scaffold the AUTORESEARCH sandbox and freeze the search-cache contract."""

from __future__ import annotations

import argparse
import csv
import json
import logging
import os
from collections import Counter
from concurrent.futures import ProcessPoolExecutor, as_completed
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Mapping, Optional, Sequence

import joblib
import numpy as np

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch import build_contract
from lyzortx.pipeline.autoresearch import derive_phage_stats_features
from lyzortx.pipeline.autoresearch import runtime_contract
from lyzortx.pipeline.deployment_paired_features import derive_host_stats_features
from lyzortx.pipeline.deployment_paired_features import derive_host_typing_features
from lyzortx.pipeline.deployment_paired_features import run_all_host_surface
from lyzortx.pipeline.deployment_paired_features.derive_host_defense_features import (
    DEFAULT_PANEL_DEFENSE_SUBTYPES_PATH,
    PER_HOST_COUNTS_FILENAME,
)
from lyzortx.pipeline.deployment_paired_features.run_all_host_defense import (
    aggregate_host_defense_csvs,
)
from lyzortx.pipeline.steel_thread_v0.io.write_outputs import ensure_directory, write_csv, write_json
from lyzortx.pipeline.track_l.steps.run_novel_host_defense_finder import (
    DEFAULT_MODELS_DIR,
    MODEL_INSTALL_MODE_FORBID,
    ensure_defense_finder_models,
)
from lyzortx.pipeline.track_l.steps import build_tl17_phage_compatibility_preprocessor as tl17_preprocessor
from lyzortx.pipeline.track_l.steps import deployable_tl17_runtime as tl17_runtime

LOGGER = logging.getLogger(__name__)

DEFAULT_OUTPUT_ROOT = runtime_contract.DEFAULT_OUTPUT_ROOT
DEFAULT_CACHE_DIR = runtime_contract.DEFAULT_CACHE_DIR

CACHE_MANIFEST_FILENAME = runtime_contract.CACHE_MANIFEST_FILENAME
SCHEMA_MANIFEST_FILENAME = runtime_contract.SCHEMA_MANIFEST_FILENAME
PROVENANCE_MANIFEST_FILENAME = runtime_contract.PROVENANCE_MANIFEST_FILENAME
TRAIN_PAIR_TABLE_FILENAME = runtime_contract.TRAIN_PAIR_TABLE_FILENAME
INNER_VAL_PAIR_TABLE_FILENAME = runtime_contract.INNER_VAL_PAIR_TABLE_FILENAME
ENTITY_INDEX_FILENAME = runtime_contract.ENTITY_INDEX_FILENAME
SLOT_SCHEMA_FILENAME = runtime_contract.SLOT_SCHEMA_FILENAME
SLOT_FEATURE_TABLE_FILENAME = runtime_contract.SLOT_FEATURE_TABLE_FILENAME
SLOT_FEATURES_FILENAME = runtime_contract.SLOT_FEATURES_FILENAME
HOST_DEFENSE_BUILD_MANIFEST_FILENAME = "host_defense_build_manifest.json"
HOST_TYPING_BUILD_MANIFEST_FILENAME = "host_typing_build_manifest.json"
HOST_STATS_BUILD_MANIFEST_FILENAME = "host_stats_build_manifest.json"
PHAGE_PROJECTION_BUILD_MANIFEST_FILENAME = "phage_projection_build_manifest.json"
PHAGE_STATS_BUILD_MANIFEST_FILENAME = "phage_stats_build_manifest.json"
PHAGE_KMER_BUILD_MANIFEST_FILENAME = "phage_kmer_build_manifest.json"
PHAGE_KMER_K = 4
PHAGE_KMER_DIM = 4**PHAGE_KMER_K  # 256
DEFAULT_PHAROKKA_ANNOTATION_DIR = Path("data/annotations/pharokka")

TASK_ID = "AR02"
CACHE_CONTRACT_ID = runtime_contract.CACHE_CONTRACT_ID
SCHEMA_MANIFEST_ID = runtime_contract.SCHEMA_MANIFEST_ID
SEARCH_PAIR_TABLE_ID = runtime_contract.SEARCH_PAIR_TABLE_ID
HOLDOUT_HANDLING_RULE = runtime_contract.HOLDOUT_HANDLING_RULE
SUPPORTED_SEARCH_SPLITS = runtime_contract.SUPPORTED_SEARCH_SPLITS
DISALLOWED_SEARCH_SPLITS = runtime_contract.DISALLOWED_SEARCH_SPLITS
PAIR_KEY = runtime_contract.PAIR_KEY
HOST_SURFACE_BUILD_DIRNAME = "host_surface_cache_build"
HOST_TYPING_BUILD_DIRNAME = "host_typing_cache_build"
HOST_STATS_BUILD_DIRNAME = "host_stats_cache_build"
PHAGE_PROJECTION_BUILD_DIRNAME = "phage_projection_cache_build"


SlotSpec = runtime_contract.SlotSpec
SLOT_SPECS = runtime_contract.SLOT_SPECS
SLOT_SPEC_BY_NAME = runtime_contract.SLOT_SPEC_BY_NAME


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-interactions-path",
        type=Path,
        default=build_contract.DEFAULT_RAW_INTERACTIONS_PATH,
        help="Semicolon-delimited raw interaction table.",
    )
    parser.add_argument(
        "--host-assembly-dir",
        type=Path,
        default=build_contract.DEFAULT_ASSEMBLY_DIR,
        help="Directory containing Picard host FASTAs.",
    )
    parser.add_argument(
        "--phage-fasta-dir",
        type=Path,
        default=build_contract.DEFAULT_PHAGE_FASTA_DIR,
        help="Directory containing phage FASTA files.",
    )
    parser.add_argument(
        "--output-root",
        type=Path,
        default=DEFAULT_OUTPUT_ROOT,
        help="AUTORESEARCH generated-output root directory.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=DEFAULT_CACHE_DIR,
        help="Directory where the search cache should be written.",
    )
    parser.add_argument(
        "--warm-cache-manifest-path",
        type=Path,
        default=None,
        help="Optional manifest describing precomputed warm-cache accelerators.",
    )
    parser.add_argument(
        "--holdout-fraction",
        type=float,
        default=build_contract.DEFAULT_HOLDOUT_FRACTION,
        help="Fraction of bacteria routed to sealed holdout.",
    )
    parser.add_argument(
        "--inner-val-fraction",
        type=float,
        default=build_contract.DEFAULT_INNER_VAL_FRACTION,
        help="Fraction of bacteria routed to inner validation.",
    )
    parser.add_argument(
        "--split-salt",
        default=build_contract.DEFAULT_SPLIT_SALT,
        help="Deterministic salt for bacteria split assignment.",
    )
    parser.add_argument(
        "--skip-host-assembly-resolution",
        action="store_true",
        help="Skip download_picard_assemblies() and trust the provided host assembly directory as-is.",
    )
    parser.add_argument(
        "--skip-comparator-lock",
        action="store_true",
        help="Skip comparator artifact checksum validation (for CI/RunPod where Track G outputs are unavailable).",
    )
    parser.add_argument(
        "--host-defense-output-dir",
        type=Path,
        default=None,
        help="Per-host raw host-defense cache root. Defaults to <output-root>/host_defense.",
    )
    parser.add_argument(
        "--host-defense-models-dir",
        type=Path,
        default=DEFAULT_MODELS_DIR,
        help="Pinned Defense Finder models directory shared across host-defense workers.",
    )
    parser.add_argument(
        "--host-defense-max-workers",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Process fan-out for host-defense cache building.",
    )
    parser.add_argument(
        "--host-defense-force-model-update",
        action="store_true",
        help="Force a one-time coordinator-side model refresh before worker fan-out.",
    )
    parser.add_argument(
        "--host-defense-force-run",
        action="store_true",
        help="Re-run Defense Finder even when a per-host systems TSV already exists.",
    )
    parser.add_argument(
        "--host-defense-preserve-raw",
        action="store_true",
        help="Preserve raw MacSyFinder outputs from Defense Finder.",
    )
    parser.add_argument(
        "--host-defense-aggregate-only",
        action="store_true",
        help="Skip host-defense worker execution and rebuild the slot artifact from existing per-host outputs only.",
    )
    parser.add_argument(
        "--host-surface-build-dir",
        type=Path,
        default=None,
        help="Optional directory for cached host-surface intermediates; defaults under the output root.",
    )
    parser.add_argument(
        "--host-surface-max-workers",
        type=int,
        default=max(1, os.cpu_count() or 1),
        help="Worker count for the host-surface pyhmmer fast path.",
    )
    parser.add_argument(
        "--tl17-output-dir",
        type=Path,
        default=tl17_preprocessor.DEFAULT_OUTPUT_DIR,
        help="Directory containing the frozen TL17 runtime payload and reference bank.",
    )
    parser.add_argument(
        "--pharokka-annotation-dir",
        type=Path,
        default=DEFAULT_PHAROKKA_ANNOTATION_DIR,
        help="Directory containing Pharokka merged CDS annotation TSVs per phage.",
    )
    return parser.parse_args(argv)


def load_csv_rows(path: Path) -> list[dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def read_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def run_ar01_contract(args: argparse.Namespace) -> tuple[Path, Path, Path]:
    LOGGER.info("AR02 starting: refreshing AR01 raw-input contract under %s", args.output_root)
    exit_code = build_contract.main(
        [
            "--raw-interactions-path",
            str(args.raw_interactions_path),
            "--host-assembly-dir",
            str(args.host_assembly_dir),
            "--phage-fasta-dir",
            str(args.phage_fasta_dir),
            "--output-dir",
            str(args.output_root),
            "--holdout-fraction",
            str(args.holdout_fraction),
            "--inner-val-fraction",
            str(args.inner_val_fraction),
            "--split-salt",
            args.split_salt,
            *(["--skip-host-assembly-resolution"] if args.skip_host_assembly_resolution else []),
            *(["--skip-comparator-lock"] if args.skip_comparator_lock else []),
        ]
    )
    if exit_code != 0:
        raise RuntimeError(f"AR01 contract build failed with exit code {exit_code}")

    pair_table_path = args.output_root / build_contract.PAIR_TABLE_FILENAME
    contract_manifest_path = args.output_root / build_contract.CONTRACT_MANIFEST_FILENAME
    input_checksums_path = args.output_root / build_contract.INPUT_CHECKSUMS_FILENAME
    for path in (pair_table_path, contract_manifest_path, input_checksums_path):
        if not path.exists():
            raise FileNotFoundError(f"Expected AR01 artifact was not written: {path}")
    return pair_table_path, contract_manifest_path, input_checksums_path


def select_search_rows(pair_rows: Sequence[Mapping[str, str]]) -> dict[str, list[dict[str, str]]]:
    selected: dict[str, list[dict[str, str]]] = {split: [] for split in SUPPORTED_SEARCH_SPLITS}
    seen_disallowed_splits: set[str] = set()

    for row in pair_rows:
        split_name = str(row["split"])
        if split_name in SUPPORTED_SEARCH_SPLITS:
            selected[split_name].append(dict(row))
        elif split_name in DISALLOWED_SEARCH_SPLITS:
            seen_disallowed_splits.add(split_name)
        else:
            raise ValueError(f"Unexpected split in AR01 pair table: {split_name}")

    for split_name in SUPPORTED_SEARCH_SPLITS:
        if not selected[split_name]:
            raise ValueError(f"AR02 search cache would have an empty required split: {split_name}")
    if not seen_disallowed_splits:
        raise ValueError("AR01 pair table did not contain the sealed holdout split expected by AR02.")
    return selected


def select_all_retained_rows(pair_rows: Sequence[Mapping[str, str]]) -> dict[str, list[dict[str, str]]]:
    """Select rows from all splits (including holdout) for feature materialization.

    Feature slots must cover holdout entities so that AR09 replication can compute
    holdout embeddings. Only *labels* are sealed — entity-level features derived from
    raw FASTAs carry no label information.
    """
    all_splits = (*SUPPORTED_SEARCH_SPLITS, *DISALLOWED_SEARCH_SPLITS)
    selected: dict[str, list[dict[str, str]]] = {split: [] for split in all_splits}
    for row in pair_rows:
        split_name = str(row["split"])
        if split_name in all_splits:
            selected[split_name].append(dict(row))
    return selected


def build_slot_index_rows(
    *,
    slot_spec: SlotSpec,
    selected_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> list[dict[str, str]]:
    values = {
        str(row[slot_spec.entity_key])
        for split_rows in selected_rows.values()
        for row in split_rows
        if str(row["retained_for_autoresearch"]) == "1"
    }
    if not values:
        raise ValueError(f"Reserved slot {slot_spec.slot_name} would have zero retained entities.")
    return [{slot_spec.entity_key: value} for value in sorted(values)]


def build_slot_schema_manifest(
    slot_spec: SlotSpec,
    row_count: int,
    *,
    reserved_feature_columns: Sequence[str] = (),
) -> dict[str, Any]:
    resolved_columns = [str(column) for column in reserved_feature_columns]
    return {
        "task_id": TASK_ID,
        "schema_manifest_id": SCHEMA_MANIFEST_ID,
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "entity_key": slot_spec.entity_key,
        "join_keys": slot_spec.join_keys,
        "column_family_prefix": slot_spec.column_prefix,
        "block_role": slot_spec.block_role,
        "reserved_feature_columns": resolved_columns,
        "reserved_feature_column_count": len(resolved_columns),
        "entity_index_row_count": row_count,
        "composability_contract": {
            "join_type": "left",
            "row_granularity": f"one_row_per_{slot_spec.entity_key}",
            "column_ownership": (
                f"Future columns for {slot_spec.slot_name} must start with {slot_spec.column_prefix} "
                f"and may only be added inside this slot."
            ),
        },
        "description": slot_spec.description,
    }


def build_top_level_schema_manifest(
    *,
    slot_feature_columns: Mapping[str, Sequence[str]] | None = None,
) -> dict[str, Any]:
    resolved_columns = {} if slot_feature_columns is None else slot_feature_columns
    resolved_slot_entries = {
        spec.slot_name: [str(column) for column in resolved_columns.get(spec.slot_name, ())] for spec in SLOT_SPECS
    }
    return {
        "task_id": TASK_ID,
        "schema_manifest_id": SCHEMA_MANIFEST_ID,
        "cache_contract_id": CACHE_CONTRACT_ID,
        "pair_keys": list(PAIR_KEY),
        "pair_table_id": SEARCH_PAIR_TABLE_ID,
        "supported_search_splits": list(SUPPORTED_SEARCH_SPLITS),
        "disallowed_search_splits": list(DISALLOWED_SEARCH_SPLITS),
        "holdout_handling_rule": HOLDOUT_HANDLING_RULE,
        "pair_table_contract": {
            "row_granularity": "one_row_per_bacteria_phage_pair",
            "pair_join_keys": list(PAIR_KEY),
            "labels_read_only": True,
            "required_columns": [
                "pair_id",
                "bacteria",
                "phage",
                "split",
                "label_any_lysis",
                "training_weight_v3",
                "retained_for_autoresearch",
                "host_fasta_path",
                "phage_fasta_path",
            ],
        },
        "slot_order": [spec.slot_name for spec in SLOT_SPECS],
        "feature_slots": {
            spec.slot_name: {
                "entity_key": spec.entity_key,
                "join_keys": spec.join_keys,
                "column_family_prefix": spec.column_prefix,
                "block_role": spec.block_role,
                "reserved_feature_columns": resolved_slot_entries[spec.slot_name],
                "reserved_feature_column_count": len(resolved_slot_entries[spec.slot_name]),
            }
            for spec in SLOT_SPECS
        },
        "composability_contract": {
            "training_code_mutation_boundary": "train.py may consume the cache but must not rewrite schema or manifests",
            "cache_building_boundary": "prepare.py is the only supported path from raw inputs to the search cache",
            "split_visibility_rule": (
                "Only train and inner_val pair tables are exported into the search workspace; sealed holdout rows are "
                "kept outside the cache entirely."
            ),
            "feature_block_composition_rule": (
                "Each slot is joined independently to the pair tables by its declared key; slot files may add columns "
                "later but may not change slot names, join keys, or prefixes."
            ),
            "warm_cache_rule": (
                "Optional warm-cache artifacts are accelerators only. They must declare the same schema_manifest_id and "
                "match the fixed slot contract exactly."
            ),
        },
    }


def parse_warm_cache_manifest(path: Path) -> dict[str, Any]:
    manifest = read_json(path)
    if "slot_artifacts" not in manifest:
        raise ValueError(f"Warm-cache manifest is missing slot_artifacts: {path}")
    if not isinstance(manifest["slot_artifacts"], Mapping):
        raise ValueError(f"Warm-cache slot_artifacts must be a mapping: {path}")
    return manifest


def validate_warm_cache_manifest(path: Path, *, schema_manifest: Mapping[str, Any]) -> dict[str, Any]:
    manifest = parse_warm_cache_manifest(path)
    schema_manifest_id = schema_manifest["schema_manifest_id"]
    if manifest.get("schema_manifest_id") != schema_manifest_id:
        raise ValueError(
            "Warm-cache manifest schema mismatch: "
            f"expected {schema_manifest_id}, got {manifest.get('schema_manifest_id')}"
        )

    validated_slots: dict[str, Any] = {}
    for slot_name, descriptor in manifest["slot_artifacts"].items():
        if slot_name not in SLOT_SPEC_BY_NAME:
            raise ValueError(f"Warm-cache manifest declares unknown slot: {slot_name}")
        if not isinstance(descriptor, Mapping):
            raise ValueError(f"Warm-cache descriptor must be a mapping for slot {slot_name}")

        slot_spec = SLOT_SPEC_BY_NAME[slot_name]
        join_keys = [str(value) for value in descriptor.get("join_keys", [])]
        if join_keys != slot_spec.join_keys:
            raise ValueError(
                f"Warm-cache join keys do not match frozen contract for {slot_name}: "
                f"expected {slot_spec.join_keys}, got {join_keys}"
            )

        prefix = str(descriptor.get("column_family_prefix", ""))
        if prefix != slot_spec.column_prefix:
            raise ValueError(
                f"Warm-cache column prefix does not match frozen contract for {slot_name}: "
                f"expected {slot_spec.column_prefix}, got {prefix}"
            )

        columns = [str(value) for value in descriptor.get("columns", [])]
        if any(not column.startswith(slot_spec.column_prefix) for column in columns):
            raise ValueError(f"Warm-cache columns must stay inside slot prefix {slot_spec.column_prefix}: {slot_name}")

        artifact_path = path.parent / str(descriptor.get("path", ""))
        if not artifact_path.exists():
            raise FileNotFoundError(f"Warm-cache artifact not found for slot {slot_name}: {artifact_path}")

        header = load_csv_header(artifact_path)
        expected_header = slot_spec.join_keys + columns
        if header != expected_header:
            raise ValueError(
                f"Warm-cache artifact header mismatch for {slot_name}: expected {expected_header}, got {header}"
            )

        validated_slots[slot_name] = {
            "path": str(artifact_path),
            "join_keys": join_keys,
            "column_family_prefix": prefix,
            "columns": columns,
            "column_count": len(columns),
            "sha256": build_contract.sha256_file(artifact_path),
        }

    return {
        "warm_cache_manifest_path": str(path),
        "warm_cache_manifest_sha256": build_contract.sha256_file(path),
        "warm_cache_manifest_id": manifest.get("warm_cache_manifest_id", ""),
        "schema_manifest_id": manifest.get("schema_manifest_id"),
        "source_kind": manifest.get("source_kind", ""),
        "source_notes": manifest.get("source_notes", ""),
        "validated_slots": validated_slots,
    }


def load_csv_header(path: Path) -> list[str]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle)
        return next(reader)


def build_entity_path_map(
    *,
    selected_rows: Mapping[str, Sequence[Mapping[str, str]]],
    entity_key: str,
    path_key: str,
) -> dict[str, Path]:
    entity_paths: dict[str, Path] = {}
    for split_rows in selected_rows.values():
        for row in split_rows:
            if str(row["retained_for_autoresearch"]) != "1":
                continue
            entity_id = str(row[entity_key])
            path = Path(str(row[path_key]))
            existing = entity_paths.get(entity_id)
            if existing is not None and existing != path:
                raise ValueError(f"Entity {entity_id!r} resolved to multiple {path_key} values: {existing} vs {path}")
            entity_paths[entity_id] = path
    if not entity_paths:
        raise ValueError(f"No retained entities resolved for {entity_key}/{path_key}.")
    return entity_paths


def _load_tl17_frozen_runtime_assets(tl17_output_dir: Path) -> dict[str, Any]:
    manifest_path = tl17_output_dir / tl17_preprocessor.MANIFEST_FILENAME
    runtime_payload_path = tl17_output_dir / tl17_preprocessor.RUNTIME_FILENAME
    reference_fasta_path = tl17_output_dir / tl17_preprocessor.REFERENCE_FASTA_FILENAME
    reference_metadata_path = tl17_output_dir / tl17_preprocessor.REFERENCE_METADATA_FILENAME
    family_metadata_path = tl17_output_dir / tl17_preprocessor.FAMILY_METADATA_FILENAME
    schema_manifest_path = tl17_output_dir / tl17_runtime.SCHEMA_MANIFEST_FILENAME
    required_paths = (
        manifest_path,
        runtime_payload_path,
        reference_fasta_path,
        reference_metadata_path,
        family_metadata_path,
        schema_manifest_path,
    )
    for path in required_paths:
        if not path.exists():
            raise FileNotFoundError(f"Missing frozen TL17 runtime asset: {path}")

    manifest = read_json(manifest_path)
    runtime_payload = joblib.load(runtime_payload_path)
    if not isinstance(runtime_payload, Mapping):
        raise TypeError(f"TL17 runtime payload must be a mapping, got {type(runtime_payload)!r}")
    schema_manifest = read_json(schema_manifest_path)
    if str(schema_manifest.get("feature_block", "")) != tl17_runtime.TL17_BLOCK_ID:
        raise ValueError(
            f"TL17 schema manifest at {schema_manifest_path} declares unexpected feature block "
            f"{schema_manifest.get('feature_block')!r}"
        )

    reference_bank_provenance = {
        "task_id": str(manifest.get("task_id", "")),
        "tl17_output_dir": str(tl17_output_dir),
        "manifest_path": str(manifest_path),
        "manifest_sha256": build_contract.sha256_file(manifest_path),
        "runtime_payload_path": str(runtime_payload_path),
        "runtime_payload_sha256": build_contract.sha256_file(runtime_payload_path),
        "reference_fasta_path": str(reference_fasta_path),
        "reference_fasta_sha256": build_contract.sha256_file(reference_fasta_path),
        "reference_metadata_path": str(reference_metadata_path),
        "reference_metadata_sha256": build_contract.sha256_file(reference_metadata_path),
        "family_metadata_path": str(family_metadata_path),
        "family_metadata_sha256": build_contract.sha256_file(family_metadata_path),
        "schema_manifest_path": str(schema_manifest_path),
        "schema_manifest_sha256": build_contract.sha256_file(schema_manifest_path),
        "retained_family_count": int(manifest.get("counts", {}).get("retained_family_count", 0)),
        "retained_reference_protein_count": int(manifest.get("counts", {}).get("retained_reference_protein_count", 0)),
        "matching_policy": dict(manifest.get("matching_policy", {})),
        "projected_feature_csv_path_not_used": str(manifest.get("outputs", {}).get("projected_feature_csv", "")),
    }
    return {
        "manifest": manifest,
        "runtime_payload": dict(runtime_payload),
        "schema_manifest": schema_manifest,
        "reference_fasta_path": reference_fasta_path,
        "reference_bank_provenance": reference_bank_provenance,
    }


def namespace_slot_feature_rows(
    *,
    rows: Sequence[Mapping[str, object]],
    slot_spec: SlotSpec,
) -> tuple[list[str], list[dict[str, object]]]:
    namespaced_rows: list[dict[str, object]] = []
    feature_columns: list[str] = []
    for row in rows:
        namespaced_row: dict[str, object] = {slot_spec.entity_key: row[slot_spec.entity_key]}
        for key, value in row.items():
            if key == slot_spec.entity_key:
                continue
            column_name = f"{slot_spec.column_prefix}{key}"
            namespaced_row[column_name] = value
            if column_name not in feature_columns:
                feature_columns.append(column_name)
        namespaced_rows.append(namespaced_row)
    feature_columns.sort()
    return feature_columns, namespaced_rows


def try_reuse_slot(
    *,
    cache_dir: Path,
    slot_spec: SlotSpec,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
    path_key: str,
) -> Optional[dict[str, Any]]:
    """Validate an existing slot features.csv and return a summary dict if reusable, else None."""
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    index_path = slot_dir / ENTITY_INDEX_FILENAME

    if not feature_path.is_file():
        return None

    expected_entities = set(
        build_entity_path_map(
            selected_rows=split_rows,
            entity_key=slot_spec.entity_key,
            path_key=path_key,
        )
    )

    with feature_path.open("r", newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        headers = list(reader.fieldnames or [])
        found_entities: set[str] = set()
        row_count = 0
        for row in reader:
            found_entities.add(str(row[slot_spec.entity_key]))
            row_count += 1

    if found_entities != expected_entities:
        LOGGER.info(
            "Slot %s features.csv entity mismatch: expected %d, found %d — will rebuild",
            slot_spec.slot_name,
            len(expected_entities),
            len(found_entities),
        )
        return None

    feature_columns = [c for c in headers if c != slot_spec.entity_key]
    bad_columns = [c for c in feature_columns if not c.startswith(slot_spec.column_prefix)]
    if bad_columns:
        LOGGER.info(
            "Slot %s features.csv has columns outside prefix %s — will rebuild",
            slot_spec.slot_name,
            slot_spec.column_prefix,
        )
        return None

    # Rewrite the slot-level schema manifest so it reflects the actual columns.
    # Earlier partial runs may have written an empty manifest before features were materialized.
    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=row_count,
        reserved_feature_columns=feature_columns,
    )
    write_json(schema_path, schema_manifest)

    LOGGER.info(
        "Reusing existing %s slot: %d entities, %d feature columns",
        slot_spec.slot_name,
        row_count,
        len(feature_columns),
    )
    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(index_path),
        "schema_manifest_path": str(schema_path),
        "entity_count": row_count,
        "sha256": build_contract.sha256_file(index_path) if index_path.is_file() else "",
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
    }


def write_split_pair_tables(cache_dir: Path, split_rows: Mapping[str, Sequence[Mapping[str, str]]]) -> dict[str, Any]:
    pair_table_dir = cache_dir / "search_pairs"
    ensure_directory(pair_table_dir)
    pair_table_summaries: dict[str, Any] = {}

    for split_name, filename in (
        (build_contract.TRAIN_SPLIT, TRAIN_PAIR_TABLE_FILENAME),
        (build_contract.INNER_VAL_SPLIT, INNER_VAL_PAIR_TABLE_FILENAME),
    ):
        rows = [dict(row) for row in split_rows[split_name]]
        path = pair_table_dir / filename
        write_csv(path, fieldnames=list(rows[0].keys()), rows=rows)
        pair_table_summaries[split_name] = summarize_pair_table(path, rows)

    return pair_table_summaries


def summarize_pair_table(path: Path, rows: Sequence[Mapping[str, str]]) -> dict[str, Any]:
    retained_rows = [row for row in rows if str(row["retained_for_autoresearch"]) == "1"]
    return {
        "path": str(path),
        "sha256": build_contract.sha256_file(path),
        "row_count": len(rows),
        "retained_row_count": len(retained_rows),
        "bacteria_count": len({str(row["bacteria"]) for row in retained_rows}),
        "phage_count": len({str(row["phage"]) for row in retained_rows}),
        "label_counts": dict(Counter(str(row["label_any_lysis"]) for row in retained_rows)),
    }


def _retained_host_fastas(split_rows: Mapping[str, Sequence[Mapping[str, str]]]) -> dict[str, Path]:
    host_paths: dict[str, Path] = {}
    for rows in split_rows.values():
        for row in rows:
            if str(row["retained_for_autoresearch"]) != "1":
                continue
            bacteria = str(row["bacteria"])
            host_fasta_path = Path(str(row["host_fasta_path"]))
            if bacteria in host_paths and host_paths[bacteria] != host_fasta_path:
                raise ValueError(
                    f"Host {bacteria} resolved to multiple FASTA paths: {host_paths[bacteria]} vs {host_fasta_path}"
                )
            host_paths[bacteria] = host_fasta_path
    if not host_paths:
        raise ValueError("AUTORESEARCH retained set has zero host FASTAs for host-defense cache building.")
    return dict(sorted(host_paths.items()))


def _process_one_host_defense_cache_entry(
    assembly_path: Path,
    bacteria_id: str,
    output_dir: Path,
    panel_path: Path,
    models_dir: Path,
    force_run: bool,
    preserve_raw: bool,
) -> tuple[str, bool, str]:
    try:
        from lyzortx.pipeline.deployment_paired_features.derive_host_defense_features import (
            derive_host_defense_features,
        )

        derive_host_defense_features(
            assembly_path,
            bacteria_id=bacteria_id,
            output_dir=output_dir / bacteria_id,
            panel_defense_subtypes_path=panel_path,
            models_dir=models_dir,
            workers=1,
            force_model_update=False,
            model_install_mode=MODEL_INSTALL_MODE_FORBID,
            force_run=force_run,
            preserve_raw=preserve_raw,
        )
        return bacteria_id, True, "ok"
    except Exception as exc:
        return bacteria_id, False, str(exc)


def _build_host_defense_slot_artifact(
    *,
    args: argparse.Namespace,
    cache_dir: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, Any]:
    slot_dir = cache_dir / "feature_slots" / "host_defense"
    ensure_directory(slot_dir)
    per_host_output_dir = (
        args.host_defense_output_dir if args.host_defense_output_dir is not None else args.output_root / "host_defense"
    )
    host_fastas = _retained_host_fastas(split_rows)
    bacteria_ids = list(host_fastas.keys())
    coordinator_model_status = "not_requested"
    pending_host_count = 0
    skipped_host_count = 0
    failure_messages: list[dict[str, str]] = []
    cold_cache_elapsed_seconds = None

    if not args.host_defense_aggregate_only:
        LOGGER.info(
            "AR03 host defense: ensuring pinned release models once in %s before %d-worker fan-out",
            args.host_defense_models_dir,
            args.host_defense_max_workers,
        )
        coordinator_model_status = ensure_defense_finder_models(
            models_dir=args.host_defense_models_dir,
            force_update=args.host_defense_force_model_update,
        )
        pending: list[tuple[str, Path]] = []
        for bacteria_id, assembly_path in host_fastas.items():
            counts_path = per_host_output_dir / bacteria_id / PER_HOST_COUNTS_FILENAME
            if counts_path.exists() and not args.host_defense_force_run:
                skipped_host_count += 1
                continue
            pending.append((bacteria_id, assembly_path))

        pending_host_count = len(pending)
        if pending:
            start = datetime.now(timezone.utc)
            with ProcessPoolExecutor(max_workers=args.host_defense_max_workers) as pool:
                futures = {
                    pool.submit(
                        _process_one_host_defense_cache_entry,
                        assembly_path,
                        bacteria_id,
                        per_host_output_dir,
                        DEFAULT_PANEL_DEFENSE_SUBTYPES_PATH,
                        args.host_defense_models_dir,
                        args.host_defense_force_run,
                        args.host_defense_preserve_raw,
                    ): bacteria_id
                    for bacteria_id, assembly_path in pending
                }
                for future in as_completed(futures):
                    bacteria_id = futures[future]
                    try:
                        _, ok, message = future.result()
                    except Exception as exc:
                        ok, message = False, f"worker process error: {exc}"
                    if not ok:
                        failure_messages.append({"bacteria": bacteria_id, "message": message})
                        LOGGER.error("AR03 host defense failed for %s: %s", bacteria_id, message)
            cold_cache_elapsed_seconds = round((datetime.now(timezone.utc) - start).total_seconds(), 3)
            if failure_messages:
                raise RuntimeError(
                    "AR03 host-defense cache build failed for "
                    + ", ".join(entry["bacteria"] for entry in failure_messages)
                )
        else:
            LOGGER.info("AR03 host defense: all retained hosts already have per-host outputs; skipping worker fan-out")

    aggregate_path = slot_dir / SLOT_FEATURES_FILENAME
    aggregate_host_defense_csvs(
        per_host_output_dir,
        aggregate_path,
        DEFAULT_PANEL_DEFENSE_SUBTYPES_PATH,
        bacteria_ids=bacteria_ids,
    )

    raw_rows = load_csv_rows(aggregate_path)
    if len(raw_rows) != len(bacteria_ids):
        raise ValueError(
            f"AR03 host-defense aggregation row count mismatch: expected {len(bacteria_ids)}, got {len(raw_rows)}"
        )
    exported_columns = [
        f"{SLOT_SPEC_BY_NAME['host_defense'].column_prefix}{column}" for column in raw_rows[0] if column != "bacteria"
    ]
    projected_rows = []
    for row in raw_rows:
        projected_row = {"bacteria": str(row["bacteria"])}
        for column in raw_rows[0]:
            if column == "bacteria":
                continue
            projected_row[f"{SLOT_SPEC_BY_NAME['host_defense'].column_prefix}{column}"] = row[column]
        projected_rows.append(projected_row)
    write_csv(aggregate_path, fieldnames=["bacteria", *exported_columns], rows=projected_rows)

    build_manifest = {
        "task_id": "AR03",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": "host_defense",
        "per_host_output_dir": str(per_host_output_dir),
        "slot_artifact_path": str(aggregate_path),
        "retained_host_count": len(bacteria_ids),
        "retained_hosts": bacteria_ids,
        "runtime_controls": {
            "aggregate_only": args.host_defense_aggregate_only,
            "max_workers": args.host_defense_max_workers,
            "force_model_update": args.host_defense_force_model_update,
            "force_run": args.host_defense_force_run,
            "preserve_raw": args.host_defense_preserve_raw,
            "worker_model_install_mode": MODEL_INSTALL_MODE_FORBID,
        },
        "coordinator_model_status": coordinator_model_status,
        "cold_cache_elapsed_seconds": cold_cache_elapsed_seconds,
        "pending_host_count": pending_host_count,
        "skipped_host_count": skipped_host_count,
        "failure_messages": failure_messages,
        "guardrails": {
            "source_of_truth": "raw host FASTAs plus pinned Defense Finder release models",
            "forbidden_regressions": [
                "per-host model installation inside parallel workers",
                "source-checkout model directory confusion",
                "hidden model downloads during worker execution",
            ],
            "interpretation_limits": (
                "Detected defense hits are useful positive evidence. Zero counts are annotation-limited absences and "
                "must not be interpreted as clean biological absence."
            ),
        },
    }
    write_json(slot_dir / HOST_DEFENSE_BUILD_MANIFEST_FILENAME, build_manifest)
    return {
        "artifact_path": str(aggregate_path),
        "artifact_sha256": build_contract.sha256_file(aggregate_path),
        "columns": exported_columns,
        "column_count": len(exported_columns),
        "entity_count": len(projected_rows),
        "build_manifest_path": str(slot_dir / HOST_DEFENSE_BUILD_MANIFEST_FILENAME),
    }


def materialize_host_surface_slot(
    *,
    cache_dir: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
    max_workers: int,
    build_dir: Path,
) -> dict[str, Any]:
    slot_spec = SLOT_SPEC_BY_NAME["host_surface"]
    host_fasta_by_bacteria = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="host_fasta_path",
    )
    assemblies = [host_fasta_by_bacteria[bacteria] for bacteria in sorted(host_fasta_by_bacteria)]
    fast_path_result = run_all_host_surface.build_host_surface_rows_fast_path(
        assemblies=assemblies,
        output_dir=build_dir,
        max_workers=max(1, min(len(assemblies), max_workers)),
        include_lps_core_type=False,
    )
    if len(fast_path_result["rows"]) != len(assemblies):
        raise ValueError(
            "Host-surface fast path returned the wrong row count: "
            f"expected {len(assemblies)}, got {len(fast_path_result['rows'])}"
        )

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=fast_path_result["rows"], slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "source_runtime_id": fast_path_result["runtime_metadata"]["runtime_id"],
        "legacy_nhmmer_path_forbidden": fast_path_result["runtime_metadata"]["legacy_nhmmer_path_forbidden"],
        "source_schema_includes_lps_core_type": fast_path_result["schema"]["includes_lps_core_type"],
        "source_columns_dropped_for_autoresearch": ["host_lps_core_type"],
        "rebuildable_from_raw_fastas": True,
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_runtime": dict(fast_path_result["runtime_metadata"]),
        "build_dir": str(build_dir),
    }


def materialize_host_typing_slot(
    *,
    cache_dir: Path,
    output_root: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, Any]:
    slot_spec = SLOT_SPEC_BY_NAME["host_typing"]
    typing_schema = derive_host_typing_features.build_host_typing_schema()
    host_fasta_by_bacteria = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="host_fasta_path",
    )
    build_dir = output_root / HOST_TYPING_BUILD_DIRNAME
    rows: list[dict[str, object]] = []
    runtime_caveats: list[dict[str, str]] = []
    for bacteria in sorted(host_fasta_by_bacteria):
        assembly_path = host_fasta_by_bacteria[bacteria]
        result = derive_host_typing_features.derive_host_typing_features(
            assembly_path,
            bacteria_id=bacteria,
            output_dir=build_dir / bacteria,
            picard_metadata_path=None,
        )
        rows.append(dict(result["feature_row"]))
        runtime_caveats.extend(result["manifest"]["runtime_caveats"])
    if len(rows) != len(host_fasta_by_bacteria):
        raise ValueError(
            f"Host-typing materialization row count mismatch: expected {len(host_fasta_by_bacteria)}, got {len(rows)}"
        )

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=rows, slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "per_host_output_dir": str(build_dir),
        "caller_envs": typing_schema["caller_envs"],
        "panel_metadata_used_for_feature_construction": False,
        "rebuildable_from_raw_fastas": True,
        "runtime_caveat_count": len(runtime_caveats),
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    build_manifest = {
        "task_id": "AR05",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "per_host_output_dir": str(build_dir),
        "slot_artifact_path": str(feature_path),
        "retained_host_count": len(rows),
        "retained_hosts": sorted(host_fasta_by_bacteria),
        "caller_envs": typing_schema["caller_envs"],
        "guardrails": {
            "source_of_truth": "raw host FASTAs plus pinned phylogroup, serotype, and sequence-type caller envs",
            "panel_metadata_used_for_feature_construction": False,
            "comparison_paths_may_use_panel_metadata": True,
            "unresolved_calls_are_left_blank": True,
        },
        "runtime_caveats": runtime_caveats,
    }
    build_manifest_path = slot_dir / HOST_TYPING_BUILD_MANIFEST_FILENAME
    write_json(build_manifest_path, build_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_manifest_path": str(build_manifest_path),
        "build_dir": str(build_dir),
        "runtime_caveat_count": len(runtime_caveats),
    }


def materialize_host_stats_slot(
    *,
    cache_dir: Path,
    output_root: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, Any]:
    slot_spec = SLOT_SPEC_BY_NAME["host_stats"]
    host_fasta_by_bacteria = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="host_fasta_path",
    )
    build_dir = output_root / HOST_STATS_BUILD_DIRNAME
    n_hosts = len(host_fasta_by_bacteria)
    LOGGER.info("Host stats: computing sequence statistics for %d hosts", n_hosts)
    rows: list[dict[str, object]] = []
    for i, bacteria in enumerate(sorted(host_fasta_by_bacteria), 1):
        assembly_path = host_fasta_by_bacteria[bacteria]
        result = derive_host_stats_features.derive_host_stats_features(
            assembly_path,
            bacteria_id=bacteria,
            output_dir=build_dir / bacteria,
        )
        rows.append(dict(result["feature_row"]))
        if i % 50 == 0 or i == n_hosts:
            LOGGER.info("Host stats: %d/%d hosts completed", i, n_hosts)
    if len(rows) != len(host_fasta_by_bacteria):
        raise ValueError(
            f"Host-stats materialization row count mismatch: expected {len(host_fasta_by_bacteria)}, got {len(rows)}"
        )

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=rows, slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "per_host_output_dir": str(build_dir),
        "rebuildable_from_raw_fastas": True,
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    build_manifest = {
        "task_id": "AR05",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "per_host_output_dir": str(build_dir),
        "slot_artifact_path": str(feature_path),
        "retained_host_count": len(rows),
        "retained_hosts": sorted(host_fasta_by_bacteria),
        "guardrails": {
            "source_of_truth": "raw host FASTAs only",
            "panel_metadata_used": False,
            "low_cost_baseline_feature_family": True,
        },
        "exported_numeric_columns": feature_columns,
    }
    build_manifest_path = slot_dir / HOST_STATS_BUILD_MANIFEST_FILENAME
    write_json(build_manifest_path, build_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_manifest_path": str(build_manifest_path),
        "build_dir": str(build_dir),
    }


def materialize_phage_projection_slot(
    *,
    cache_dir: Path,
    output_root: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
    tl17_output_dir: Path,
) -> dict[str, Any]:
    slot_spec = SLOT_SPEC_BY_NAME["phage_projection"]
    phage_fasta_by_phage = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="phage_fasta_path",
    )
    retained_phages = sorted(phage_fasta_by_phage)
    tl17_assets = _load_tl17_frozen_runtime_assets(tl17_output_dir)
    build_dir = output_root / PHAGE_PROJECTION_BUILD_DIRNAME

    LOGGER.info(
        "Phage projection: starting batched mmseqs search for %d phages",
        len(retained_phages),
    )
    start = datetime.now(timezone.utc)
    rows = tl17_runtime.project_phage_feature_rows_batched(
        [phage_fasta_by_phage[phage] for phage in retained_phages],
        runtime_payload=tl17_assets["runtime_payload"],
        reference_fasta_path=tl17_assets["reference_fasta_path"],
        scratch_root=build_dir,
    )
    elapsed_seconds = round((datetime.now(timezone.utc) - start).total_seconds(), 3)
    LOGGER.info("Phage projection: batched mmseqs search finished in %.1fs", elapsed_seconds)
    if len(rows) != len(retained_phages):
        raise ValueError(
            f"Phage-projection materialization row count mismatch: expected {len(retained_phages)}, got {len(rows)}"
        )

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=rows, slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "build_dir": str(build_dir),
        "frozen_runtime_task_id": tl17_assets["reference_bank_provenance"]["task_id"],
        "batched_projection_path_used": True,
        "rebuildable_from_raw_fastas": True,
        "checked_in_projection_csv_used": False,
        "reference_bank_provenance": dict(tl17_assets["reference_bank_provenance"]),
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    build_manifest = {
        "task_id": "AR06",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "build_dir": str(build_dir),
        "slot_artifact_path": str(feature_path),
        "retained_phage_count": len(rows),
        "retained_phages": retained_phages,
        "projection_elapsed_seconds": elapsed_seconds,
        "build_runtime": {
            "projection_mode": "batched_mmseqs_easy_search",
            "frozen_runtime_task_id": tl17_assets["reference_bank_provenance"]["task_id"],
            "accepted_reference_hit_column": tl17_runtime.SUMMARY_HIT_COUNT_COLUMN,
        },
        "guardrails": {
            "source_of_truth": "raw phage FASTAs plus the frozen TL17 runtime payload and reference bank",
            "checked_in_projection_csv_used": False,
            "batched_projection_path_used": True,
            "per_phage_reference_rebuild_forbidden": True,
            "panel_only_host_metadata_used": False,
            "label_derived_pair_features_used": False,
        },
        "reference_bank_provenance": dict(tl17_assets["reference_bank_provenance"]),
    }
    build_manifest_path = slot_dir / PHAGE_PROJECTION_BUILD_MANIFEST_FILENAME
    write_json(build_manifest_path, build_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_manifest_path": str(build_manifest_path),
        "build_dir": str(build_dir),
        "build_runtime": dict(build_manifest["build_runtime"]),
        "reference_bank_provenance": dict(tl17_assets["reference_bank_provenance"]),
    }


def materialize_phage_stats_slot(
    *,
    cache_dir: Path,
    output_root: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, Any]:
    slot_spec = SLOT_SPEC_BY_NAME["phage_stats"]
    phage_fasta_by_phage = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="phage_fasta_path",
    )
    rows: list[dict[str, object]] = []
    n_phages = len(phage_fasta_by_phage)
    LOGGER.info("Phage stats: computing sequence statistics for %d phages", n_phages)
    for i, phage in enumerate(sorted(phage_fasta_by_phage), 1):
        fasta_path = phage_fasta_by_phage[phage]
        feature_row = derive_phage_stats_features.build_phage_stats_feature_row(
            fasta_path,
            phage_id=phage,
        )
        rows.append(dict(feature_row))
    LOGGER.info("Phage stats: %d/%d phages completed", n_phages, n_phages)
    if len(rows) != len(phage_fasta_by_phage):
        raise ValueError(
            f"Phage-stats materialization row count mismatch: expected {len(phage_fasta_by_phage)}, got {len(rows)}"
        )

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=rows, slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "rebuildable_from_raw_fastas": True,
        "low_cost_baseline_feature_family": True,
        "direct_feature_row_path_used": True,
        "per_phage_intermediate_artifacts_written": False,
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    build_manifest = {
        "task_id": "AR06",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "slot_artifact_path": str(feature_path),
        "retained_phage_count": len(rows),
        "retained_phages": sorted(phage_fasta_by_phage),
        "guardrails": {
            "source_of_truth": "raw phage FASTAs only",
            "panel_metadata_used": False,
            "low_cost_baseline_feature_family": True,
            "per_phage_intermediate_artifacts_written": False,
        },
        "exported_numeric_columns": feature_columns,
    }
    build_manifest_path = slot_dir / PHAGE_STATS_BUILD_MANIFEST_FILENAME
    write_json(build_manifest_path, build_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_manifest_path": str(build_manifest_path),
    }


def materialize_phage_kmer_slot(
    *,
    cache_dir: Path,
    output_root: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
) -> dict[str, Any]:
    from lyzortx.pipeline.track_d.steps.build_phage_genome_kmer_features import (
        compute_kmer_frequency_vector,
        read_fasta_records,
    )

    slot_spec = SLOT_SPEC_BY_NAME["phage_kmer"]
    phage_fasta_by_phage = build_entity_path_map(
        selected_rows=split_rows,
        entity_key=slot_spec.entity_key,
        path_key="phage_fasta_path",
    )
    k = PHAGE_KMER_K
    kmer_dim = PHAGE_KMER_DIM
    feature_names = [f"tetra_freq_{i:03d}" for i in range(kmer_dim)]

    rows: list[dict[str, object]] = []
    n_phages = len(phage_fasta_by_phage)
    LOGGER.info("Phage kmer: computing %d-mer frequency vectors for %d phages", k, n_phages)
    for i, phage in enumerate(sorted(phage_fasta_by_phage), 1):
        fasta_path = phage_fasta_by_phage[phage]
        records = read_fasta_records(fasta_path, protein=False)
        sequences = [record.sequence for record in records]

        combined_vector = np.zeros(kmer_dim, dtype=np.float64)
        total_windows = 0
        for sequence in sequences:
            vector = compute_kmer_frequency_vector(sequence, k=k)
            window_count = sum(
                1 for start in range(len(sequence) - k + 1) if all(c in "ACGT" for c in sequence[start : start + k])
            )
            if window_count:
                combined_vector += vector * window_count
            total_windows += window_count
        if total_windows == 0:
            raise ValueError(f"No valid {k}-mer windows found for phage {phage}")
        combined_vector /= total_windows

        row: dict[str, object] = {"phage": phage}
        for j, name in enumerate(feature_names):
            row[name] = round(float(combined_vector[j]), 6)
        rows.append(row)
        if i % 20 == 0 or i == n_phages:
            LOGGER.info("Phage kmer: %d/%d phages completed", i, n_phages)

    feature_columns, namespaced_rows = namespace_slot_feature_rows(rows=rows, slot_spec=slot_spec)
    slot_dir = cache_dir / "feature_slots" / slot_spec.slot_name
    feature_path = slot_dir / SLOT_FEATURES_FILENAME
    write_csv(feature_path, fieldnames=[slot_spec.entity_key, *feature_columns], rows=namespaced_rows)

    schema_manifest = build_slot_schema_manifest(
        slot_spec,
        row_count=len(namespaced_rows),
        reserved_feature_columns=feature_columns,
    )
    schema_manifest["materialization"] = {
        "feature_csv_path": str(feature_path),
        "rebuildable_from_raw_fastas": True,
        "kmer_k": k,
        "kmer_dim": kmer_dim,
        "svd_applied": False,
    }
    schema_path = slot_dir / SLOT_SCHEMA_FILENAME
    write_json(schema_path, schema_manifest)

    build_manifest = {
        "task_id": "AR06",
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "slot_name": slot_spec.slot_name,
        "slot_artifact_path": str(feature_path),
        "retained_phage_count": len(rows),
        "retained_phages": sorted(phage_fasta_by_phage),
        "guardrails": {
            "source_of_truth": "raw phage FASTAs only",
            "panel_metadata_used": False,
            "svd_applied": False,
            "kmer_k": k,
            "kmer_dim": kmer_dim,
        },
        "exported_numeric_columns": feature_columns,
    }
    build_manifest_path = slot_dir / PHAGE_KMER_BUILD_MANIFEST_FILENAME
    write_json(build_manifest_path, build_manifest)

    return {
        "entity_key": slot_spec.entity_key,
        "index_path": str(slot_dir / ENTITY_INDEX_FILENAME),
        "schema_manifest_path": str(schema_path),
        "entity_count": len(namespaced_rows),
        "sha256": build_contract.sha256_file(slot_dir / ENTITY_INDEX_FILENAME),
        "columns": feature_columns,
        "column_count": len(feature_columns),
        "feature_csv_path": str(feature_path),
        "feature_csv_sha256": build_contract.sha256_file(feature_path),
        "materialized_feature_columns": feature_columns,
        "materialized_feature_column_count": len(feature_columns),
        "build_manifest_path": str(build_manifest_path),
    }


def write_slot_indexes(
    cache_dir: Path,
    split_rows: Mapping[str, Sequence[Mapping[str, str]]],
    *,
    args: argparse.Namespace,
) -> dict[str, Any]:
    slot_root = cache_dir / "feature_slots"
    ensure_directory(slot_root)
    slot_summaries: dict[str, Any] = {}

    for slot_spec in SLOT_SPECS:
        slot_dir = slot_root / slot_spec.slot_name
        ensure_directory(slot_dir)
        rows = build_slot_index_rows(slot_spec=slot_spec, selected_rows=split_rows)
        index_path = slot_dir / ENTITY_INDEX_FILENAME
        write_csv(index_path, fieldnames=slot_spec.join_keys, rows=rows)

        slot_artifact_summary = {}
        if slot_spec.slot_name == "host_defense":
            slot_artifact_summary = _build_host_defense_slot_artifact(
                args=args,
                cache_dir=cache_dir,
                split_rows=split_rows,
            )

        schema_manifest = build_slot_schema_manifest(
            slot_spec,
            row_count=len(rows),
            reserved_feature_columns=slot_artifact_summary.get("columns", ()),
        )
        schema_path = slot_dir / SLOT_SCHEMA_FILENAME
        write_json(schema_path, schema_manifest)

        slot_summaries[slot_spec.slot_name] = {
            "entity_key": slot_spec.entity_key,
            "index_path": str(index_path),
            "schema_manifest_path": str(schema_path),
            "entity_count": len(rows),
            "sha256": build_contract.sha256_file(index_path),
            "columns": list(schema_manifest["reserved_feature_columns"]),
            "column_count": schema_manifest["reserved_feature_column_count"],
            **slot_artifact_summary,
        }

    return slot_summaries


def build_provenance_manifest(
    *,
    output_root: Path,
    cache_dir: Path,
    contract_manifest_path: Path,
    input_checksums_path: Path,
    pair_rows: Sequence[Mapping[str, str]],
    split_pair_tables: Mapping[str, Any],
    slot_summaries: Mapping[str, Any],
    warm_cache_validation: Optional[Mapping[str, Any]],
) -> dict[str, Any]:
    holdout_rows = [row for row in pair_rows if str(row["split"]) == build_contract.HOLDOUT_SPLIT]
    holdout_retained_rows = [row for row in holdout_rows if str(row["retained_for_autoresearch"]) == "1"]

    return {
        "task_id": TASK_ID,
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "schema_manifest_id": SCHEMA_MANIFEST_ID,
        "build_mode": "raw_inputs_plus_optional_warm_cache" if warm_cache_validation else "raw_inputs_only",
        "source_contract": {
            "pair_contract_manifest_path": str(contract_manifest_path),
            "pair_contract_manifest_sha256": build_contract.sha256_file(contract_manifest_path),
            "input_checksums_manifest_path": str(input_checksums_path),
            "input_checksums_manifest_sha256": build_contract.sha256_file(input_checksums_path),
            "output_root": str(output_root),
        },
        "search_workspace": {
            "cache_dir": str(cache_dir),
            "exported_splits": list(SUPPORTED_SEARCH_SPLITS),
            "disallowed_splits": list(DISALLOWED_SEARCH_SPLITS),
            "pair_tables": dict(split_pair_tables),
            "feature_slots": dict(slot_summaries),
        },
        "sealed_holdout": {
            "split_name": build_contract.HOLDOUT_SPLIT,
            "row_count": len(holdout_rows),
            "retained_row_count": len(holdout_retained_rows),
            "exported_to_search_cache": False,
            "rule": "Holdout labels and holdout-ready evaluation tables stay outside the RunPod workspace entirely.",
        },
        "warm_cache_validation": None if warm_cache_validation is None else dict(warm_cache_validation),
    }


def build_cache_manifest(
    *,
    cache_dir: Path,
    schema_manifest_path: Path,
    provenance_manifest_path: Path,
    split_pair_tables: Mapping[str, Any],
    slot_summaries: Mapping[str, Any],
) -> dict[str, Any]:
    return {
        "task_id": TASK_ID,
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "cache_contract_id": CACHE_CONTRACT_ID,
        "schema_manifest_id": SCHEMA_MANIFEST_ID,
        "cache_dir": str(cache_dir),
        "schema_manifest_path": str(schema_manifest_path),
        "provenance_manifest_path": str(provenance_manifest_path),
        "pair_tables": dict(split_pair_tables),
        "feature_slots": dict(slot_summaries),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    setup_logging()
    args = parse_args(argv)
    ensure_directory(args.output_root)
    ensure_directory(args.cache_dir)

    pair_table_path, contract_manifest_path, input_checksums_path = run_ar01_contract(args)
    pair_rows = load_csv_rows(pair_table_path)
    split_rows = select_search_rows(pair_rows)
    all_retained_rows = select_all_retained_rows(pair_rows)

    split_pair_tables = write_split_pair_tables(args.cache_dir, split_rows)
    slot_summaries = write_slot_indexes(args.cache_dir, all_retained_rows, args=args)

    host_surface_build_dir = args.host_surface_build_dir or (args.output_root / HOST_SURFACE_BUILD_DIRNAME)

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["host_surface"],
        split_rows=all_retained_rows,
        path_key="host_fasta_path",
    )
    if reused:
        slot_summaries["host_surface"] = reused
    else:
        LOGGER.info("Materializing host_surface slot")
        slot_summaries["host_surface"] = materialize_host_surface_slot(
            cache_dir=args.cache_dir,
            split_rows=all_retained_rows,
            max_workers=args.host_surface_max_workers,
            build_dir=host_surface_build_dir,
        )
        LOGGER.info("Completed host_surface slot")

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["host_typing"],
        split_rows=all_retained_rows,
        path_key="host_fasta_path",
    )
    if reused:
        slot_summaries["host_typing"] = reused
    else:
        LOGGER.info("Materializing host_typing slot")
        slot_summaries["host_typing"] = materialize_host_typing_slot(
            cache_dir=args.cache_dir,
            output_root=args.output_root,
            split_rows=all_retained_rows,
        )
        LOGGER.info("Completed host_typing slot")

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["host_stats"],
        split_rows=all_retained_rows,
        path_key="host_fasta_path",
    )
    if reused:
        slot_summaries["host_stats"] = reused
    else:
        LOGGER.info("Materializing host_stats slot")
        slot_summaries["host_stats"] = materialize_host_stats_slot(
            cache_dir=args.cache_dir,
            output_root=args.output_root,
            split_rows=all_retained_rows,
        )
        LOGGER.info("Completed host_stats slot")

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["phage_projection"],
        split_rows=all_retained_rows,
        path_key="phage_fasta_path",
    )
    if reused:
        slot_summaries["phage_projection"] = reused
    else:
        LOGGER.info("Materializing phage_projection slot")
        slot_summaries["phage_projection"] = materialize_phage_projection_slot(
            cache_dir=args.cache_dir,
            output_root=args.output_root,
            split_rows=all_retained_rows,
            tl17_output_dir=args.tl17_output_dir,
        )
        LOGGER.info("Completed phage_projection slot")

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["phage_stats"],
        split_rows=all_retained_rows,
        path_key="phage_fasta_path",
    )
    if reused:
        slot_summaries["phage_stats"] = reused
    else:
        LOGGER.info("Materializing phage_stats slot")
        slot_summaries["phage_stats"] = materialize_phage_stats_slot(
            cache_dir=args.cache_dir,
            output_root=args.output_root,
            split_rows=all_retained_rows,
        )
        LOGGER.info("Completed phage_stats slot")

    reused = try_reuse_slot(
        cache_dir=args.cache_dir,
        slot_spec=SLOT_SPEC_BY_NAME["phage_kmer"],
        split_rows=all_retained_rows,
        path_key="phage_fasta_path",
    )
    if reused:
        slot_summaries["phage_kmer"] = reused
    else:
        LOGGER.info("Materializing phage_kmer slot")
        slot_summaries["phage_kmer"] = materialize_phage_kmer_slot(
            cache_dir=args.cache_dir,
            output_root=args.output_root,
            split_rows=all_retained_rows,
        )
        LOGGER.info("Completed phage_kmer slot")

    schema_manifest = build_top_level_schema_manifest(
        slot_feature_columns={slot_name: summary["columns"] for slot_name, summary in slot_summaries.items()}
    )
    schema_manifest_path = args.cache_dir / SCHEMA_MANIFEST_FILENAME
    write_json(schema_manifest_path, schema_manifest)

    warm_cache_validation = None
    if args.warm_cache_manifest_path is not None:
        LOGGER.info("Validating optional warm-cache manifest: %s", args.warm_cache_manifest_path)
        warm_cache_validation = validate_warm_cache_manifest(
            args.warm_cache_manifest_path,
            schema_manifest=schema_manifest,
        )

    provenance_manifest = build_provenance_manifest(
        output_root=args.output_root,
        cache_dir=args.cache_dir,
        contract_manifest_path=contract_manifest_path,
        input_checksums_path=input_checksums_path,
        pair_rows=pair_rows,
        split_pair_tables=split_pair_tables,
        slot_summaries=slot_summaries,
        warm_cache_validation=warm_cache_validation,
    )
    provenance_manifest_path = args.cache_dir / PROVENANCE_MANIFEST_FILENAME
    write_json(provenance_manifest_path, provenance_manifest)

    cache_manifest = build_cache_manifest(
        cache_dir=args.cache_dir,
        schema_manifest_path=schema_manifest_path,
        provenance_manifest_path=provenance_manifest_path,
        split_pair_tables=split_pair_tables,
        slot_summaries=slot_summaries,
    )
    write_json(args.cache_dir / CACHE_MANIFEST_FILENAME, cache_manifest)

    LOGGER.info("AR02 completed: wrote search cache to %s", args.cache_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
