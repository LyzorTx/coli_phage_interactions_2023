"""Shared AUTORESEARCH runtime contract helpers and constants."""

from __future__ import annotations

import hashlib
from dataclasses import dataclass
from pathlib import Path
from typing import Iterable

DEFAULT_OUTPUT_ROOT = Path("lyzortx/generated_outputs/autoresearch")
DEFAULT_CACHE_DIR = DEFAULT_OUTPUT_ROOT / "search_cache_v1"

CACHE_MANIFEST_FILENAME = "ar02_search_cache_manifest_v1.json"
SCHEMA_MANIFEST_FILENAME = "ar02_schema_manifest_v1.json"
PROVENANCE_MANIFEST_FILENAME = "ar02_provenance_manifest_v1.json"
TRAIN_PAIR_TABLE_FILENAME = "train_pairs.csv"
INNER_VAL_PAIR_TABLE_FILENAME = "inner_val_pairs.csv"
ENTITY_INDEX_FILENAME = "entity_index.csv"
SLOT_SCHEMA_FILENAME = "schema_manifest.json"
SLOT_FEATURE_TABLE_FILENAME = "feature_table.csv"
SLOT_FEATURES_FILENAME = "features.csv"

CACHE_CONTRACT_ID = "autoresearch_search_cache_v1"
SCHEMA_MANIFEST_ID = "autoresearch_feature_schema_v1"
SEARCH_PAIR_TABLE_ID = "autoresearch_search_pair_tables_v1"
HOLDOUT_HANDLING_RULE = "sealed_holdout_outside_workspace"

TRAIN_SPLIT = "train"
INNER_VAL_SPLIT = "inner_val"
HOLDOUT_SPLIT = "holdout"
SPLIT_ORDER = (TRAIN_SPLIT, INNER_VAL_SPLIT, HOLDOUT_SPLIT)
SUPPORTED_SEARCH_SPLITS = (TRAIN_SPLIT, INNER_VAL_SPLIT)
DISALLOWED_SEARCH_SPLITS = (HOLDOUT_SPLIT,)
PAIR_KEY = ("pair_id", "bacteria", "phage")


@dataclass(frozen=True)
class SlotSpec:
    slot_name: str
    entity_key: str
    column_prefix: str
    block_role: str
    description: str

    @property
    def join_keys(self) -> list[str]:
        return [self.entity_key]


SLOT_SPECS = (
    SlotSpec(
        slot_name="host_defense",
        entity_key="bacteria",
        column_prefix="host_defense__",
        block_role="host",
        description="Reserved host defense-system features derived from raw host assemblies.",
    ),
    SlotSpec(
        slot_name="host_surface",
        entity_key="bacteria",
        column_prefix="host_surface__",
        block_role="host",
        description="Reserved host surface and adsorption-related features derived from raw host assemblies.",
    ),
    SlotSpec(
        slot_name="host_typing",
        entity_key="bacteria",
        column_prefix="host_typing__",
        block_role="host",
        description="Reserved host typing calls derived from raw host assemblies.",
    ),
    SlotSpec(
        slot_name="host_stats",
        entity_key="bacteria",
        column_prefix="host_stats__",
        block_role="host",
        description="Reserved low-cost host sequence statistics derived from raw host assemblies.",
    ),
    SlotSpec(
        slot_name="phage_projection",
        entity_key="phage",
        column_prefix="phage_projection__",
        block_role="phage",
        description="Reserved phage projection features derived from raw phage genomes.",
    ),
    SlotSpec(
        slot_name="phage_stats",
        entity_key="phage",
        column_prefix="phage_stats__",
        block_role="phage",
        description="Reserved low-cost phage sequence statistics derived from raw phage genomes.",
    ),
    SlotSpec(
        slot_name="phage_kmer",
        entity_key="phage",
        column_prefix="phage_kmer__",
        block_role="phage",
        description="Raw tetranucleotide (k=4) frequency vectors from phage genomes. No SVD reduction.",
    ),
)
SLOT_SPEC_BY_NAME = {spec.slot_name: spec for spec in SLOT_SPECS}


def sha256_file(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def normalized_hash_01(text: str) -> float:
    digest = hashlib.sha256(text.encode("utf-8")).digest()
    value = int.from_bytes(digest[:8], byteorder="big", signed=False)
    return value / float(2**64 - 1)


def sha256_strings(values: Iterable[str]) -> str:
    digest = hashlib.sha256()
    for value in values:
        digest.update(value.encode("utf-8"))
        digest.update(b"\n")
    return digest.hexdigest()
