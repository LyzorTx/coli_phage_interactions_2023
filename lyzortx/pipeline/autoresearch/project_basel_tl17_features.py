#!/usr/bin/env python3
"""SX06: Project BASEL phage RBP proteomes through the deployable TL17 runtime.

SX02 zero-filled the `phage_projection` slot for BASEL phages despite the TL17
reference bank already existing, which compromised SX03 Arm B/C. This script calls the
existing TL17 runtime against BASEL FASTAs and merges the results into the extended
`phage_projection` slot CSV at `.scratch/basel/feature_slots/phage_projection/features.csv`.

Reuses:
  - `project_phage_feature_rows_batched` (lyzortx/pipeline/track_l/steps/deployable_tl17_runtime.py)
  - tl17_rbp_runtime.joblib + tl17_rbp_reference_bank.faa
    (lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/)

Usage:
    python -m lyzortx.pipeline.autoresearch.project_basel_tl17_features
"""

from __future__ import annotations

import argparse
import logging
from pathlib import Path
from typing import Sequence

import joblib
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.track_l.steps.deployable_tl17_runtime import (
    FAMILY_COLUMN_PREFIX,
    SUMMARY_HIT_COUNT_COLUMN,
    project_phage_feature_rows_batched,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_BASEL_GENOMES_DIR = Path(".scratch/basel/genomes")
DEFAULT_RUNTIME_PAYLOAD_PATH = Path(
    "lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/tl17_rbp_runtime.joblib"
)
DEFAULT_REFERENCE_FASTA_PATH = Path(
    "lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/tl17_rbp_reference_bank.faa"
)
DEFAULT_SLOT_FEATURES_PATH = Path(".scratch/basel/feature_slots/phage_projection/features.csv")
DEFAULT_SCRATCH_ROOT = Path(".scratch/basel/tl17_projection")

SLOT_HIT_COUNT_COLUMN = f"phage_projection__{SUMMARY_HIT_COUNT_COLUMN}"


def list_basel_phages(genomes_dir: Path) -> list[Path]:
    paths = sorted(p for p in genomes_dir.glob("*.fna") if p.stem.startswith("Bas"))
    if not paths:
        raise FileNotFoundError(f"No BASEL FASTAs under {genomes_dir} matching Bas*.fna")
    return paths


def merge_basel_features_into_slot(
    slot_features_path: Path,
    basel_feature_rows: Sequence[dict[str, object]],
) -> pd.DataFrame:
    """Overwrite BASEL rows in the existing slot CSV with projected feature values."""
    slot_df = pd.read_csv(slot_features_path)
    slot_columns = [c for c in slot_df.columns if c != "phage"]

    # Map runtime column names → slot column names.
    rename_map: dict[str, str] = {}
    for runtime_col in next(iter(basel_feature_rows), {}):
        if runtime_col == "phage":
            continue
        if runtime_col == SUMMARY_HIT_COUNT_COLUMN:
            slot_col = SLOT_HIT_COUNT_COLUMN
        else:
            slot_col = f"phage_projection__{runtime_col}"
        if slot_col not in slot_columns:
            raise KeyError(
                f"Runtime column {runtime_col!r} has no matching slot column {slot_col!r} in {slot_features_path}"
            )
        rename_map[runtime_col] = slot_col

    basel_df = pd.DataFrame(basel_feature_rows).rename(columns=rename_map)
    basel_df = basel_df[["phage", *[rename_map[c] for c in rename_map if c != "phage"]]]
    missing_slot_cols = set(slot_columns) - set(basel_df.columns)
    for col in missing_slot_cols:
        basel_df[col] = 0.0
    basel_df = basel_df[["phage", *slot_columns]]

    slot_df = slot_df.set_index("phage")
    basel_df = basel_df.set_index("phage")

    overlap = slot_df.index.intersection(basel_df.index)
    if len(overlap) != len(basel_df):
        missing = set(basel_df.index) - set(slot_df.index)
        raise KeyError(f"BASEL phages not present in slot CSV: {sorted(missing)}")
    slot_df.loc[overlap, slot_columns] = basel_df.loc[overlap, slot_columns].values
    return slot_df.reset_index()


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--basel-genomes-dir", type=Path, default=DEFAULT_BASEL_GENOMES_DIR)
    parser.add_argument("--runtime-payload-path", type=Path, default=DEFAULT_RUNTIME_PAYLOAD_PATH)
    parser.add_argument("--reference-fasta-path", type=Path, default=DEFAULT_REFERENCE_FASTA_PATH)
    parser.add_argument("--slot-features-path", type=Path, default=DEFAULT_SLOT_FEATURES_PATH)
    parser.add_argument("--scratch-root", type=Path, default=DEFAULT_SCRATCH_ROOT)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)

    phage_paths = list_basel_phages(args.basel_genomes_dir)
    LOGGER.info("Projecting %d BASEL phages through the deployable TL17 runtime", len(phage_paths))

    LOGGER.info("Loading TL17 runtime payload from %s", args.runtime_payload_path)
    runtime_payload = joblib.load(args.runtime_payload_path)

    feature_rows = project_phage_feature_rows_batched(
        phage_paths=phage_paths,
        runtime_payload=runtime_payload,
        reference_fasta_path=args.reference_fasta_path,
        scratch_root=args.scratch_root,
    )
    nonzero_count = sum(
        1 for row in feature_rows if any(v and v != 0 for k, v in row.items() if k.startswith(FAMILY_COLUMN_PREFIX))
    )
    LOGGER.info(
        "Projected %d BASEL phages; %d have at least one non-zero TL17 family score",
        len(feature_rows),
        nonzero_count,
    )

    merged = merge_basel_features_into_slot(args.slot_features_path, feature_rows)
    merged.to_csv(args.slot_features_path, index=False)

    # Sanity: count BASEL nonzero rows after merge.
    num_cols = [c for c in merged.columns if c != "phage"]
    basel_mask = merged["phage"].str.startswith("Bas")
    basel_nonzero = int((merged.loc[basel_mask, num_cols].abs().sum(axis=1) > 0).sum())
    LOGGER.info("Post-merge: BASEL rows with nonzero TL17 features = %d / %d", basel_nonzero, int(basel_mask.sum()))


if __name__ == "__main__":
    main()
