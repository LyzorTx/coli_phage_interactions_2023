#!/usr/bin/env python3
"""CH02: Emit ch02_fixed_folds.csv + diff report vs SPANDEX fold assignment.

Loads the SX10 clean-label training frame (same construction sx01_eval uses) and
reports, for each of the 369 training bacteria:
  - SPANDEX fold (hash on bacterium name) — buggy, allowed near-identical bacteria
    within a cv_group to split across folds
  - CHISEL fold (hash on cv_group, deterministic, same salt) — all bacteria in a
    cv_group share one fold
  - cv_group id

Also reports how many cv_groups were split across folds under the SPANDEX scheme
and how many bacteria changed fold under the CHISEL scheme.
"""

from __future__ import annotations

import hashlib
import logging
from collections import Counter, defaultdict
from pathlib import Path

import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    build_st03_training_frame,
    load_st03_holdout_frame,
)
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.sx01_eval import (
    FOLD_SALT,
    N_FOLDS,
    RAW_INTERACTIONS_PATH,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)

LOGGER = logging.getLogger(__name__)
OUTPUT_DIR = Path("lyzortx/generated_outputs/ch02_cv_group_fix")


def spandex_fold_for_bacterium(name: str) -> int:
    """Reproduce the SPANDEX-era buggy hashing for comparison only."""
    h = hashlib.sha256(f"{FOLD_SALT}:{name}".encode()).hexdigest()
    return int(h, 16) % N_FOLDS


def main() -> None:
    setup_logging()
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    training = build_st03_training_frame()
    holdout = load_st03_holdout_frame()
    full_frame = pd.concat([training, holdout], ignore_index=True)

    ambiguous = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean = full_frame[~full_frame["pair_id"].isin(ambiguous)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean), clean["bacteria"].nunique())

    mapping = bacteria_to_cv_group_map(clean)
    chisel_folds = assign_bacteria_folds(mapping)

    rows = []
    for bac, cv_group in sorted(mapping.items()):
        rows.append(
            {
                "bacteria": bac,
                "cv_group": cv_group,
                "fold_spandex": spandex_fold_for_bacterium(bac),
                "fold_chisel": chisel_folds[bac],
            }
        )
    df = pd.DataFrame(rows).sort_values(["cv_group", "bacteria"])
    df.to_csv(OUTPUT_DIR / "ch02_fixed_folds.csv", index=False)
    LOGGER.info("Wrote %s (%d bacteria)", OUTPUT_DIR / "ch02_fixed_folds.csv", len(df))

    split_cv_groups = df.groupby("cv_group")["fold_spandex"].nunique().reset_index(name="n_folds_spandex")
    n_split_spandex = (split_cv_groups["n_folds_spandex"] > 1).sum()
    LOGGER.info(
        "SPANDEX (bacterium-hashed) folds split %d / %d cv_groups across multiple folds",
        int(n_split_spandex),
        len(split_cv_groups),
    )

    chisel_split = df.groupby("cv_group")["fold_chisel"].nunique()
    assert (chisel_split == 1).all(), "CHISEL fold assignment leaked a cv_group across folds"

    changed = (df["fold_spandex"] != df["fold_chisel"]).sum()
    LOGGER.info("Bacteria whose fold changed: %d / %d", int(changed), len(df))

    fold_counts_spandex = Counter(df["fold_spandex"])
    fold_counts_chisel = Counter(df["fold_chisel"])
    LOGGER.info("Fold sizes (SPANDEX hashing): %s", dict(sorted(fold_counts_spandex.items())))
    LOGGER.info("Fold sizes (CHISEL hashing):  %s", dict(sorted(fold_counts_chisel.items())))

    multi_bac = defaultdict(list)
    for bac, cv_group in mapping.items():
        multi_bac[cv_group].append(bac)
    multi_bac = {g: b for g, b in multi_bac.items() if len(b) > 1}
    LOGGER.info("cv_groups with >1 bacterium: %d / %d", len(multi_bac), len(mapping) and len(set(mapping.values())))


if __name__ == "__main__":
    main()
