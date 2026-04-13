#!/usr/bin/env python3
"""SX05: Align data/interactions/interaction_matrix.csv with the paper's own MLC protocol.

The paper's committed interaction_matrix.csv encodes MLC 0-4 where MLC=4 is a morphological
distinction ("entire lysis of the bacterial lawn at 5 x 10^6") derived from plaque images.
Our binary 0/1/n spot-test data cannot reconstruct that distinction, and our pipeline's old
DILUTION_WEIGHT_MAP repurposed MLC=4 to mean "lysis at 5 x 10^4" — contradicting the paper's
own protocol, which explicitly excludes the unreplicated 5 x 10^4 pfu/ml observation.

What this script does:
  - cap every MLC=4 cell at MLC=3 (we cannot distinguish individualized plaques from entire
    lawn lysis from binary spot data; MLC=3 and MLC=4 collapse into MLC=3)
  - leave MLC ∈ {0, 1, 2, 3} cells untouched (they encode the paper's strict replicated
    judgment and our raw data cannot reproduce it without access to the paper's plaque imagery)
  - keep the matrix dimensions so downstream consumers keep their bacteria/phage universe

What this script does NOT do:
  - re-derive MLC from raw_interactions.csv: the paper's MLC protocol requires morphology
    (MLC=3 vs MLC=4) and likely replicated-positive thresholds that our binary raw data cannot
    reproduce. See track_SPANDEX.md for the discrepancy audit.

Paper quotes (Gaborieau 2024 Methods, "Evaluating phage-bacteria interaction outcomes by
plaque assay experiments"):

    "The outcome of interaction at 5 x 10^4 pfu/ml was not taken into account in the
    calculation of the MLC score because it was not verified by a replicate."

    "Interactions observed at the lowest phage concentration (5 x 10^6 pfu/ml) are
    distinguished between 3 (individualized lysis plaque) and 4 (entire lysis of the
    bacterial lawn)."

Usage:
    python -m lyzortx.pipeline.autoresearch.regenerate_interaction_matrix
"""

from __future__ import annotations

import argparse
import csv
import logging
from collections import Counter, defaultdict
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Sequence

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.track_a.steps.build_track_a_foundation import (
    DILUTION_WEIGHT_MAP,
    EXCLUDED_LOG_DILUTIONS,
    find_best_dilution_any_lysis,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")
DEFAULT_INTERACTION_MATRIX_PATH = Path("data/interactions/interaction_matrix.csv")
MAX_MLC = max(DILUTION_WEIGHT_MAP.values())


def read_raw_rows(path: Path) -> List[Dict[str, str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=";")
        return [{k: (v.strip() if isinstance(v, str) else "") for k, v in row.items()} for row in reader]


def build_dilution_counts_by_pair(
    raw_rows: Iterable[Mapping[str, str]],
) -> Dict[tuple, Dict[int, Counter]]:
    """Group raw rows by (bacteria, phage) and count score occurrences per non-excluded dilution."""
    pair_counts: Dict[tuple, Dict[int, Counter]] = defaultdict(lambda: defaultdict(Counter))
    for row in raw_rows:
        dilution = int(row["log_dilution"])
        if dilution in EXCLUDED_LOG_DILUTIONS:
            continue
        pair_counts[(row["bacteria"], row["phage"])][dilution][row["score"]] += 1
    return pair_counts


def compute_mlc_for_pair(dilution_counts: Mapping[int, Counter]) -> int:
    """Return MLC {0, 1, 2, 3} for one pair given per-dilution score counters.

    Uses `find_best_dilution_any_lysis` (lysis = any score='1' observation at the given
    dilution). This is a permissive rule that differs from the paper's protocol, so callers
    that need paper-consistent MLC should use the paper's committed matrix instead.
    """
    best_dilution = find_best_dilution_any_lysis(dilution_counts)
    if best_dilution is None:
        return 0
    return DILUTION_WEIGHT_MAP[best_dilution]


def cap_matrix_cells(cells: Mapping[tuple, str], cap: int) -> Dict[tuple, str]:
    """Cap numeric MLC cells at `cap`; preserve non-numeric (blank / missing) cells.

    Preserves the original formatting (e.g., "1.0" stays "1.0") for cells that do not need
    capping, so downstream diff output reflects only the MLC=4 → MLC=cap change.
    """
    capped: Dict[tuple, str] = {}
    cap_str = f"{float(cap):.1f}"  # match "x.0" formatting of the paper's matrix
    for key, value in cells.items():
        if value == "":
            capped[key] = value
            continue
        try:
            parsed = int(float(value))
        except ValueError:
            capped[key] = value
            continue
        if parsed > cap:
            capped[key] = cap_str
        else:
            capped[key] = value
    return capped


def read_existing_matrix(matrix_path: Path) -> tuple[List[str], List[str], Dict[tuple, str]]:
    with matrix_path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.reader(handle, delimiter=";")
        header = next(reader)
        phages = header[1:]
        bacteria: List[str] = []
        cells: Dict[tuple, str] = {}
        for row in reader:
            b = row[0]
            bacteria.append(b)
            for phage, value in zip(phages, row[1:]):
                cells[(b, phage)] = value
    return bacteria, phages, cells


def write_matrix(
    matrix_path: Path,
    bacteria: Sequence[str],
    phages: Sequence[str],
    cells: Mapping[tuple, str],
    line_terminator: str = "\r\n",
) -> Dict[str, int]:
    """Write wide-format matrix; return string-value distribution for logging.

    Defaults to CRLF to match the committed paper matrix, keeping the git diff focused on the
    MLC=4 → MLC=3 value change instead of a cosmetic line-ending flip.
    """
    distribution: Counter = Counter()
    matrix_path.parent.mkdir(parents=True, exist_ok=True)
    with matrix_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.writer(handle, delimiter=";", lineterminator=line_terminator)
        writer.writerow(["bacteria", *phages])
        for b in bacteria:
            row: List[str] = [b]
            for p in phages:
                value = cells.get((b, p), "")
                distribution[value] += 1
                row.append(value)
            writer.writerow(row)
    return dict(distribution)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--raw-interactions-path", type=Path, default=DEFAULT_RAW_INTERACTIONS_PATH)
    parser.add_argument("--matrix-path", type=Path, default=DEFAULT_INTERACTION_MATRIX_PATH)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)

    if not args.matrix_path.exists():
        raise FileNotFoundError(
            f"Expected committed interaction_matrix.csv at {args.matrix_path}; SX05 caps existing cells rather than "
            "reconstructing from raw, because raw binary data cannot reproduce the paper's morphological MLC split."
        )

    bacteria, phages, original_cells = read_existing_matrix(args.matrix_path)
    before_counts = Counter(original_cells.values())
    LOGGER.info("Existing matrix: %d bacteria x %d phages", len(bacteria), len(phages))
    for value, count in sorted(before_counts.items()):
        LOGGER.info("  before MLC=%s: %d cells", value or "<blank>", count)

    capped_cells = cap_matrix_cells(original_cells, MAX_MLC)
    after_counts = Counter(capped_cells.values())

    LOGGER.info("Writing capped MLC matrix (MLC=4 → MLC=3) to %s", args.matrix_path)
    write_matrix(args.matrix_path, bacteria, phages, capped_cells)
    for value, count in sorted(after_counts.items()):
        LOGGER.info("  after MLC=%s: %d cells", value or "<blank>", count)


if __name__ == "__main__":
    main()
