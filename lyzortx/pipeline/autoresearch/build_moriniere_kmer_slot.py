#!/usr/bin/env python3
"""SX12: Build a direct phage-side k-mer feature slot from Moriniere Dataset S6.

GenoPHI Dataset S6 (Moriniere 2026) lists 815 amino-acid k-mers (typically 5-mers) that are
predictive of phage receptor class at AUROC 0.99. GT06 used these as an intermediate to
predict a receptor class, then cross-termed the prediction with host OMP score — which
failed because of `omp-score-homogeneity`. This script emits the raw k-mer presence vector
as a direct phage-side feature slot so the k-mer signal can be fed straight to the model.

Inputs:
  - Moriniere Dataset S6 (column containing the k-mer sequences), default path
    `.scratch/genophi/Supplementary_Datasets_S1-S7.xlsx`
  - Guelin phage protein FASTAs under `lyzortx/generated_outputs/track_d/phage_protein_sets/protein_fastas/`
  - BASEL phage protein FASTAs under `.scratch/basel/pharokka_output/<phage>/phanotate.faa`

Output:
  - `.scratch/moriniere_kmer/features.csv` with schema
    `phage, phage_moriniere_kmer__<kmer1>, phage_moriniere_kmer__<kmer2>, ...`
    value = 1 if the k-mer occurs anywhere in the phage's predicted proteome, else 0.

Usage:
    python -m lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot
"""

from __future__ import annotations

import argparse
import csv
import logging
from pathlib import Path
from typing import Iterable, Sequence

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.predict_receptor_from_kmers import load_receptor_kmers

LOGGER = logging.getLogger(__name__)

DEFAULT_DATASET_PATH = Path(".scratch/genophi/Supplementary_Datasets_S1-S7.xlsx")
DEFAULT_GUELIN_PROTEIN_DIR = Path("lyzortx/generated_outputs/track_d/phage_protein_sets/protein_fastas")
DEFAULT_BASEL_PHAROKKA_DIR = Path(".scratch/basel/pharokka_output")
DEFAULT_OUTPUT_PATH = Path(".scratch/moriniere_kmer/features.csv")
FEATURE_PREFIX = "phage_moriniere_kmer__"


def read_fasta_sequences(path: Path) -> Iterable[str]:
    """Yield protein sequences from a FASTA file."""
    current: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current:
                    yield "".join(current)
                current = []
            else:
                current.append(line)
        if current:
            yield "".join(current)


def collect_guelin_phage_proteomes(protein_dir: Path) -> dict[str, str]:
    """Return phage -> concatenated proteome string (joined with '*' separator to avoid false cross-boundary kmers)."""
    proteomes: dict[str, str] = {}
    for fasta in sorted(protein_dir.glob("*.faa")):
        phage = fasta.stem
        proteomes[phage] = "*".join(read_fasta_sequences(fasta))
    LOGGER.info("Loaded %d Guelin phage proteomes from %s", len(proteomes), protein_dir)
    return proteomes


def collect_basel_phage_proteomes(pharokka_dir: Path) -> dict[str, str]:
    proteomes: dict[str, str] = {}
    for subdir in sorted(p for p in pharokka_dir.iterdir() if p.is_dir()):
        faa = subdir / "phanotate.faa"
        if not faa.exists():
            LOGGER.warning("Missing phanotate.faa for BASEL phage %s", subdir.name)
            continue
        proteomes[subdir.name] = "*".join(read_fasta_sequences(faa))
    LOGGER.info("Loaded %d BASEL phage proteomes from %s", len(proteomes), pharokka_dir)
    return proteomes


def enumerate_kmers(receptor_kmer_table: dict[str, set[str]]) -> list[str]:
    """Deduplicate k-mers across receptors; keep deterministic sorted order for reproducibility."""
    all_kmers: set[str] = set()
    for kmers in receptor_kmer_table.values():
        all_kmers.update(kmers)
    return sorted(all_kmers)


def build_feature_matrix(
    proteomes: dict[str, str],
    kmers: Sequence[str],
) -> list[dict[str, object]]:
    """For each phage, emit binary presence-absence per k-mer."""
    rows: list[dict[str, object]] = []
    for phage in sorted(proteomes):
        proteome = proteomes[phage]
        row: dict[str, object] = {"phage": phage}
        for kmer in kmers:
            row[f"{FEATURE_PREFIX}{kmer}"] = 1.0 if kmer in proteome else 0.0
        rows.append(row)
    return rows


def write_feature_csv(rows: list[dict[str, object]], kmers: Sequence[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["phage", *[f"{FEATURE_PREFIX}{k}" for k in kmers]]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-path", type=Path, default=DEFAULT_DATASET_PATH)
    parser.add_argument("--guelin-protein-dir", type=Path, default=DEFAULT_GUELIN_PROTEIN_DIR)
    parser.add_argument("--basel-pharokka-dir", type=Path, default=DEFAULT_BASEL_PHAROKKA_DIR)
    parser.add_argument("--output-path", type=Path, default=DEFAULT_OUTPUT_PATH)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)

    if not args.dataset_path.exists():
        raise FileNotFoundError(
            f"Moriniere Dataset S6 not found at {args.dataset_path}; re-download GenoPHI supplementary datasets"
        )

    receptor_kmers = load_receptor_kmers(args.dataset_path)
    kmers = enumerate_kmers(receptor_kmers)
    LOGGER.info("Deduplicated k-mers across receptors: %d unique k-mer strings", len(kmers))

    guelin = collect_guelin_phage_proteomes(args.guelin_protein_dir)
    basel = collect_basel_phage_proteomes(args.basel_pharokka_dir)
    proteomes = {**guelin, **basel}
    LOGGER.info("Combined proteomes: %d phages (%d Guelin + %d BASEL)", len(proteomes), len(guelin), len(basel))

    rows = build_feature_matrix(proteomes, kmers)

    # Light diagnostics: per-phage k-mer hit counts.
    hit_counts = [sum(1 for k in kmers if row[f"{FEATURE_PREFIX}{k}"]) for row in rows]
    LOGGER.info(
        "K-mer hit distribution across phages: min=%d, median=%d, max=%d",
        min(hit_counts),
        sorted(hit_counts)[len(hit_counts) // 2],
        max(hit_counts),
    )

    write_feature_csv(rows, kmers, args.output_path)
    LOGGER.info("Wrote %d phages × %d k-mer features to %s", len(rows), len(kmers), args.output_path)


if __name__ == "__main__":
    main()
