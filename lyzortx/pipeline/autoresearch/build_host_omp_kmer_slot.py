#!/usr/bin/env python3
"""SX13: Build host-side OMP 5-mer feature slot.

For each of 369 clinical E. coli hosts, extract the best phmmer-matched protein per core OMP
(BTUB, FADL, FHUA, LAMB, LPTD, NFRA, OMPA, OMPC, OMPF, TOLC, TSX, PQQU). Enumerate amino-acid
5-mers per (host, OMP), build a global per-OMP dictionary, prune to k-mers with 5-95%
across-host frequency, and emit a binary presence-absence slot CSV.

This inverts the Moriniere trick onto the host side: whole-gene HMM scores are CV 0.01-0.17 across
clinical strains (omp-score-homogeneity), but loop-level allelic variation is preserved in the
coding sequence and discoverable as 5-mer presence-absence.

Inputs:
  - Per-host predicted proteomes at lyzortx/generated_outputs/autoresearch/host_surface_cache_build/{host}/predicted_proteins.faa
  - OMP reference at lyzortx/pipeline/track_l/assets/tl15_omp_reference_proteins.faa

Outputs:
  - .scratch/host_omp_proteins/{host}/{omp}.faa (per-host best-hit OMP protein)
  - .scratch/host_omp_kmer/features.csv (presence-absence slot)

Usage:
    python -m lyzortx.pipeline.autoresearch.build_host_omp_kmer_slot
"""

from __future__ import annotations

import argparse
import csv
import logging
import re
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

from lyzortx.log_config import setup_logging

LOGGER = logging.getLogger(__name__)

DEFAULT_HOST_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/host_surface_cache_build")
DEFAULT_OMP_REFERENCE = Path("lyzortx/pipeline/track_l/assets/tl15_omp_reference_proteins.faa")
DEFAULT_PROTEIN_OUTPUT_DIR = Path(".scratch/host_omp_proteins")
DEFAULT_SLOT_OUTPUT_PATH = Path(".scratch/host_omp_kmer/features.csv")
FEATURE_PREFIX = "host_omp_kmer__"
KMER_LEN = 5
MIN_FREQ = 0.05
MAX_FREQ = 0.95


def parse_fasta(path: Path) -> dict[str, str]:
    """Return name -> sequence dict from a FASTA file. Names are the first whitespace-token of header."""
    seqs: dict[str, str] = {}
    name: str | None = None
    parts: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if name is not None:
                    seqs[name] = "".join(parts)
                name = line[1:].split()[0]
                parts = []
            else:
                parts.append(line)
        if name is not None:
            seqs[name] = "".join(parts)
    return seqs


def load_omp_reference_names(path: Path) -> list[str]:
    """Return the OMP gene names from the reference FASTA in source order.

    Accepts headers like '>sp|P06129|BTUB_ECOLI ...'. Returns ['BTUB', 'FADL', ...].
    """
    names: list[str] = []
    with path.open("r", encoding="utf-8") as handle:
        for line in handle:
            if not line.startswith(">"):
                continue
            header = line[1:].strip()
            # Pattern: sp|ACCESSION|GENE_ECOLI ...
            m = re.match(r"sp\|[^|]+\|([A-Z0-9]+)_ECOLI", header)
            if m:
                names.append(m.group(1))
            else:
                # Fall back to first token
                names.append(header.split()[0])
    return names


def extract_best_hits_for_host(
    *,
    host: str,
    proteome_path: Path,
    omp_reference_path: Path,
    output_dir: Path,
    omp_names: Sequence[str],
) -> dict[str, str]:
    """Run phmmer and write best-hit-per-OMP FASTAs to output_dir/{host}/{omp}.faa.

    Returns dict mapping OMP name to the matched protein sequence (empty if no hit).
    """
    host_out_dir = output_dir / host
    host_out_dir.mkdir(parents=True, exist_ok=True)

    # Skip if all per-OMP outputs already exist (idempotent re-runs).
    expected_files = [host_out_dir / f"{omp}.faa" for omp in omp_names]
    if all(p.exists() for p in expected_files):
        return {omp: _read_first_seq(p) for omp, p in zip(omp_names, expected_files)}

    proteome_seqs = parse_fasta(proteome_path)

    # phmmer: queries = OMP reference, db = host proteome
    with tempfile.NamedTemporaryFile("w", suffix=".tbl", delete=False) as tbl:
        tbl_path = Path(tbl.name)
    try:
        subprocess.run(
            [
                "phmmer",
                "--noali",
                "--cpu",
                "1",
                "--tblout",
                str(tbl_path),
                str(omp_reference_path),
                str(proteome_path),
            ],
            check=True,
            capture_output=True,
        )

        # Parse tblout: cols are target_name, target_acc, query_name, query_acc, evalue, score, ...
        # Pick best score per query (each query is one OMP reference).
        best_per_query: dict[str, tuple[str, float]] = {}
        with tbl_path.open("r", encoding="utf-8") as f:
            for line in f:
                if line.startswith("#") or not line.strip():
                    continue
                parts = line.split()
                target = parts[0]
                query = parts[2]
                score = float(parts[5])
                if query not in best_per_query or score > best_per_query[query][1]:
                    best_per_query[query] = (target, score)

        result: dict[str, str] = {}
        for omp in omp_names:
            # Find query header matching this OMP name
            matching_query = next((q for q in best_per_query if re.match(rf"sp\|[^|]+\|{omp}_ECOLI", q)), None)
            if matching_query is None:
                # No hit — write empty placeholder
                (host_out_dir / f"{omp}.faa").write_text("", encoding="utf-8")
                result[omp] = ""
                continue
            target, _ = best_per_query[matching_query]
            seq = proteome_seqs.get(target, "")
            (host_out_dir / f"{omp}.faa").write_text(f">{host}__{omp}__{target}\n{seq}\n", encoding="utf-8")
            result[omp] = seq

        return result
    finally:
        tbl_path.unlink(missing_ok=True)


def _read_first_seq(path: Path) -> str:
    """Return the first sequence in a single-record FASTA, or empty string if file is empty."""
    if not path.exists() or path.stat().st_size == 0:
        return ""
    seqs = parse_fasta(path)
    return next(iter(seqs.values()), "")


def enumerate_kmers(seq: str, k: int) -> set[str]:
    """Return all unique k-mers in seq."""
    if len(seq) < k:
        return set()
    return {seq[i : i + k] for i in range(len(seq) - k + 1)}


def discover_hosts(host_cache_dir: Path) -> list[str]:
    """Return sorted host names that have a predicted_proteins.faa file."""
    hosts: list[str] = []
    for d in sorted(host_cache_dir.iterdir()):
        if d.is_dir() and (d / "predicted_proteins.faa").exists():
            hosts.append(d.name)
    return hosts


def build_kmer_matrix(
    host_omp_seqs: dict[str, dict[str, str]],
    omp_names: Sequence[str],
    *,
    k: int = KMER_LEN,
    min_freq: float = MIN_FREQ,
    max_freq: float = MAX_FREQ,
) -> tuple[list[str], list[dict[str, object]]]:
    """Build the pruned k-mer presence-absence matrix.

    Per OMP, enumerate the union of host k-mers, count host coverage, keep only k-mers with
    coverage in [min_freq, max_freq] of hosts. Returns (column_order, rows).
    """
    hosts = sorted(host_omp_seqs.keys())
    n_hosts = len(hosts)
    columns: list[str] = []
    # per-OMP retained k-mers
    retained_per_omp: dict[str, list[str]] = {}

    for omp in omp_names:
        host_kmers: dict[str, set[str]] = {}
        all_kmers: set[str] = set()
        for h in hosts:
            seq = host_omp_seqs[h].get(omp, "")
            kmers = enumerate_kmers(seq, k)
            host_kmers[h] = kmers
            all_kmers.update(kmers)

        # Count coverage per k-mer
        coverage: dict[str, int] = {km: 0 for km in all_kmers}
        for kmers in host_kmers.values():
            for km in kmers:
                coverage[km] += 1

        retained = [km for km in sorted(all_kmers) if min_freq * n_hosts <= coverage[km] <= max_freq * n_hosts]
        retained_per_omp[omp] = retained
        columns.extend(f"{FEATURE_PREFIX}{omp}__{km}" for km in retained)
        LOGGER.info(
            "  %s: %d unique k-mers, %d retained (%.0f%%-%.0f%% frequency band)",
            omp,
            len(all_kmers),
            len(retained),
            min_freq * 100,
            max_freq * 100,
        )

    LOGGER.info("Total retained k-mer features across %d OMPs: %d", len(omp_names), len(columns))

    rows: list[dict[str, object]] = []
    for h in hosts:
        row: dict[str, object] = {"bacteria": h}
        for omp in omp_names:
            seq = host_omp_seqs[h].get(omp, "")
            kmers = enumerate_kmers(seq, k)
            for km in retained_per_omp[omp]:
                row[f"{FEATURE_PREFIX}{omp}__{km}"] = 1.0 if km in kmers else 0.0
        rows.append(row)
    return columns, rows


def write_slot_csv(rows: list[dict[str, object]], columns: list[str], output_path: Path) -> None:
    output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["bacteria", *columns]
    with output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow(row)


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--host-cache-dir", type=Path, default=DEFAULT_HOST_CACHE_DIR)
    parser.add_argument("--omp-reference", type=Path, default=DEFAULT_OMP_REFERENCE)
    parser.add_argument("--protein-output-dir", type=Path, default=DEFAULT_PROTEIN_OUTPUT_DIR)
    parser.add_argument("--slot-output-path", type=Path, default=DEFAULT_SLOT_OUTPUT_PATH)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)

    if not args.omp_reference.exists():
        raise FileNotFoundError(f"OMP reference FASTA missing: {args.omp_reference}")
    if not args.host_cache_dir.exists():
        raise FileNotFoundError(f"Host cache dir missing: {args.host_cache_dir}")

    omp_names = load_omp_reference_names(args.omp_reference)
    LOGGER.info("Core OMPs (%d): %s", len(omp_names), ", ".join(omp_names))

    hosts = discover_hosts(args.host_cache_dir)
    LOGGER.info("Discovered %d hosts with predicted_proteins.faa", len(hosts))

    # Step 1: per-host phmmer extraction
    host_omp_seqs: dict[str, dict[str, str]] = {}
    for i, host in enumerate(hosts, 1):
        proteome_path = args.host_cache_dir / host / "predicted_proteins.faa"
        host_omp_seqs[host] = extract_best_hits_for_host(
            host=host,
            proteome_path=proteome_path,
            omp_reference_path=args.omp_reference,
            output_dir=args.protein_output_dir,
            omp_names=omp_names,
        )
        if i % 25 == 0 or i == len(hosts):
            covered = sum(1 for s in host_omp_seqs[host].values() if s)
            LOGGER.info("phmmer %d/%d hosts (last: %s, %d/%d OMPs hit)", i, len(hosts), host, covered, len(omp_names))

    # Step 2: build pruned k-mer matrix
    LOGGER.info(
        "Building k-mer presence-absence matrix (k=%d, freq band %.0f%%-%.0f%%)",
        KMER_LEN,
        MIN_FREQ * 100,
        MAX_FREQ * 100,
    )
    columns, rows = build_kmer_matrix(host_omp_seqs, omp_names)

    write_slot_csv(rows, columns, args.slot_output_path)
    LOGGER.info("Wrote %d hosts × %d k-mer features to %s", len(rows), len(columns), args.slot_output_path)


if __name__ == "__main__":
    main()
