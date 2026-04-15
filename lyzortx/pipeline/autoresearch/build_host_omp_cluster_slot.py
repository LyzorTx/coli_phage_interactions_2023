#!/usr/bin/env python3
"""SX13 Path 1 fallback: per-OMP MMseqs2 cluster IDs as host-side categorical features.

For each of 12 core OMPs, gather all per-host best-hit protein sequences (from
build_host_omp_kmer_slot's extraction step), MMseqs2-cluster them at 99% identity, and emit each
host's cluster ID per OMP as a categorical feature.

Inputs:
  - .scratch/host_omp_proteins/{host}/{omp}.faa (built by build_host_omp_kmer_slot)

Outputs:
  - .scratch/host_omp_cluster/features.csv with columns:
    bacteria, host_omp_cluster__BTUB, host_omp_cluster__OMPC, ...

Usage:
    python -m lyzortx.pipeline.autoresearch.build_host_omp_cluster_slot
"""

from __future__ import annotations

import argparse
import csv
import logging
import shutil
import subprocess
import tempfile
from pathlib import Path
from typing import Sequence

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.build_host_omp_kmer_slot import (
    DEFAULT_OMP_REFERENCE,
    DEFAULT_PROTEIN_OUTPUT_DIR,
    load_omp_reference_names,
    parse_fasta,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_SLOT_OUTPUT_PATH = Path(".scratch/host_omp_cluster/features.csv")
FEATURE_PREFIX = "host_omp_cluster__"
CLUSTER_IDENTITY = 0.99
NO_HIT_CLUSTER_ID = "NO_HIT"


def cluster_omp_sequences(
    omp: str,
    host_seqs: dict[str, str],
    *,
    identity: float = CLUSTER_IDENTITY,
) -> dict[str, str]:
    """MMseqs2-cluster host sequences for one OMP. Returns host -> cluster_id mapping.

    Hosts with empty sequence (no phmmer hit) get NO_HIT_CLUSTER_ID.
    """
    hosts_with_seq = {h: s for h, s in host_seqs.items() if s}
    no_hit_hosts = [h for h, s in host_seqs.items() if not s]
    LOGGER.info(
        "%s: clustering %d host sequences (%d hosts had no phmmer hit)",
        omp,
        len(hosts_with_seq),
        len(no_hit_hosts),
    )

    if not hosts_with_seq:
        return {h: NO_HIT_CLUSTER_ID for h in host_seqs}

    with tempfile.TemporaryDirectory(prefix=f"mmseqs_{omp}_") as tmp:
        tmp_path = Path(tmp)
        input_fasta = tmp_path / "input.faa"
        with input_fasta.open("w", encoding="utf-8") as f:
            for h, s in hosts_with_seq.items():
                f.write(f">{h}\n{s}\n")

        # MMseqs2 easy-cluster
        result_prefix = tmp_path / "cluster_result"
        scratch_dir = tmp_path / "scratch"
        scratch_dir.mkdir()
        subprocess.run(
            [
                "mmseqs",
                "easy-cluster",
                str(input_fasta),
                str(result_prefix),
                str(scratch_dir),
                "--min-seq-id",
                f"{identity}",
                "-c",
                "0.8",
                "--cov-mode",
                "0",
                "-v",
                "1",
            ],
            check=True,
            capture_output=True,
        )

        # Output: <prefix>_cluster.tsv with columns rep_id, member_id (one row per member)
        cluster_tsv = Path(f"{result_prefix}_cluster.tsv")
        if not cluster_tsv.exists():
            raise FileNotFoundError(f"MMseqs2 output missing: {cluster_tsv}")

        host_cluster: dict[str, str] = {}
        # Use cluster representative as the cluster ID, prefixed for readability
        with cluster_tsv.open("r", encoding="utf-8") as f:
            for line in f:
                parts = line.rstrip("\n").split("\t")
                if len(parts) != 2:
                    continue
                rep, member = parts
                host_cluster[member] = f"{omp}_C_{rep}"

        # Add no-hit hosts
        for h in no_hit_hosts:
            host_cluster[h] = NO_HIT_CLUSTER_ID

        n_clusters = len(set(v for v in host_cluster.values() if v != NO_HIT_CLUSTER_ID))
        LOGGER.info("  %s: %d clusters across %d hosts", omp, n_clusters, len(host_cluster))
        return host_cluster


def parse_args(argv: Sequence[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--omp-reference", type=Path, default=DEFAULT_OMP_REFERENCE)
    parser.add_argument("--protein-input-dir", type=Path, default=DEFAULT_PROTEIN_OUTPUT_DIR)
    parser.add_argument("--slot-output-path", type=Path, default=DEFAULT_SLOT_OUTPUT_PATH)
    return parser.parse_args(argv)


def main(argv: Sequence[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)

    if not shutil.which("mmseqs"):
        raise RuntimeError("mmseqs CLI not on PATH; activate the phage_env environment first")
    if not args.protein_input_dir.exists():
        raise FileNotFoundError(
            f"Per-host OMP proteins missing at {args.protein_input_dir}; run build_host_omp_kmer_slot first"
        )

    omp_names = load_omp_reference_names(args.omp_reference)
    hosts = sorted(d.name for d in args.protein_input_dir.iterdir() if d.is_dir())
    LOGGER.info("Discovered %d hosts under %s", len(hosts), args.protein_input_dir)
    LOGGER.info("Clustering %d OMPs at %.0f%% identity", len(omp_names), CLUSTER_IDENTITY * 100)

    # Gather per-OMP host sequences
    omp_host_clusters: dict[str, dict[str, str]] = {}
    for omp in omp_names:
        host_seqs: dict[str, str] = {}
        for h in hosts:
            faa = args.protein_input_dir / h / f"{omp}.faa"
            if not faa.exists() or faa.stat().st_size == 0:
                host_seqs[h] = ""
                continue
            seqs = parse_fasta(faa)
            host_seqs[h] = next(iter(seqs.values()), "")
        omp_host_clusters[omp] = cluster_omp_sequences(omp, host_seqs)

    # Write slot CSV
    args.slot_output_path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = ["bacteria", *[f"{FEATURE_PREFIX}{omp}" for omp in omp_names]]
    with args.slot_output_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for h in hosts:
            row: dict[str, object] = {"bacteria": h}
            for omp in omp_names:
                row[f"{FEATURE_PREFIX}{omp}"] = omp_host_clusters[omp].get(h, NO_HIT_CLUSTER_ID)
            writer.writerow(row)

    LOGGER.info(
        "Wrote %d hosts × %d categorical OMP cluster columns to %s", len(hosts), len(omp_names), args.slot_output_path
    )


if __name__ == "__main__":
    main()
