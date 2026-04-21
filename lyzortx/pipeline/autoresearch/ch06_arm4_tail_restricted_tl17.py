#!/usr/bin/env python3
"""CH06 Arm 4: tail-protein-restricted TL17 BLAST.

Narrows the TL17 BLAST query input to Pharokka-annotated tail / baseplate /
RBP proteins before projection against the RBP reference bank. Shang 2025
observed that tail proteins alone match full-genome phage-host prediction
performance — Arm 4 tests whether restricting the query side (not the reference
side, which is already RBP-only) carries discrimination signal missed by the
current "all query proteins against RBP bank" search.

Plan.yml caveats (these remain honest constraints):
- The existing reference bank is already `tl17_rbp_reference_bank.faa` (built
  via `classify_rbp_genes` over the same 96 Guelin phages). Query restriction
  therefore operates on the phage's *own* proteome. If non-tail proteins rarely
  cross-hit the RBP bank (which the narrow, RBP-focused reference makes likely),
  the tail restriction is a near-no-op — which is itself a reportable outcome.
- BASEL pre-flight: all 52 BASEL phages carry tail/baseplate CDS annotations
  (verified at category='tail' ∪ RBP_PATTERNS on the merged Pharokka TSV), so
  the arm is not skipped on that side.

Pipeline:

1. `extract_tail_proteins`: for each of 148 phages, parse the Pharokka merged
   TSV; keep CDS records whose `category ∈ {"tail", "baseplate"}` or whose
   `annot` matches any of the RBP/tail regex patterns from `parse_annotations`.
   For Guelin, look up protein sequences from the track_d per-phage FASTAs by
   coordinate overlap (start, end, contig match the `start=/end=/contig=`
   metadata in the FASTA header). For BASEL, look up by CDS name against the
   per-phage `phanotate.faa` since its headers use the same CDS id as the
   merged TSV gene field.

2. `write_tail_restricted_fasta`: concatenate all extracted tail proteins into
   a single combined query FASTA with `<phage>|query_prot_<idx:04d>` headers
   (matching the TL17 runtime schema so downstream projection reuses the
   batched mmseqs path without code changes).

3. `project_tail_restricted`: call MMseqs2 `easy-search` against the existing
   `tl17_rbp_reference_bank.faa`, aggregate per-family hits exactly as
   `project_phage_feature_rows_batched` does, and materialize a
   phage_projection CSV with the same 33-column TL17 schema. Phages that have
   zero tail-annotated proteins end up with a zero row (documented honest
   signal — Arm 4 cannot rescue an un-annotated phage).

4. `run_arm4_training_eval`: CH05 two-axis retraining with the tail-restricted
   phage_projection slot patched in via `patch_context_with_extended_slots`.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
import sys
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.track_l.steps.parse_annotations import (
    RBP_PATTERNS,
    matches_any_pattern,
    parse_merged_tsv,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_GUELIN_PHAROKKA_DIR = Path("data/annotations/pharokka")
DEFAULT_BASEL_PHAROKKA_DIR = Path(".scratch/basel/pharokka_output")
DEFAULT_GUELIN_PROTEIN_DIR = Path("lyzortx/generated_outputs/track_d/phage_protein_sets/protein_fastas")
DEFAULT_REFERENCE_FASTA = Path(
    "lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/tl17_rbp_reference_bank.faa"
)
DEFAULT_ORIGINAL_SLOT_CSV = Path(".scratch/basel/feature_slots/phage_projection/features.csv")
DEFAULT_SLOT_CSV = Path(".scratch/basel/feature_slots_arm4/phage_projection/features.csv")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch06_arm4_tail_restricted_tl17")
DEFAULT_SCRATCH_DIR = Path(".scratch/ch06_arm4")

TAIL_CATEGORIES = {"tail", "baseplate"}
MMSEQS_MIN_QUERY_COVERAGE = 0.70
MMSEQS_MAX_TARGET_SEQS = 20


def _is_tail_or_rbp(annot: str, category: str) -> bool:
    """Category 'tail' / 'baseplate' OR annotation matches RBP/tail regex patterns."""
    if category in TAIL_CATEGORIES:
        return True
    return matches_any_pattern(annot or "", RBP_PATTERNS)


def _parse_header_metadata(header: str) -> dict[str, str]:
    """Parse ' contig=... start=... end=... strand=...' metadata from a FASTA header."""
    meta: dict[str, str] = {}
    for token in header.split()[1:]:
        if "=" in token:
            k, v = token.split("=", 1)
            meta[k] = v
    return meta


def _read_protein_fasta_with_coords(path: Path) -> list[tuple[dict[str, str], str]]:
    """Read a protein FASTA, returning [(metadata_dict, sequence), ...] per protein."""
    entries: list[tuple[dict[str, str], str]] = []
    current_meta: dict[str, str] | None = None
    current_seq: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_meta is not None:
                    entries.append((current_meta, "".join(current_seq)))
                current_meta = _parse_header_metadata(line)
                current_meta["protein_id"] = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_meta is not None:
            entries.append((current_meta, "".join(current_seq)))
    return entries


def _read_named_fasta(path: Path) -> dict[str, str]:
    """Read FASTA as {protein_id: sequence} — protein_id = first whitespace-separated token."""
    result: dict[str, str] = {}
    current_id: str | None = None
    current_seq: list[str] = []
    with path.open("r", encoding="utf-8") as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if current_id is not None:
                    result[current_id] = "".join(current_seq)
                current_id = line[1:].split()[0]
                current_seq = []
            else:
                current_seq.append(line)
        if current_id is not None:
            result[current_id] = "".join(current_seq)
    return result


def extract_guelin_tail_proteins(
    phage: str,
    merged_tsv: Path,
    protein_fasta: Path,
) -> list[tuple[str, str]]:
    """Extract tail/RBP protein sequences from a Guelin phage via coordinate match.

    Returns [(protein_id, sequence), ...] in deterministic order.
    """
    cds_records = parse_merged_tsv(merged_tsv)
    tail_records = [r for r in cds_records if _is_tail_or_rbp(r.annot, r.category)]
    if not tail_records:
        return []

    proteins = _read_protein_fasta_with_coords(protein_fasta)
    # Index proteins by (contig, start, end) for O(1) lookup.
    coord_index: dict[tuple[str, int, int], str] = {}
    for meta, seq in proteins:
        contig = meta.get("contig", phage)
        try:
            start = int(meta["start"])
            end = int(meta["end"])
        except KeyError:
            continue
        coord_index[(contig, start, end)] = seq

    extracted: list[tuple[str, str]] = []
    for record in tail_records:
        key = (record.contig, record.start, record.stop)
        if key in coord_index:
            extracted.append((record.gene, coord_index[key]))
    return extracted


def extract_basel_tail_proteins(
    phage: str,
    pharokka_dir: Path,
) -> list[tuple[str, str]]:
    """Extract tail/RBP proteins from a BASEL phage via CDS name match."""
    merged_tsv = pharokka_dir / f"{phage}_cds_final_merged_output.tsv"
    phanotate_faa = pharokka_dir / "phanotate.faa"
    if not merged_tsv.exists() or not phanotate_faa.exists():
        return []
    cds_records = parse_merged_tsv(merged_tsv)
    tail_records = [r for r in cds_records if _is_tail_or_rbp(r.annot, r.category)]
    if not tail_records:
        return []
    name_to_seq = _read_named_fasta(phanotate_faa)
    return [(r.gene, name_to_seq[r.gene]) for r in tail_records if r.gene in name_to_seq]


def build_tail_restricted_query_fasta(
    guelin_pharokka_dir: Path,
    guelin_protein_dir: Path,
    basel_pharokka_dir: Path,
    guelin_panel: set[str],
    output_fasta: Path,
) -> dict[str, int]:
    """Write a combined FASTA of tail/RBP proteins for all 148 phages.

    Header format: `>{phage}|query_prot_{index:04d}` — matching the TL17 runtime's
    expected schema so downstream aggregation by phage id works without changes.
    """
    output_fasta.parent.mkdir(parents=True, exist_ok=True)
    per_phage_count: dict[str, int] = {}
    with output_fasta.open("w", encoding="utf-8") as out:
        # Guelin (96 panel phages).
        for phage in sorted(guelin_panel):
            merged = guelin_pharokka_dir / f"{phage}_cds_final_merged_output.tsv"
            protein = guelin_protein_dir / f"{phage}.faa"
            if not merged.exists() or not protein.exists():
                LOGGER.warning("Guelin %s missing annotation or protein FASTA; zero tail proteins", phage)
                per_phage_count[phage] = 0
                continue
            tail_prots = extract_guelin_tail_proteins(phage, merged, protein)
            for idx, (_gene, seq) in enumerate(tail_prots, 1):
                out.write(f">{phage}|query_prot_{idx:04d}\n{seq}\n")
            per_phage_count[phage] = len(tail_prots)

        # BASEL (52 phages).
        for basel_dir in sorted(basel_pharokka_dir.iterdir()):
            if not basel_dir.is_dir() or not basel_dir.name.startswith("Bas"):
                continue
            phage = basel_dir.name
            tail_prots = extract_basel_tail_proteins(phage, basel_dir)
            for idx, (_gene, seq) in enumerate(tail_prots, 1):
                out.write(f">{phage}|query_prot_{idx:04d}\n{seq}\n")
            per_phage_count[phage] = len(tail_prots)

    LOGGER.info(
        "Tail-restricted query FASTA: %d phages, %d total tail proteins, %d zero-tail phages",
        len(per_phage_count),
        sum(per_phage_count.values()),
        sum(1 for v in per_phage_count.values() if v == 0),
    )
    return per_phage_count


def run_mmseqs_search(
    query_fasta: Path,
    reference_fasta: Path,
    output_tsv: Path,
    scratch_dir: Path,
) -> None:
    """Run MMseqs2 easy-search with the same column layout as the track_l runtime."""
    tmp_dir = scratch_dir / "mmseqs_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)
    cmd = [
        "mmseqs",
        "easy-search",
        str(query_fasta),
        str(reference_fasta),
        str(output_tsv),
        str(tmp_dir),
        "--max-seqs",
        str(MMSEQS_MAX_TARGET_SEQS),
        "--format-output",
        "query,target,pident,alnlen,qstart,qend,qlen,tstart,tend,tlen,evalue,bits",
    ]
    LOGGER.info("Running tail-restricted MMseqs2: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)


def aggregate_tl17_features(
    hits_tsv: Path,
    runtime_payload_path: Path,
    all_phages: Sequence[str],
    min_query_coverage: float = MMSEQS_MIN_QUERY_COVERAGE,
) -> pd.DataFrame:
    """Aggregate protein-level mmseqs hits to phage × family feature rows.

    Uses the TL17 runtime payload to resolve reference_id → family_id, then
    per-family takes the max percent identity across passing hits (same as
    `project_phage_feature_rows_batched`).
    """
    import joblib

    from lyzortx.pipeline.track_l.steps.deployable_tl17_runtime import (
        SUMMARY_HIT_COUNT_COLUMN,
        parse_runtime_payload,
    )

    runtime_payload = joblib.load(runtime_payload_path)
    family_rows, reference_rows, _policy = parse_runtime_payload(runtime_payload)
    reference_to_family = {row.reference_id: row.family_id for row in reference_rows}
    family_schema = {row.family_id: row.column_name for row in family_rows}

    hits = pd.read_csv(
        hits_tsv,
        sep="\t",
        header=None,
        names=[
            "query",
            "target",
            "pident",
            "alnlen",
            "qstart",
            "qend",
            "qlen",
            "tstart",
            "tend",
            "tlen",
            "evalue",
            "bits",
        ],
    )
    hits["phage"] = hits["query"].str.split("|", n=1).str[0]
    hits["query_coverage"] = (hits["qend"] - hits["qstart"] + 1) / hits["qlen"]
    hits = hits[hits["query_coverage"] >= min_query_coverage]
    hits["family_id"] = hits["target"].map(reference_to_family)
    hits = hits.dropna(subset=["family_id"])

    # Per (phage, family): max pident.
    per_phage_family = hits.groupby(["phage", "family_id"])["pident"].max().unstack(fill_value=0.0)
    per_phage_hits = hits.groupby("phage").size()

    rows: list[dict[str, object]] = []
    for phage in all_phages:
        row: dict[str, object] = {"phage": phage}
        for family_id, column_name in family_schema.items():
            score = 0.0
            if phage in per_phage_family.index and family_id in per_phage_family.columns:
                score = float(per_phage_family.loc[phage, family_id])
            row[column_name] = score
        row[SUMMARY_HIT_COUNT_COLUMN] = int(per_phage_hits.get(phage, 0))
        rows.append(row)
    return pd.DataFrame(rows)


def materialize_arm4_slot(
    feature_df: pd.DataFrame,
    original_slot_csv: Path,
    slot_csv: Path,
) -> pd.DataFrame:
    """Convert runtime feature rows to phage_projection slot CSV schema.

    The runtime emits column names like `tl17_phage_rbp_family_105_present`.
    The slot CSV prepends `phage_projection__` to each numeric column.
    """
    original = pd.read_csv(original_slot_csv)
    slot_columns = [c for c in original.columns if c != "phage"]

    rename_map: dict[str, str] = {}
    for src_col in feature_df.columns:
        if src_col == "phage":
            continue
        target = f"phage_projection__{src_col}"
        if target not in slot_columns:
            raise KeyError(
                f"Runtime column {src_col!r} has no matching slot column {target!r}. "
                f"Arm 4 cannot materialize a schema mismatch."
            )
        rename_map[src_col] = target

    renamed = feature_df.rename(columns=rename_map)
    ordered = renamed[["phage", *slot_columns]]
    slot_csv.parent.mkdir(parents=True, exist_ok=True)
    ordered.to_csv(slot_csv, index=False)
    LOGGER.info("Arm 4 slot written to %s (%d phages, %d features)", slot_csv, len(ordered), len(slot_columns))
    return ordered


def run_variance_preflight(
    slot_df: pd.DataFrame,
    original_slot_csv: Path,
    interactions_path: Path,
    output_json: Path,
) -> dict[str, object]:
    """Three-subset CV + Cohen's d gate, matching Arms 2/3."""
    original = pd.read_csv(original_slot_csv)
    feature_cols = [c for c in original.columns if c != "phage"]
    original_nonzero = original[feature_cols].abs().sum(axis=1) > 0
    nonzero_phages = set(original.loc[original_nonzero, "phage"])

    guelin = [p for p in slot_df["phage"] if not p.startswith("Bas")]
    basel_all = [p for p in slot_df["phage"] if p.startswith("Bas")]
    basel_nonzero_proj = [p for p in basel_all if p in nonzero_phages]
    basel_zero_proj = [p for p in basel_all if p not in nonzero_phages]

    interactions = pd.read_csv(interactions_path, sep=";", low_memory=False)
    interactions["score"] = interactions["score"].astype(str)
    guelin_labels = (interactions["score"] == "1").astype(int)
    guelin_lysis = (
        pd.DataFrame({"phage": interactions["phage"], "label": guelin_labels}).groupby("phage")["label"].mean()
    )
    basel_path = Path(".scratch/genophi_data/BASEL_ECOR_interaction_matrix.csv")
    basel_lysis = pd.Series(dtype=float)
    if basel_path.exists():
        basel = pd.read_csv(basel_path).dropna(subset=["interaction"])
        basel_lysis = basel.groupby("phage")["interaction"].mean()
    phage_lysis_rate = pd.concat([guelin_lysis, basel_lysis])

    summary: dict[str, dict[str, object]] = {}
    subsets = {
        "guelin_n=%d" % len(guelin): guelin,
        "basel_nonzero_proj_n=%d" % len(basel_nonzero_proj): basel_nonzero_proj,
        "basel_zero_proj_n=%d" % len(basel_zero_proj): basel_zero_proj,
    }

    for label, phages in subsets.items():
        sub = slot_df[slot_df["phage"].isin(phages)]
        if sub.empty:
            summary[label] = {"status": "empty_subset"}
            continue
        mat = sub[feature_cols].to_numpy(dtype=float)
        mean = mat.mean(axis=0)
        std = mat.std(axis=0, ddof=1)
        cv = np.where(np.abs(mean) > 1e-10, std / np.abs(mean), np.nan)

        lysis = sub["phage"].map(phage_lysis_rate).fillna(phage_lysis_rate.median()).to_numpy()
        median_lysis = np.median(lysis)
        hi = mat[lysis > median_lysis]
        lo = mat[lysis <= median_lysis]
        if hi.size and lo.size:
            pooled = np.sqrt((hi.var(axis=0, ddof=1) + lo.var(axis=0, ddof=1)) / 2)
            cohens_d = np.where(pooled > 1e-10, (hi.mean(axis=0) - lo.mean(axis=0)) / pooled, 0.0)
        else:
            cohens_d = np.zeros(len(feature_cols))

        max_abs_cv = float(np.nanmax(np.abs(cv)))
        max_abs_d = float(np.max(np.abs(cohens_d)))
        summary[label] = {
            "n_phages": len(sub),
            "n_features": len(feature_cols),
            "max_abs_cv": round(max_abs_cv, 4),
            "max_abs_cohens_d": round(max_abs_d, 4),
            "n_features_with_cv_gt_0.1": int(np.nansum(np.abs(cv) > 0.1)),
            "n_features_with_cohens_d_gt_0.1": int(np.sum(np.abs(cohens_d) > 0.1)),
            "gate_status": "pass" if (max_abs_cv > 0.1 or max_abs_d > 0.1) else "fail",
        }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info("Arm 4 variance pre-flight: %s", json.dumps(summary, indent=2))
    return summary


def run_arm4_training_eval(
    *,
    slot_csv: Path,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
    device_type: str = "cpu",
    max_folds: int | None = None,
    num_workers: int = 3,
    drop_high_titer_only_positives: bool = False,
) -> dict[str, object]:
    """Run CH05 two-axis eval with the Arm 4 phage_projection slot."""
    from lyzortx.pipeline.autoresearch.candidate_replay import load_module_from_path
    from lyzortx.pipeline.autoresearch.ch04_eval import (
        BOOTSTRAP_RANDOM_STATE,
        BOOTSTRAP_SAMPLES,
        build_clean_row_training_frame,
        select_pair_max_concentration_rows,
    )
    from lyzortx.pipeline.autoresearch.ch05_eval import (
        BASEL_LOG10_PFU_ML,
        SOURCE_BASEL,
        SOURCE_GUELIN,
        _bootstrap_by_unit,
        _ci_to_dict,
        load_unified_row_frame,
        run_bacteria_axis,
        run_phage_axis,
    )
    from lyzortx.pipeline.autoresearch.sx03_eval import patch_context_with_extended_slots
    from lyzortx.pipeline.autoresearch.sx15_eval import load_unified_phage_family_map

    output_dir.mkdir(parents=True, exist_ok=True)
    if not slot_csv.exists():
        raise FileNotFoundError(f"Arm 4 slot CSV missing: {slot_csv}")

    arm4_slot_dir = slot_csv.parent.parent
    arm4_slot_dir.mkdir(parents=True, exist_ok=True)
    for slot_name in ("phage_stats", "phage_rbp_struct"):
        arm_side = arm4_slot_dir / slot_name / "features.csv"
        baseline_side = Path(".scratch/basel/feature_slots") / slot_name / "features.csv"
        if not arm_side.exists() and baseline_side.exists():
            arm_side.parent.mkdir(parents=True, exist_ok=True)
            arm_side.symlink_to(baseline_side.resolve())

    unified = load_unified_row_frame(basel_log10_pfu_ml=BASEL_LOG10_PFU_ML)
    clean_rows = build_clean_row_training_frame(unified, drop_high_titer_only_positives=drop_high_titer_only_positives)
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]
    phage_family = load_unified_phage_family_map()
    candidate_module = load_module_from_path("ch06_arm4_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context, slots_dir=arm4_slot_dir)

    bacteria_per_row = run_bacteria_axis(
        candidate_module=candidate_module,
        context=context,
        candidate_dir=candidate_dir,
        clean_rows=clean_rows,
        device_type=device_type,
        output_dir=output_dir,
        max_folds=max_folds,
        num_workers=num_workers,
    )
    phage_per_row = run_phage_axis(
        candidate_module=candidate_module,
        context=context,
        candidate_dir=candidate_dir,
        clean_rows=clean_rows,
        phage_family=phage_family,
        device_type=device_type,
        output_dir=output_dir,
        max_folds=max_folds,
        num_workers=num_workers,
    )

    bacteria_pairs = select_pair_max_concentration_rows(bacteria_per_row)
    bacteria_pairs["source"] = bacteria_pairs["pair_id"].map(pair_source)
    bacteria_pairs.to_csv(output_dir / "ch06_arm4_bacteria_axis_predictions.csv", index=False)
    bacteria_cis = _bootstrap_by_unit(
        bacteria_pairs.to_dict(orient="records"),
        unit_key="bacteria",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    phage_pairs = select_pair_max_concentration_rows(phage_per_row)
    phage_pairs["source"] = phage_pairs["pair_id"].map(pair_source)
    phage_pairs.to_csv(output_dir / "ch06_arm4_phage_axis_predictions.csv", index=False)
    phage_rows = phage_pairs.to_dict(orient="records")
    phage_cis_all = _bootstrap_by_unit(
        phage_rows,
        unit_key="phage",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    cross_source_rows: list[dict[str, object]] = []
    for source_label in (SOURCE_GUELIN, SOURCE_BASEL):
        subset = [r for r in phage_rows if r["source"] == source_label]
        if not subset:
            raise ValueError(f"Empty cross-source subset for {source_label!r}")
        subset_cis = _bootstrap_by_unit(
            subset,
            unit_key="phage",
            bootstrap_samples=BOOTSTRAP_SAMPLES,
            bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
        )
        cross_source_rows.append(
            {
                "source": source_label,
                "n_pairs": len(subset),
                "n_phages": len({r["phage"] for r in subset}),
                "auc": subset_cis["holdout_roc_auc"].point_estimate,
                "auc_low": subset_cis["holdout_roc_auc"].ci_low,
                "auc_high": subset_cis["holdout_roc_auc"].ci_high,
                "brier": subset_cis["holdout_brier_score"].point_estimate,
                "brier_low": subset_cis["holdout_brier_score"].ci_low,
                "brier_high": subset_cis["holdout_brier_score"].ci_high,
            }
        )
    pd.DataFrame(cross_source_rows).to_csv(output_dir / "ch06_arm4_cross_source_breakdown.csv", index=False)

    summary = {
        "arm": "CH06 Arm 4 — tail-restricted TL17 BLAST",
        "slot_csv": str(slot_csv),
        "bacteria_axis": {
            "aggregate": {name: _ci_to_dict(ci) for name, ci in bacteria_cis.items()},
            "n_pairs": len(bacteria_pairs),
        },
        "phage_axis": {
            "aggregate": {name: _ci_to_dict(ci) for name, ci in phage_cis_all.items()},
            "cross_source": cross_source_rows,
            "n_pairs": len(phage_pairs),
        },
    }
    with (output_dir / "ch06_arm4_metrics.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    LOGGER.info("CH06 Arm 4 summary written to %s", output_dir / "ch06_arm4_metrics.json")
    return summary


def _load_guelin_panel() -> set[str]:
    """Load the 96-phage Guelin panel list from the combined_queries.faa header order."""
    combined = Path("lyzortx/generated_outputs/autoresearch/phage_projection_cache_build/_batched/combined_queries.faa")
    panel: set[str] = set()
    with combined.open("r", encoding="utf-8") as fh:
        for line in fh:
            if line.startswith(">"):
                panel.add(line[1:].split("|", 1)[0])
    return panel


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--guelin-pharokka-dir", type=Path, default=DEFAULT_GUELIN_PHAROKKA_DIR)
    parser.add_argument("--guelin-protein-dir", type=Path, default=DEFAULT_GUELIN_PROTEIN_DIR)
    parser.add_argument("--basel-pharokka-dir", type=Path, default=DEFAULT_BASEL_PHAROKKA_DIR)
    parser.add_argument("--reference-fasta", type=Path, default=DEFAULT_REFERENCE_FASTA)
    parser.add_argument("--original-slot-csv", type=Path, default=DEFAULT_ORIGINAL_SLOT_CSV)
    parser.add_argument("--slot-csv", type=Path, default=DEFAULT_SLOT_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH_DIR)
    parser.add_argument(
        "--interactions-path",
        type=Path,
        default=Path("data/interactions/raw/raw_interactions.csv"),
    )
    parser.add_argument(
        "--runtime-payload",
        type=Path,
        default=Path("lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/tl17_rbp_runtime.joblib"),
    )
    parser.add_argument("--skip-preflight", action="store_true")
    parser.add_argument("--run-training", action="store_true")
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("lyzortx/generated_outputs/autoresearch/search_cache_v1"),
    )
    parser.add_argument("--candidate-dir", type=Path, default=Path("lyzortx/autoresearch"))
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--max-folds", type=int, default=None)
    parser.add_argument("--num-workers", type=int, default=3)
    parser.add_argument("--force", action="store_true")
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    args.scratch_dir.mkdir(parents=True, exist_ok=True)

    guelin_panel = _load_guelin_panel()
    LOGGER.info("Guelin panel: %d phages", len(guelin_panel))

    combined_fasta = args.scratch_dir / "tail_restricted_query.faa"
    per_phage_counts = build_tail_restricted_query_fasta(
        args.guelin_pharokka_dir,
        args.guelin_protein_dir,
        args.basel_pharokka_dir,
        guelin_panel,
        combined_fasta,
    )

    preflight_gate = {
        "guelin_tail_annotated": sum(1 for p, c in per_phage_counts.items() if not p.startswith("Bas") and c > 0),
        "guelin_total": sum(1 for p in per_phage_counts if not p.startswith("Bas")),
        "basel_tail_annotated": sum(1 for p, c in per_phage_counts.items() if p.startswith("Bas") and c > 0),
        "basel_total": sum(1 for p in per_phage_counts if p.startswith("Bas")),
    }
    LOGGER.info("Tail-annotation coverage: %s", preflight_gate)
    if preflight_gate["basel_total"] > 0:
        basel_coverage = preflight_gate["basel_tail_annotated"] / preflight_gate["basel_total"]
        if basel_coverage < 0.5:
            LOGGER.error(
                "BASEL pre-flight gate FAILED: only %.1f%% (%d/%d) of BASEL phages have tail annotations; "
                "arm would zero out the BASEL side. Skip.",
                basel_coverage * 100,
                preflight_gate["basel_tail_annotated"],
                preflight_gate["basel_total"],
            )
            sys.exit(2)

    hits_tsv = args.scratch_dir / "mmseqs_hits.tsv"
    if args.force or not hits_tsv.exists():
        run_mmseqs_search(combined_fasta, args.reference_fasta, hits_tsv, args.scratch_dir)
    else:
        LOGGER.info("Reusing cached mmseqs hits: %s", hits_tsv)

    feature_df = aggregate_tl17_features(
        hits_tsv,
        args.runtime_payload,
        sorted(per_phage_counts.keys()),
    )
    slot_df = materialize_arm4_slot(feature_df, args.original_slot_csv, args.slot_csv)

    if not args.skip_preflight:
        run_variance_preflight(
            slot_df,
            args.original_slot_csv,
            args.interactions_path,
            args.output_dir / "ch06_arm4_variance_preflight.json",
        )

    if args.run_training:
        run_arm4_training_eval(
            slot_csv=args.slot_csv,
            output_dir=args.output_dir,
            cache_dir=args.cache_dir,
            candidate_dir=args.candidate_dir,
            device_type=args.device_type,
            max_folds=args.max_folds,
            num_workers=args.num_workers,
        )


if __name__ == "__main__":
    main()
