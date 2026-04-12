#!/usr/bin/env python3
"""SX02: BASEL phage feature computation.

Runs the full annotation and feature pipeline on 52 BASEL phage genomes:
  1. Pharokka annotation (must be pre-run, checked at startup)
  2. DepoScope depolymerase prediction
  3. Phage stats features (GC, length, N50, record count)
  4. Extend slot CSVs with BASEL phage rows

Produces extended slot CSVs at .scratch/basel/feature_slots/ that can be
loaded alongside the existing Guelin cache for SX03 integration.

Usage:
    python -m lyzortx.pipeline.autoresearch.sx02_basel_features
"""

from __future__ import annotations

import logging
import time
from pathlib import Path

import pandas as pd
import torch

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.run_deposcope import load_model, predict_phage

LOGGER = logging.getLogger(__name__)

BASEL_GENOMES_DIR = Path(".scratch/basel/genomes")
PHAROKKA_OUTPUT_DIR = Path(".scratch/basel/pharokka_output")
OUTPUT_DIR = Path(".scratch/basel/feature_slots")
DEPOSCOPE_PREDICTIONS_PATH = Path(".scratch/basel/deposcope_predictions.csv")
EXISTING_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1/feature_slots")


def check_pharokka_complete() -> list[str]:
    """Verify Pharokka ran on all BASEL genomes. Returns list of phage names."""
    genome_files = sorted(BASEL_GENOMES_DIR.glob("*.fna"))
    if not genome_files:
        raise FileNotFoundError(f"No genome files in {BASEL_GENOMES_DIR}")

    phage_names = []
    for g in genome_files:
        name = g.stem
        outdir = PHAROKKA_OUTPUT_DIR / name
        gff = outdir / f"{name}.gff"
        if not gff.exists():
            raise FileNotFoundError(f"Pharokka output missing for {name}: expected {gff}")
        phage_names.append(name)

    LOGGER.info("Pharokka output verified for %d BASEL phages", len(phage_names))
    return phage_names


def find_protein_fasta(phage_name: str) -> Path:
    """Find the protein FASTA from Pharokka output."""
    outdir = PHAROKKA_OUTPUT_DIR / phage_name
    faa_files = list(outdir.glob("*.faa"))
    if not faa_files:
        raise FileNotFoundError(f"No .faa files for {phage_name} in {outdir}")
    # Prefer the CDS FASTA (not tRNA etc.)
    for f in faa_files:
        if "cds" in f.name.lower() or phage_name in f.name:
            return f
    return faa_files[0]


def run_deposcope_on_basel(phage_names: list[str]) -> pd.DataFrame:
    """Run DepoScope on all BASEL phage proteins. Returns predictions DataFrame."""
    if DEPOSCOPE_PREDICTIONS_PATH.exists():
        LOGGER.info("Loading cached DepoScope predictions from %s", DEPOSCOPE_PREDICTIONS_PATH)
        return pd.read_csv(DEPOSCOPE_PREDICTIONS_PATH)

    device = torch.device("mps" if torch.backends.mps.is_available() else "cpu")
    LOGGER.info("Loading DepoScope model on %s...", device)
    model = load_model(device)

    all_results: list[dict] = []
    t0 = time.time()
    total_proteins = 0

    for i, phage in enumerate(phage_names):
        fasta_path = find_protein_fasta(phage)
        proteins = predict_phage(model, fasta_path)
        for p in proteins:
            p["phage"] = phage
        all_results.extend(proteins)
        total_proteins += len(proteins)

        depo_count = sum(1 for p in proteins if p["is_depolymerase"])
        if (i + 1) % 10 == 0 or i + 1 == len(phage_names):
            elapsed = time.time() - t0
            LOGGER.info(
                "DepoScope progress: %d/%d phages, %d proteins, %.0fs elapsed",
                i + 1,
                len(phage_names),
                total_proteins,
                elapsed,
            )
        if depo_count > 0:
            LOGGER.info("  %s: %d depolymerases found", phage, depo_count)

    df = pd.DataFrame(all_results)
    df.to_csv(DEPOSCOPE_PREDICTIONS_PATH, index=False)
    LOGGER.info(
        "DepoScope complete: %d proteins, %d depolymerases across %d phages",
        len(df),
        df["is_depolymerase"].sum() if "is_depolymerase" in df.columns else 0,
        len(phage_names),
    )
    return df


def compute_phage_stats(phage_names: list[str]) -> pd.DataFrame:
    """Compute phage_stats features for BASEL phages."""
    rows = []
    for name in phage_names:
        genome_path = BASEL_GENOMES_DIR / f"{name}.fna"
        text = genome_path.read_text()
        lines = text.strip().split("\n")
        records = []
        current_seq: list[str] = []
        for line in lines:
            if line.startswith(">"):
                if current_seq:
                    records.append("".join(current_seq))
                current_seq = []
            else:
                current_seq.append(line.strip())
        if current_seq:
            records.append("".join(current_seq))

        total_seq = "".join(records)
        gc = (total_seq.upper().count("G") + total_seq.upper().count("C")) / len(total_seq) if total_seq else 0

        lengths = sorted([len(r) for r in records], reverse=True)
        cumsum = 0
        n50 = lengths[0] if lengths else 0
        for length in lengths:
            cumsum += length
            if cumsum >= len(total_seq) / 2:
                n50 = length
                break

        rows.append(
            {
                "phage": name,
                "phage_stats__phage_gc_content": round(gc, 6),
                "phage_stats__phage_genome_length_nt": len(total_seq),
                "phage_stats__phage_n50_contig_length_nt": n50,
                "phage_stats__phage_sequence_record_count": len(records),
            }
        )

    return pd.DataFrame(rows)


def build_combined_deposcope_files(depo_predictions: pd.DataFrame) -> dict[str, int]:
    """Combine Guelin + BASEL DepoScope predictions into unified files.

    The pairwise cross-term code (derive_pairwise_depo_capsule_features.py)
    reads from .scratch/deposcope/predictions.csv. We create a combined file
    so cross-terms work for both Guelin and BASEL phages at training time.

    Returns per-phage depo counts for BASEL phages.
    """
    guelin_preds_path = Path(".scratch/deposcope/predictions.csv")
    combined_dir = Path(".scratch/deposcope_combined")
    combined_dir.mkdir(parents=True, exist_ok=True)

    # Combine predictions.
    guelin_preds = pd.read_csv(guelin_preds_path)
    combined_preds = pd.concat([guelin_preds, depo_predictions], ignore_index=True)
    combined_preds.to_csv(combined_dir / "predictions.csv", index=False)
    LOGGER.info(
        "Combined DepoScope predictions: %d Guelin + %d BASEL = %d total",
        len(guelin_preds),
        len(depo_predictions),
        len(combined_preds),
    )

    # Copy cluster file (BASEL depos don't have cluster assignments yet —
    # cross-terms will use has_depo and depo_count for BASEL phages).
    guelin_clusters_path = Path(".scratch/deposcope/depo_clusters_cluster.tsv")
    if guelin_clusters_path.exists():
        import shutil

        shutil.copy2(guelin_clusters_path, combined_dir / "depo_clusters_cluster.tsv")

    # Per-phage depo counts for BASEL.
    basel_depos = (
        depo_predictions[depo_predictions["is_depolymerase"]].copy() if len(depo_predictions) > 0 else pd.DataFrame()
    )
    phage_depo_counts = basel_depos.groupby("phage").size().to_dict() if len(basel_depos) > 0 else {}
    return phage_depo_counts


def compute_rbp_struct_features(phage_names: list[str]) -> pd.DataFrame:
    """Compute phage_rbp_struct features for BASEL phages.

    Sets has_annotated_rbp and rbp_count from Pharokka merged output.
    PLM PCA features are zero-filled (would need ProstT5+SaProt inference).
    """
    rows = []
    for name in phage_names:
        outdir = PHAROKKA_OUTPUT_DIR / name
        tsv_files = list(outdir.glob("*_cds_final_merged_output.tsv"))
        rbp_count = 0
        if tsv_files:
            import csv

            with open(tsv_files[0], encoding="utf-8") as f:
                reader = csv.DictReader(f, delimiter="\t")
                for row in reader:
                    annot = (row.get("annot", "") + " " + row.get("phrog", "")).lower()
                    if any(kw in annot for kw in ("tail fiber", "tail spike", "receptor binding", "rbp")):
                        rbp_count += 1
        rows.append({"phage": name, "rbp_count": rbp_count, "has_annotated_rbp": int(rbp_count > 0)})

    df = pd.DataFrame(rows)
    LOGGER.info("RBP struct: %d phages, %d with annotated RBPs", len(df), df["has_annotated_rbp"].sum())
    return df


def verify_feature_distributions(guelin_df: pd.DataFrame, basel_df: pd.DataFrame, slot_name: str) -> None:
    """Log CV and unique value counts comparing BASEL vs Guelin features."""

    numeric_cols = [c for c in guelin_df.columns if c != "phage" and guelin_df[c].dtype in ("float64", "int64")]
    if not numeric_cols:
        LOGGER.info("  %s: no numeric columns to compare", slot_name)
        return

    LOGGER.info("Feature distribution verification for %s (%d numeric columns):", slot_name, len(numeric_cols))
    for col in numeric_cols[:10]:  # Sample first 10 columns.
        g_vals = guelin_df[col].dropna()
        b_vals = basel_df[col].dropna() if col in basel_df.columns else pd.Series(dtype=float)
        g_cv = float(g_vals.std() / g_vals.mean()) if g_vals.mean() != 0 and len(g_vals) > 1 else 0
        b_cv = float(b_vals.std() / b_vals.mean()) if len(b_vals) > 1 and b_vals.mean() != 0 else 0
        g_unique = int(g_vals.nunique())
        b_unique = int(b_vals.nunique()) if len(b_vals) > 0 else 0
        LOGGER.info("  %s: Guelin CV=%.3f (%d unique), BASEL CV=%.3f (%d unique)", col, g_cv, g_unique, b_cv, b_unique)

    # Overall summary.
    guelin_nonzero = sum(1 for c in numeric_cols if guelin_df[c].abs().sum() > 0)
    basel_nonzero = sum(1 for c in numeric_cols if c in basel_df.columns and basel_df[c].abs().sum() > 0)
    LOGGER.info(
        "  %s summary: %d/%d columns nonzero in Guelin, %d/%d in BASEL",
        slot_name,
        guelin_nonzero,
        len(numeric_cols),
        basel_nonzero,
        len(numeric_cols),
    )


def write_extended_slot(
    slot_name: str,
    guelin_df: pd.DataFrame,
    basel_df: pd.DataFrame,
    entity_key: str = "phage",
) -> Path:
    """Write extended slot CSV combining Guelin + BASEL rows."""
    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)
    slot_dir = OUTPUT_DIR / slot_name
    slot_dir.mkdir(parents=True, exist_ok=True)

    # Align columns — BASEL may have fewer columns, fill with 0.
    for col in guelin_df.columns:
        if col not in basel_df.columns and col != entity_key:
            basel_df[col] = 0.0

    # Only keep columns present in Guelin.
    common_cols = [entity_key] + [c for c in guelin_df.columns if c != entity_key]
    basel_aligned = basel_df.reindex(columns=common_cols, fill_value=0.0)

    combined = pd.concat([guelin_df, basel_aligned], ignore_index=True)
    out_path = slot_dir / "features.csv"
    combined.to_csv(out_path, index=False)
    LOGGER.info(
        "Extended %s: %d Guelin + %d BASEL = %d total rows", slot_name, len(guelin_df), len(basel_df), len(combined)
    )

    # Also write entity index.
    combined[[entity_key]].to_csv(slot_dir / "entity_index.csv", index=False)
    return out_path


def main() -> None:
    setup_logging()
    LOGGER.info("SX02: BASEL phage feature computation starting")

    # Step 1: Verify Pharokka annotations.
    phage_names = check_pharokka_complete()

    # Step 2: Run DepoScope.
    depo_predictions = run_deposcope_on_basel(phage_names)

    # Step 3: Compute phage_stats.
    LOGGER.info("Computing phage_stats features...")
    stats_df = compute_phage_stats(phage_names)
    LOGGER.info("Phage stats: %d phages, columns: %s", len(stats_df), list(stats_df.columns))

    # Step 4: Build combined DepoScope files for pairwise cross-terms.
    phage_depo_counts = build_combined_deposcope_files(depo_predictions)
    LOGGER.info(
        "BASEL phages with depolymerases: %d/%d", sum(1 for v in phage_depo_counts.values() if v > 0), len(phage_names)
    )

    # Step 5: Compute RBP struct features.
    LOGGER.info("Computing phage_rbp_struct features...")
    rbp_df = compute_rbp_struct_features(phage_names)

    # Step 6: Extend slot CSVs.
    # phage_stats
    guelin_stats = pd.read_csv(EXISTING_CACHE_DIR / "phage_stats" / "features.csv")
    verify_feature_distributions(guelin_stats, stats_df, "phage_stats")
    write_extended_slot("phage_stats", guelin_stats, stats_df)

    # phage_projection — BASEL phages get zero-filled RBP family features.
    guelin_proj = pd.read_csv(EXISTING_CACHE_DIR / "phage_projection" / "features.csv")
    basel_proj = pd.DataFrame({"phage": phage_names})
    for col in guelin_proj.columns:
        if col != "phage":
            basel_proj[col] = 0.0
    verify_feature_distributions(guelin_proj, basel_proj, "phage_projection")
    write_extended_slot("phage_projection", guelin_proj, basel_proj)

    # phage_rbp_struct — rbp_count and has_rbp from Pharokka, zero-fill PLM PCA.
    guelin_rbp = pd.read_csv(EXISTING_CACHE_DIR / "phage_rbp_struct" / "features.csv")
    # Map BASEL rbp features to slot column names.
    basel_rbp = pd.DataFrame({"phage": rbp_df["phage"]})
    basel_rbp["phage_rbp_struct__has_annotated_rbp"] = rbp_df["has_annotated_rbp"]
    basel_rbp["phage_rbp_struct__rbp_count"] = rbp_df["rbp_count"]
    verify_feature_distributions(guelin_rbp, basel_rbp, "phage_rbp_struct")
    write_extended_slot("phage_rbp_struct", guelin_rbp, basel_rbp)

    # Summary.
    depo_phages = sum(1 for v in phage_depo_counts.values() if v > 0)
    LOGGER.info("=" * 60)
    LOGGER.info("SX02 BASEL Feature Computation Summary")
    LOGGER.info("  Phages annotated: %d", len(phage_names))
    LOGGER.info("  Phages with DepoScope depolymerases: %d/%d", depo_phages, len(phage_names))
    LOGGER.info("  Phages with annotated RBPs: %d/%d", rbp_df["has_annotated_rbp"].sum(), len(phage_names))
    LOGGER.info("  Combined DepoScope predictions at: .scratch/deposcope_combined/")
    LOGGER.info("  Extended slots (3): phage_stats, phage_projection, phage_rbp_struct")
    LOGGER.info("  Extended slots written to: %s", OUTPUT_DIR)
    LOGGER.info("=" * 60)


if __name__ == "__main__":
    main()
