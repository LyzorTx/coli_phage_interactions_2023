#!/usr/bin/env python3
"""CH06 Arm 2: MMseqs2 pairwise proteome similarity as phage-side feature slot.

Replaces the Guelin-only TL17 categorical projection (33-col phage_projection slot)
with a continuous 148-dim pairwise proteome similarity vector reduced to 32 PCA
components. Each phage's feature row is its MMseqs2 all-vs-all bit-score profile
against every other phage in the unified 148-phage panel (96 Guelin + 52 BASEL),
then projected onto PCA axes fit on the same matrix.

A BASEL phage that has no close Guelin neighbour still produces a real (mostly
low) similarity profile rather than TL17's zero-vector failure mode — that is the
plan.yml hypothesis Arm 2 tests.

Two-phase computation:

1. Precompute (`build_similarity_matrix`): extract proteins from all 148 phages
   into a single FASTA with `<phage_id>|<protein_id>` headers, run MMseqs2
   easy-search all-vs-all (e-value 1e-5), aggregate hits to a 148×148 phage-level
   sum-bit-score matrix. Self-similarity (diagonal) is zeroed so `sim_to_<self>`
   is not a trivial constant feature that could encode phage identity.
   Cached to `.scratch/ch06_arm2/pairwise_bitscore.npz`.

2. Slot materializer (`materialize_arm2_slot`): PCA(n_components=32) fit on the
   148 phage rows, write a feature CSV with 32 `phage_projection__arm2_pc_<k>`
   columns (reusing the `phage_projection` slot name so the existing CH05
   patch-in-place path picks it up without any feature-bundle refactor).

3. Variance pre-flight (`run_variance_preflight`): for each PCA component, compute
   CV across (a) Guelin phages, (b) non-zero-projection BASEL phages, (c) other
   BASEL phages, plus Cohen's d for discriminating lysed vs non-lysed pairs at
   the phage level. Fails loudly if every component has CV < 0.1 AND Cohen's
   d < 0.1 on all three subsets per plan.yml variance-preflight gate.

4. Eval driver (`run_arm2_eval`): patches the `phage_projection` slot with the
   Arm 2 CSV and runs CH05's `run_bacteria_axis` + `run_phage_axis` under the
   same parallel loop as the baseline.
"""

from __future__ import annotations

import argparse
import json
import logging
import subprocess
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd
from sklearn.decomposition import PCA

from lyzortx.log_config import setup_logging

LOGGER = logging.getLogger(__name__)

GUELIN_COMBINED_FAA = Path(
    "lyzortx/generated_outputs/autoresearch/phage_projection_cache_build/_batched/combined_queries.faa"
)
BASEL_PHAROKKA_DIR = Path(".scratch/basel/pharokka_output")
DEFAULT_SCRATCH_DIR = Path(".scratch/ch06_arm2")
DEFAULT_SLOT_CSV = Path(".scratch/basel/feature_slots_arm2/phage_projection/features.csv")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch06_arm2_mmseqs_proteome")
DEFAULT_ORIGINAL_SLOT_CSV = Path(".scratch/basel/feature_slots/phage_projection/features.csv")

MMSEQS_E_VALUE = 1e-5
MMSEQS_SENSITIVITY = 7.5
PCA_COMPONENTS = 32
COLUMN_PREFIX = "phage_projection__arm2_pc_"


def build_combined_proteome(output_faa: Path) -> dict[str, int]:
    """Concatenate Guelin + BASEL proteins into a single FASTA with `<phage>|<prot>` headers.

    Guelin source (`combined_queries.faa`) already has that header shape (`409_P1|query_prot_0001`).
    BASEL sources (`.scratch/basel/pharokka_output/<Bas##>/phanotate.faa`) use pharokka's
    `<contig_id>_CDS_<idx>` headers, which we rewrite to `<Bas##>|<contig_id>_CDS_<idx>`.
    """
    output_faa.parent.mkdir(parents=True, exist_ok=True)
    per_phage_counts: dict[str, int] = {}
    with output_faa.open("w") as out:
        for line in GUELIN_COMBINED_FAA.read_text().splitlines():
            out.write(line + "\n")
            if line.startswith(">"):
                phage_id = line[1:].split("|", 1)[0]
                per_phage_counts[phage_id] = per_phage_counts.get(phage_id, 0) + 1

        for basel_dir in sorted(BASEL_PHAROKKA_DIR.iterdir()):
            if not basel_dir.is_dir() or not basel_dir.name.startswith("Bas"):
                continue
            phage_id = basel_dir.name
            faa = basel_dir / "phanotate.faa"
            if not faa.exists():
                raise FileNotFoundError(f"BASEL pharokka faa missing: {faa}")
            n = 0
            for line in faa.read_text().splitlines():
                if line.startswith(">"):
                    protein_id = line[1:].split()[0]
                    out.write(f">{phage_id}|{protein_id}\n")
                    n += 1
                else:
                    out.write(line + "\n")
            per_phage_counts[phage_id] = n

    LOGGER.info(
        "Combined proteome: %d phages, %d proteins total, written to %s",
        len(per_phage_counts),
        sum(per_phage_counts.values()),
        output_faa,
    )
    return per_phage_counts


def run_mmseqs_all_vs_all(
    combined_faa: Path,
    output_tsv: Path,
    scratch_dir: Path,
    *,
    e_value: float = MMSEQS_E_VALUE,
    sensitivity: float = MMSEQS_SENSITIVITY,
) -> None:
    """Run MMseqs2 `easy-search` self-vs-self on `combined_faa` → protein-pair TSV."""
    output_tsv.parent.mkdir(parents=True, exist_ok=True)
    tmp_dir = scratch_dir / "mmseqs_tmp"
    tmp_dir.mkdir(parents=True, exist_ok=True)

    # Output columns: query, target, bits (we don't need the full BLAST 12).
    cmd = [
        "mmseqs",
        "easy-search",
        str(combined_faa),
        str(combined_faa),
        str(output_tsv),
        str(tmp_dir),
        "--format-output",
        "query,target,bits",
        "-e",
        str(e_value),
        "-s",
        str(sensitivity),
        "--threads",
        str(max(1, (_cpu_count()))),
    ]
    LOGGER.info("Running MMseqs2 easy-search: %s", " ".join(cmd))
    subprocess.run(cmd, check=True)
    LOGGER.info("MMseqs2 hits written to %s", output_tsv)


def _cpu_count() -> int:
    import os

    return os.cpu_count() or 4


def aggregate_protein_hits_to_phage_matrix(
    hits_tsv: Path,
    phage_order: Sequence[str],
) -> np.ndarray:
    """Sum bit-scores of all protein-level hits per (query_phage, target_phage) pair.

    MMseqs2 `easy-search` returns top-N hits per query (default 300). Summing all
    reported hits is a proxy for "how much protein-space does phage A share with
    phage B"; any single mutually-conserved core gene contributes once to each of
    its per-query rows. Diagonal (self-similarity) is zeroed below to avoid
    encoding phage identity.
    """
    LOGGER.info("Reading MMseqs2 hits from %s", hits_tsv)
    hits = pd.read_csv(
        hits_tsv,
        sep="\t",
        header=None,
        names=["query", "target", "bits"],
        dtype={"query": str, "target": str, "bits": float},
    )
    LOGGER.info("Parsed %d hits", len(hits))

    hits["query_phage"] = hits["query"].str.split("|", n=1).str[0]
    hits["target_phage"] = hits["target"].str.split("|", n=1).str[0]

    pairwise = hits.groupby(["query_phage", "target_phage"])["bits"].sum().reset_index(name="bits_sum")

    phage_to_idx = {p: i for i, p in enumerate(phage_order)}
    n = len(phage_order)
    matrix = np.zeros((n, n), dtype=np.float32)

    known_mask = pairwise["query_phage"].isin(phage_to_idx) & pairwise["target_phage"].isin(phage_to_idx)
    dropped = (~known_mask).sum()
    if dropped:
        LOGGER.warning("Dropped %d hits with unknown phage ids", dropped)
    pairwise = pairwise[known_mask]

    q_idx = pairwise["query_phage"].map(phage_to_idx).to_numpy()
    t_idx = pairwise["target_phage"].map(phage_to_idx).to_numpy()
    matrix[q_idx, t_idx] = pairwise["bits_sum"].to_numpy(dtype=np.float32)

    np.fill_diagonal(matrix, 0.0)

    nonzero_fraction = (matrix > 0).mean()
    LOGGER.info(
        "Pairwise bit-score matrix: %d x %d, %.1f%% nonzero (excluding diagonal), max=%.0f",
        n,
        n,
        nonzero_fraction * 100,
        matrix.max(),
    )
    return matrix


def materialize_arm2_slot(
    phage_order: Sequence[str],
    pairwise_matrix: np.ndarray,
    slot_csv: Path,
    *,
    n_components: int = PCA_COMPONENTS,
) -> pd.DataFrame:
    """PCA to 32 components, write a phage_projection-shaped CSV with 32 PCA columns.

    The slot reuses the `phage_projection__` column prefix so CH05's existing
    `patch_context_with_extended_slots` path picks it up without any
    feature-bundle refactor. Column names are stable and named so LightGBM's
    feature-importance dump distinguishes Arm 2 from Guelin TL17 slots when
    compared across runs.
    """
    slot_csv.parent.mkdir(parents=True, exist_ok=True)

    pca = PCA(n_components=n_components, random_state=0)
    projected = pca.fit_transform(pairwise_matrix)

    columns = [f"{COLUMN_PREFIX}{k:02d}" for k in range(n_components)]
    df = pd.DataFrame(projected, columns=columns)
    df.insert(0, "phage", list(phage_order))
    df.to_csv(slot_csv, index=False)

    explained = pca.explained_variance_ratio_
    LOGGER.info(
        "Arm 2 PCA slot written to %s; top-5 components explain %.1f%%, top-32 explain %.1f%%",
        slot_csv,
        explained[:5].sum() * 100,
        explained.sum() * 100,
    )
    return df


def run_variance_preflight(
    arm2_slot_df: pd.DataFrame,
    original_slot_csv: Path,
    interactions_path: Path,
    output_json: Path,
) -> dict[str, object]:
    """Per plan.yml VARIANCE PRE-FLIGHT — CV and Cohen's d on three phage subsets.

    Three subsets (per plan.yml):
      (1) Guelin phages (n=96)
      (2) BASEL phages with non-zero TL17 phage_projection (n≈39, picks up Guelin-
          neighbour-carrying BASEL phages)
      (3) BASEL phages with zero-vector TL17 projection (n≈13)

    For each subset, compute per-component CV across phages. For each phage,
    aggregate the interaction matrix to its lysis rate across 369 clinical
    bacteria, then compute Cohen's d on each PCA component split at median lysis
    rate within the subset.

    Gate: if ALL 32 components have CV < 0.1 AND Cohen's d < 0.1 on ALL three
    subsets, the arm is killed before training. Because Arm 2 feature values
    literally encode proteome similarity to specific phages, it is extremely
    unlikely to have degenerate variance — the pre-flight exists to catch the
    pathological case, not to screen out wash arms.
    """
    original = pd.read_csv(original_slot_csv)
    feature_cols = [c for c in original.columns if c != "phage"]
    original_nonzero = original[feature_cols].abs().sum(axis=1) > 0
    nonzero_phages = set(original.loc[original_nonzero, "phage"])

    guelin = [p for p in arm2_slot_df["phage"] if not p.startswith("Bas")]
    basel_all = [p for p in arm2_slot_df["phage"] if p.startswith("Bas")]
    basel_nonzero_proj = [p for p in basel_all if p in nonzero_phages]
    basel_zero_proj = [p for p in basel_all if p not in nonzero_phages]

    arm2_cols = [c for c in arm2_slot_df.columns if c.startswith(COLUMN_PREFIX)]
    interactions = pd.read_csv(interactions_path, sep=";", low_memory=False)
    interactions["score"] = interactions["score"].astype(str)
    guelin_labels = (interactions["score"] == "1").astype(int)
    guelin_lysis = (
        pd.DataFrame({"phage": interactions["phage"], "label": guelin_labels}).groupby("phage")["label"].mean()
    )
    basel_path = Path(".scratch/genophi_data/BASEL_ECOR_interaction_matrix.csv")
    basel_lysis = pd.Series(dtype=float)
    if basel_path.exists():
        basel = pd.read_csv(basel_path)
        basel = basel.dropna(subset=["interaction"])
        basel_lysis = basel.groupby("phage")["interaction"].mean()
    phage_lysis_rate = pd.concat([guelin_lysis, basel_lysis])

    summary: dict[str, dict[str, object]] = {}
    subsets = {
        "guelin_n=%d" % len(guelin): guelin,
        "basel_nonzero_proj_n=%d" % len(basel_nonzero_proj): basel_nonzero_proj,
        "basel_zero_proj_n=%d" % len(basel_zero_proj): basel_zero_proj,
    }

    for label, phages in subsets.items():
        sub = arm2_slot_df[arm2_slot_df["phage"].isin(phages)]
        if sub.empty:
            summary[label] = {"status": "empty_subset"}
            continue
        mat = sub[arm2_cols].to_numpy(dtype=float)
        mean = mat.mean(axis=0)
        std = mat.std(axis=0, ddof=1)
        cv = np.where(np.abs(mean) > 1e-10, std / np.abs(mean), np.nan)

        # Cohen's d per component across lysed vs non-lysed at phage level.
        lysis = sub["phage"].map(phage_lysis_rate).fillna(phage_lysis_rate.median()).to_numpy()
        median_lysis = np.median(lysis)
        hi = mat[lysis > median_lysis]
        lo = mat[lysis <= median_lysis]
        if hi.size and lo.size:
            pooled = np.sqrt((hi.var(axis=0, ddof=1) + lo.var(axis=0, ddof=1)) / 2)
            cohens_d = np.where(pooled > 1e-10, (hi.mean(axis=0) - lo.mean(axis=0)) / pooled, 0.0)
        else:
            cohens_d = np.zeros(len(arm2_cols))

        max_abs_cv = float(np.nanmax(np.abs(cv)))
        max_abs_d = float(np.max(np.abs(cohens_d)))
        summary[label] = {
            "n_phages": len(sub),
            "n_components": len(arm2_cols),
            "max_abs_cv_across_components": round(max_abs_cv, 4),
            "max_abs_cohens_d_across_components": round(max_abs_d, 4),
            "n_components_with_cv_gt_0.1": int(np.nansum(np.abs(cv) > 0.1)),
            "n_components_with_cohens_d_gt_0.1": int(np.sum(np.abs(cohens_d) > 0.1)),
            "gate_status": ("pass" if (max_abs_cv > 0.1 or max_abs_d > 0.1) else "fail"),
        }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info("Variance pre-flight summary: %s", json.dumps(summary, indent=2))

    any_fail = any(
        s.get("gate_status") == "fail" for s in summary.values() if isinstance(s, dict) and "gate_status" in s
    )
    if any_fail:
        LOGGER.warning("At least one subset failed the CV/Cohen's d gate — see pre-flight JSON")
    return summary


def build_similarity_matrix(
    scratch_dir: Path = DEFAULT_SCRATCH_DIR,
    *,
    force: bool = False,
) -> tuple[list[str], np.ndarray]:
    """End-to-end precompute: FASTA → MMseqs2 → 148×148 bit-score matrix."""
    scratch_dir.mkdir(parents=True, exist_ok=True)
    cache_path = scratch_dir / "pairwise_bitscore.npz"
    if cache_path.exists() and not force:
        cached = np.load(cache_path, allow_pickle=True)
        phage_order = list(cached["phage_order"])
        matrix = cached["matrix"]
        LOGGER.info("Loaded cached similarity matrix from %s (%d phages)", cache_path, len(phage_order))
        return phage_order, matrix

    combined_faa = scratch_dir / "combined_proteome.faa"
    hits_tsv = scratch_dir / "mmseqs_hits.tsv"
    per_phage = build_combined_proteome(combined_faa)
    phage_order = sorted(per_phage.keys())
    run_mmseqs_all_vs_all(combined_faa, hits_tsv, scratch_dir)
    matrix = aggregate_protein_hits_to_phage_matrix(hits_tsv, phage_order)
    np.savez(cache_path, phage_order=np.array(phage_order), matrix=matrix)
    LOGGER.info("Cached similarity matrix to %s", cache_path)
    return phage_order, matrix


def run_arm2_training_eval(
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
    """Run CH05's two-axis evaluation with the Arm 2 phage_projection slot.

    Patches `phage_projection` in the loaded context with the Arm 2 CSV, then
    invokes `run_bacteria_axis` + `run_phage_axis` directly (reusing CH05's
    parallel loop). Results land under `output_dir` with Arm 2 naming.
    """
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
        raise FileNotFoundError(f"Arm 2 slot CSV missing: {slot_csv}")

    # The Arm 2 slot directory must hold both the new phage_projection CSV AND the
    # original phage_stats / phage_rbp_struct CSVs (unchanged). Materialize symlinks
    # into the Arm 2 slot dir so patch_context_with_extended_slots finds all three.
    arm2_slot_dir = slot_csv.parent.parent
    arm2_slot_dir.mkdir(parents=True, exist_ok=True)
    for slot_name in ("phage_stats", "phage_rbp_struct"):
        arm2_side = arm2_slot_dir / slot_name / "features.csv"
        baseline_side = Path(".scratch/basel/feature_slots") / slot_name / "features.csv"
        if not arm2_side.exists() and baseline_side.exists():
            arm2_side.parent.mkdir(parents=True, exist_ok=True)
            arm2_side.symlink_to(baseline_side.resolve())
            LOGGER.info("Linked baseline slot %s → %s", baseline_side, arm2_side)

    unified = load_unified_row_frame(basel_log10_pfu_ml=BASEL_LOG10_PFU_ML)
    clean_rows = build_clean_row_training_frame(unified, drop_high_titer_only_positives=drop_high_titer_only_positives)
    LOGGER.info(
        "CH06 Arm 2 clean row frame: %d rows, %d pairs, %d bacteria, %d phages",
        len(clean_rows),
        clean_rows["pair_id"].nunique(),
        clean_rows["bacteria"].nunique(),
        clean_rows["phage"].nunique(),
    )
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]

    phage_family = load_unified_phage_family_map()
    candidate_module = load_module_from_path("ch06_arm2_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context, slots_dir=arm2_slot_dir)

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
    bacteria_pairs.to_csv(output_dir / "ch06_arm2_bacteria_axis_predictions.csv", index=False)
    bacteria_cis = _bootstrap_by_unit(
        bacteria_pairs.to_dict(orient="records"),
        unit_key="bacteria",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    phage_pairs = select_pair_max_concentration_rows(phage_per_row)
    phage_pairs["source"] = phage_pairs["pair_id"].map(pair_source)
    phage_pairs.to_csv(output_dir / "ch06_arm2_phage_axis_predictions.csv", index=False)
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
    pd.DataFrame(cross_source_rows).to_csv(output_dir / "ch06_arm2_cross_source_breakdown.csv", index=False)

    summary = {
        "arm": "CH06 Arm 2 — MMseqs2 pairwise proteome similarity",
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
    with (output_dir / "ch06_arm2_metrics.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    LOGGER.info("CH06 Arm 2 summary written to %s", output_dir / "ch06_arm2_metrics.json")
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scratch-dir", type=Path, default=DEFAULT_SCRATCH_DIR)
    parser.add_argument("--slot-csv", type=Path, default=DEFAULT_SLOT_CSV)
    parser.add_argument("--original-slot-csv", type=Path, default=DEFAULT_ORIGINAL_SLOT_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--interactions-path",
        type=Path,
        default=Path("data/interactions/raw/raw_interactions.csv"),
    )
    parser.add_argument("--force", action="store_true", help="Recompute MMseqs2 even if cached")
    parser.add_argument(
        "--skip-preflight",
        action="store_true",
        help="Skip the variance pre-flight (useful for re-running slot materialization only)",
    )
    parser.add_argument(
        "--run-training",
        action="store_true",
        help="After precompute+preflight, run the full CH05 two-axis evaluation with the Arm 2 slot.",
    )
    parser.add_argument(
        "--cache-dir",
        type=Path,
        default=Path("lyzortx/generated_outputs/autoresearch/search_cache_v1"),
    )
    parser.add_argument("--candidate-dir", type=Path, default=Path("lyzortx/autoresearch"))
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--max-folds", type=int, default=None)
    parser.add_argument("--num-workers", type=int, default=3)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    phage_order, matrix = build_similarity_matrix(args.scratch_dir, force=args.force)
    slot_df = materialize_arm2_slot(phage_order, matrix, args.slot_csv)

    if not args.skip_preflight:
        run_variance_preflight(
            slot_df,
            args.original_slot_csv,
            args.interactions_path,
            args.output_dir / "ch06_arm2_variance_preflight.json",
        )

    if args.run_training:
        run_arm2_training_eval(
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
