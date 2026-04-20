#!/usr/bin/env python3
"""CH06 Arm 3: Moriniere receptor-class probabilities as phage-side feature slot.

Replaces the TL17 categorical projection with a per-receptor k-mer-presence
fraction vector. For each phage P and each receptor class R, the feature value
is `|kmers(R) ∩ kmers(P)| / |kmers(R)|` — the fraction of receptor R's
characteristic k-mers that appear anywhere in phage P's proteome. Moriniere 2026
Dataset S6 provides 13 receptor classes (plan.yml says "19-dim", but the actual
Moriniere supplementary dataset has 13 — tsx, ompA/C/F, fhuA, btuB, lptD, lamB,
NGR, Kdo, HepI, HepII, GluI).

The raw 815 k-mers were null in SX12 (see `kmer-receptor-expansion-neutral`).
This arm tests whether aggregating to per-receptor normalized fractions (which
controls for class imbalance: HepI has 168 k-mers, Kdo has 9) carries
discrimination signal missed by the raw-k-mer representation. Plan.yml pre-
registers this as **expected null** per:

  (a) same-receptor-uncorrelated-hosts — Tsx phages have Jaccard 0.091 on host
      ranges; receptor class is far from sufficient for predicting lysis.
  (b) Moriniere's classifier trained on BW25113/BL21, which lack capsule /
      O-antigen, so polysaccharide-mediated specificity (Gate 1 mechanism —
      depo-capsule-validated) is not represented.

Normalizing per-receptor (vs raw count or softmax across receptors) is a small
refinement on SX12 at best.

Two-phase:

1. `build_receptor_fraction_matrix`: for each of 148 phages, compute the
   13-dim receptor-fraction vector. Uses `load_receptor_kmers` and
   `collect_*_phage_proteomes` already in the repo (SX12 infrastructure).

2. `run_variance_preflight`: CV and Cohen's d on the three standard subsets
   (Guelin, BASEL non-zero TL17, BASEL zero TL17). If the arm fails the gate,
   we skip training and document the kill.

3. `run_arm3_training_eval`: CH05 two-axis retraining via the same patch-in-place
   path as Arm 2 but with `.scratch/basel/feature_slots_arm3/` as the slot root.
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Sequence

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot import (
    DEFAULT_BASEL_PHAROKKA_DIR,
    DEFAULT_DATASET_PATH,
    DEFAULT_GUELIN_PROTEIN_DIR,
    collect_basel_phage_proteomes,
    collect_guelin_phage_proteomes,
)
from lyzortx.pipeline.autoresearch.predict_receptor_from_kmers import load_receptor_kmers

LOGGER = logging.getLogger(__name__)

DEFAULT_SLOT_CSV = Path(".scratch/basel/feature_slots_arm3/phage_projection/features.csv")
DEFAULT_ORIGINAL_SLOT_CSV = Path(".scratch/basel/feature_slots/phage_projection/features.csv")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch06_arm3_moriniere_receptor")
COLUMN_PREFIX = "phage_projection__recep_frac_"


def build_receptor_fraction_matrix(
    dataset_path: Path,
    guelin_protein_dir: Path,
    basel_pharokka_dir: Path,
) -> tuple[list[str], list[str], np.ndarray]:
    """Compute per-(phage, receptor) fraction of characteristic k-mers present.

    Returns (phage_order, receptor_order, matrix[n_phages, n_receptors]).
    """
    receptor_kmers = load_receptor_kmers(dataset_path)
    receptor_order = sorted(receptor_kmers.keys())
    guelin = collect_guelin_phage_proteomes(guelin_protein_dir)
    basel = collect_basel_phage_proteomes(basel_pharokka_dir)
    proteomes = {**guelin, **basel}
    phage_order = sorted(proteomes.keys())

    matrix = np.zeros((len(phage_order), len(receptor_order)), dtype=np.float32)
    for i, phage in enumerate(phage_order):
        proteome = proteomes[phage]
        for j, receptor in enumerate(receptor_order):
            kmers = receptor_kmers[receptor]
            hits = sum(1 for k in kmers if k in proteome)
            matrix[i, j] = hits / len(kmers) if kmers else 0.0

    LOGGER.info(
        "Receptor-fraction matrix: %d phages × %d receptors; overall mean=%.3f",
        matrix.shape[0],
        matrix.shape[1],
        float(matrix.mean()),
    )
    return phage_order, receptor_order, matrix


def materialize_arm3_slot(
    phage_order: Sequence[str],
    receptor_order: Sequence[str],
    matrix: np.ndarray,
    slot_csv: Path,
) -> pd.DataFrame:
    slot_csv.parent.mkdir(parents=True, exist_ok=True)
    columns = [f"{COLUMN_PREFIX}{r}" for r in receptor_order]
    df = pd.DataFrame(matrix, columns=columns)
    df.insert(0, "phage", list(phage_order))
    df.to_csv(slot_csv, index=False)
    LOGGER.info("Arm 3 receptor-fraction slot written to %s", slot_csv)
    return df


def run_variance_preflight(
    slot_df: pd.DataFrame,
    original_slot_csv: Path,
    interactions_path: Path,
    output_json: Path,
) -> dict[str, object]:
    """CV and Cohen's d on three phage subsets (Guelin / BASEL non-zero / BASEL zero).

    Same structure as Arm 2's pre-flight — see `ch06_arm2_mmseqs_proteome.py` for
    rationale. Gate passes if any subset has CV > 0.1 OR Cohen's d > 0.1.
    """
    original = pd.read_csv(original_slot_csv)
    feature_cols = [c for c in original.columns if c != "phage"]
    original_nonzero = original[feature_cols].abs().sum(axis=1) > 0
    nonzero_phages = set(original.loc[original_nonzero, "phage"])

    guelin = [p for p in slot_df["phage"] if not p.startswith("Bas")]
    basel_all = [p for p in slot_df["phage"] if p.startswith("Bas")]
    basel_nonzero_proj = [p for p in basel_all if p in nonzero_phages]
    basel_zero_proj = [p for p in basel_all if p not in nonzero_phages]

    feature_cols = [c for c in slot_df.columns if c.startswith(COLUMN_PREFIX)]

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
            "n_receptors": len(feature_cols),
            "max_abs_cv": round(max_abs_cv, 4),
            "max_abs_cohens_d": round(max_abs_d, 4),
            "n_receptors_with_cv_gt_0.1": int(np.nansum(np.abs(cv) > 0.1)),
            "n_receptors_with_cohens_d_gt_0.1": int(np.sum(np.abs(cohens_d) > 0.1)),
            "gate_status": "pass" if (max_abs_cv > 0.1 or max_abs_d > 0.1) else "fail",
        }

    output_json.parent.mkdir(parents=True, exist_ok=True)
    with output_json.open("w") as f:
        json.dump(summary, f, indent=2)
    LOGGER.info("Variance pre-flight: %s", json.dumps(summary, indent=2))
    return summary


def run_arm3_training_eval(
    *,
    slot_csv: Path,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
    device_type: str = "cpu",
    max_folds: int | None = None,
    num_workers: int = 3,
    drop_high_titer_only_positives: bool = True,
) -> dict[str, object]:
    """Run CH05's two-axis evaluation with the Arm 3 phage_projection slot."""
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
        raise FileNotFoundError(f"Arm 3 slot CSV missing: {slot_csv}")

    arm3_slot_dir = slot_csv.parent.parent
    arm3_slot_dir.mkdir(parents=True, exist_ok=True)
    for slot_name in ("phage_stats", "phage_rbp_struct"):
        arm3_side = arm3_slot_dir / slot_name / "features.csv"
        baseline_side = Path(".scratch/basel/feature_slots") / slot_name / "features.csv"
        if not arm3_side.exists() and baseline_side.exists():
            arm3_side.parent.mkdir(parents=True, exist_ok=True)
            arm3_side.symlink_to(baseline_side.resolve())

    unified = load_unified_row_frame(basel_log10_pfu_ml=BASEL_LOG10_PFU_ML)
    clean_rows = build_clean_row_training_frame(unified, drop_high_titer_only_positives=drop_high_titer_only_positives)
    pair_source = clean_rows[["pair_id", "source"]].drop_duplicates(subset=["pair_id"]).set_index("pair_id")["source"]
    phage_family = load_unified_phage_family_map()
    candidate_module = load_module_from_path("ch06_arm3_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)
    patch_context_with_extended_slots(context, slots_dir=arm3_slot_dir)

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
    bacteria_pairs.to_csv(output_dir / "ch06_arm3_bacteria_axis_predictions.csv", index=False)
    bacteria_cis = _bootstrap_by_unit(
        bacteria_pairs.to_dict(orient="records"),
        unit_key="bacteria",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    phage_pairs = select_pair_max_concentration_rows(phage_per_row)
    phage_pairs["source"] = phage_pairs["pair_id"].map(pair_source)
    phage_pairs.to_csv(output_dir / "ch06_arm3_phage_axis_predictions.csv", index=False)
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
    pd.DataFrame(cross_source_rows).to_csv(output_dir / "ch06_arm3_cross_source_breakdown.csv", index=False)

    summary = {
        "arm": "CH06 Arm 3 — Moriniere receptor-class probabilities",
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
    with (output_dir / "ch06_arm3_metrics.json").open("w") as f:
        json.dump(summary, f, indent=2, default=str)
    LOGGER.info("CH06 Arm 3 summary written to %s", output_dir / "ch06_arm3_metrics.json")
    return summary


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--dataset-path", type=Path, default=DEFAULT_DATASET_PATH)
    parser.add_argument("--guelin-protein-dir", type=Path, default=DEFAULT_GUELIN_PROTEIN_DIR)
    parser.add_argument("--basel-pharokka-dir", type=Path, default=DEFAULT_BASEL_PHAROKKA_DIR)
    parser.add_argument("--slot-csv", type=Path, default=DEFAULT_SLOT_CSV)
    parser.add_argument("--original-slot-csv", type=Path, default=DEFAULT_ORIGINAL_SLOT_CSV)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument(
        "--interactions-path",
        type=Path,
        default=Path("data/interactions/raw/raw_interactions.csv"),
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
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    phage_order, receptor_order, matrix = build_receptor_fraction_matrix(
        args.dataset_path, args.guelin_protein_dir, args.basel_pharokka_dir
    )
    slot_df = materialize_arm3_slot(phage_order, receptor_order, matrix, args.slot_csv)

    if not args.skip_preflight:
        run_variance_preflight(
            slot_df,
            args.original_slot_csv,
            args.interactions_path,
            args.output_dir / "ch06_arm3_variance_preflight.json",
        )

    if args.run_training:
        run_arm3_training_eval(
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
