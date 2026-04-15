#!/usr/bin/env python3
"""SX14: Wave-2 consolidation + stratified evaluation.

All three wave-2 tickets (SX11 potency loss, SX12 phage Moriniere k-mers, SX13 host OMP k-mers)
produced null results on the +2 pp acceptance gate. Per the SX14 spec, arms that fail their gate
are excluded from the consolidated model, so **the consolidated wave-2 final = SX10 baseline
unchanged** (no retraining). The real deliverable of SX14 is the stratified evaluation layer
applied to existing wave-2 prediction outputs.

This evaluator loads per-arm predictions from SX10, SX11, SX12, SX13, attaches four stratum
labels per (holdout bacterium, phage) pair, and reports bootstrap CIs per stratum × arm.

Stratum definitions:
  1. within_family     — holdout phage's family has >=3 training-positive pairs on this bacterium's
                         cv_group (rich family-specific signal available)
  2. cross_family      — holdout phage's family has 0 training-positive pairs on this bacterium's
                         cv_group (family cold-start)
  3. narrow_host_phage — phage's panel-wide lysis rate <30% (the hard specialists)
  4. phylogroup_orphan — holdout bacterium has <=2 training-phylogroup-siblings in its CV fold

Outputs:
  - stratified_metrics.csv — per arm × stratum × metric with point + 95% bootstrap CI
  - all_predictions.csv    — all arms' per-pair predictions with stratum label columns attached
  - notebook_table.md      — side-by-side SX10 vs each wave-2 arm (aggregate + per-stratum)

Usage:
    python -m lyzortx.pipeline.autoresearch.sx14_eval
"""

from __future__ import annotations

import argparse
import json
import logging
from pathlib import Path
from typing import Iterable

import numpy as np
import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.sx01_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    bootstrap_spandex_cis,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx14_eval")

# Prediction sources: (arm_id, csv_path).
# SX10 canonical predictions live in sx05_sx01_eval (post-SX05 labels, per-phage blending).
# SX11 arms ran without per-phage blending; we include them for comparison but flag the difference.
# SX12 has a single arm; SX13 has four (baseline/marginal/cross_term/path1_cluster).
PREDICTION_SOURCES: dict[str, dict[str, object]] = {
    "sx10_baseline": {
        "path": Path("lyzortx/generated_outputs/sx05_sx01_eval/kfold_predictions.csv"),
        "note": "Wave-1 canonical baseline (post-SX05 labels, per-phage blending)",
    },
    "sx11_binary_baseline": {
        "path": Path("lyzortx/generated_outputs/sx11_eval/binary_baseline_predictions.csv"),
        "note": "SX11 baseline arm (binary, NO per-phage blending — loss-function ablation reference)",
    },
    "sx11_hurdle": {
        "path": Path("lyzortx/generated_outputs/sx11_eval/hurdle_two_stage_predictions.csv"),
        "note": "SX11 hurdle two-stage arm",
    },
    "sx11_lambdarank": {
        "path": Path("lyzortx/generated_outputs/sx11_eval/lambdarank_predictions.csv"),
        "note": "SX11 LambdaRank arm",
    },
    "sx11_ordinal": {
        "path": Path("lyzortx/generated_outputs/sx11_eval/ordinal_all_threshold_predictions.csv"),
        "note": "SX11 ordinal all-threshold arm (best SX11 performer)",
    },
    "sx12_moriniere_kmer": {
        "path": Path("lyzortx/generated_outputs/sx12_eval/within_panel_predictions.csv"),
        "note": "SX12 phage Moriniere k-mers (with per-phage blending)",
    },
    "sx13_baseline": {
        "path": Path("lyzortx/generated_outputs/sx13_eval/baseline_predictions.csv"),
        "note": "SX13 baseline (mirrors SX10; reruns for internal consistency)",
    },
    "sx13_marginal": {
        "path": Path("lyzortx/generated_outputs/sx13_eval/marginal_predictions.csv"),
        "note": "SX13 host_omp_kmer marginal arm",
    },
    "sx13_cross_term": {
        "path": Path("lyzortx/generated_outputs/sx13_eval/cross_term_predictions.csv"),
        "note": "SX13 phage_kmer × host_omp_kmer cross-term arm",
    },
    "sx13_path1_cluster": {
        "path": Path("lyzortx/generated_outputs/sx13_eval/path1_cluster_predictions.csv"),
        "note": "SX13 MMseqs2 cluster-ID fallback arm",
    },
}

# Metadata sources (read once).
GUELIN_PHAGE_CSV = Path("data/genomics/phages/guelin_collection.csv")
HOST_FEATURES_CSV = Path("data/genomics/bacteria/isolation_strains/panacota/tree/itol/370+host_features.csv")
HOST_CV_GROUP_CSV = Path("data/metadata/370+host_cross_validation_groups_1e-4.csv")
INTERACTION_MATRIX_CSV = Path("data/interactions/interaction_matrix.csv")

NARROW_HOST_LYSIS_THRESHOLD = 0.30
WITHIN_FAMILY_MIN_POSITIVES = 3
PHYLOGROUP_ORPHAN_MAX_SIBLINGS = 2


def load_phage_family_map() -> dict[str, str]:
    """Return phage -> family mapping from the Guelin collection CSV (BASEL phages not present)."""
    df = pd.read_csv(GUELIN_PHAGE_CSV, sep=";")
    if "phage" not in df.columns or "Family" not in df.columns:
        raise ValueError(f"{GUELIN_PHAGE_CSV}: expected 'phage' and 'Family' columns, got {list(df.columns)}")
    mapping = dict(zip(df["phage"].astype(str), df["Family"].astype(str)))
    LOGGER.info("Loaded %d phage family assignments from %s", len(mapping), GUELIN_PHAGE_CSV)
    return mapping


def load_host_phylogroup_map() -> dict[str, str]:
    """Return bacteria -> phylogroup mapping from the PanACoTA host-features CSV."""
    df = pd.read_csv(HOST_FEATURES_CSV, sep="\t")
    if "bacteria" not in df.columns or "Clermont_Phylo" not in df.columns:
        raise ValueError(
            f"{HOST_FEATURES_CSV}: expected 'bacteria' and 'Clermont_Phylo' columns, got {list(df.columns)[:15]}"
        )
    mapping = dict(zip(df["bacteria"].astype(str), df["Clermont_Phylo"].astype(str)))
    LOGGER.info("Loaded %d host phylogroup assignments from %s", len(mapping), HOST_FEATURES_CSV)
    return mapping


def load_host_cv_group_map() -> dict[str, int]:
    df = pd.read_csv(HOST_CV_GROUP_CSV, sep=";")
    if "bacteria" not in df.columns or "group" not in df.columns:
        raise ValueError(f"{HOST_CV_GROUP_CSV}: expected 'bacteria' and 'group' columns, got {list(df.columns)}")
    mapping = {str(b): int(g) for b, g in zip(df["bacteria"], df["group"])}
    LOGGER.info("Loaded %d host cv_group assignments from %s", len(mapping), HOST_CV_GROUP_CSV)
    return mapping


def load_interaction_matrix() -> pd.DataFrame:
    """Return the long-format MLC matrix: rows = (bacteria, phage, mlc_score)."""
    wide = pd.read_csv(INTERACTION_MATRIX_CSV, sep=";", index_col=0)
    long = wide.reset_index().melt(id_vars=wide.index.name or "bacteria", var_name="phage", value_name="mlc_score")
    long = long.rename(columns={long.columns[0]: "bacteria"})
    long["mlc_score"] = pd.to_numeric(long["mlc_score"], errors="coerce").fillna(0)
    return long


def compute_phage_lysis_rates(long_matrix: pd.DataFrame) -> dict[str, float]:
    """Panel-wide lysis rate per phage: fraction of hosts with MLC >= 1."""
    any_lysis = (long_matrix["mlc_score"] >= 1).astype(int)
    rates = long_matrix.assign(any_lysis=any_lysis).groupby("phage")["any_lysis"].mean()
    return {str(k): float(v) for k, v in rates.items()}


def attach_stratum_labels(
    predictions: pd.DataFrame,
    *,
    phage_family: dict[str, str],
    host_phylogroup: dict[str, str],
    host_cv_group: dict[str, int],
    phage_lysis: dict[str, float],
    long_matrix: pd.DataFrame,
) -> pd.DataFrame:
    """Return predictions df with 4 boolean stratum columns attached.

    Requires predictions to have: bacteria, phage, fold_id. The fold_id is used for
    phylogroup-orphan computation (training = pairs NOT in this fold).
    """
    required = {"bacteria", "phage", "fold_id"}
    missing = required - set(predictions.columns)
    if missing:
        raise ValueError(f"predictions missing columns: {missing}")

    out = predictions.copy()
    out["phage_family"] = out["phage"].map(phage_family).fillna("UNKNOWN")
    out["host_phylogroup"] = out["bacteria"].map(host_phylogroup).fillna("UNKNOWN")
    out["host_cv_group"] = out["bacteria"].map(host_cv_group).fillna(-1).astype(int)
    out["phage_lysis_rate"] = out["phage"].map(phage_lysis).fillna(0.0)

    # Narrow-host phage: direct threshold
    out["stratum_narrow_host_phage"] = out["phage_lysis_rate"] < NARROW_HOST_LYSIS_THRESHOLD

    # Within/cross-family: for each holdout (bacterium, phage) pair, count training positives
    # where the phage family matches AND the training bacterium shares the holdout's cv_group.
    # Training set per fold = bacteria NOT in this fold_id.
    long_matrix = long_matrix.copy()
    long_matrix["bacteria"] = long_matrix["bacteria"].astype(str)
    long_matrix["phage"] = long_matrix["phage"].astype(str)
    long_matrix["any_lysis"] = (long_matrix["mlc_score"] >= 1).astype(int)
    long_matrix["phage_family"] = long_matrix["phage"].map(phage_family).fillna("UNKNOWN")
    long_matrix["host_cv_group"] = long_matrix["bacteria"].map(host_cv_group).fillna(-1).astype(int)

    fold_to_bacteria: dict[int, set[str]] = (
        out.groupby("fold_id")["bacteria"].agg(lambda x: set(x.astype(str))).to_dict()
    )

    # Index: (phage_family, host_cv_group) -> count of training positives per fold
    family_cvg_positive_counts: dict[int, dict[tuple[str, int], int]] = {}
    for fold_id, holdout_bacteria in fold_to_bacteria.items():
        train_mask = ~long_matrix["bacteria"].isin(holdout_bacteria) & (long_matrix["any_lysis"] == 1)
        train_long = long_matrix.loc[train_mask]
        counts = train_long.groupby(["phage_family", "host_cv_group"]).size().to_dict()
        family_cvg_positive_counts[fold_id] = counts

    def within_cross_family(row: pd.Series) -> tuple[bool, bool]:
        fold_id = int(row["fold_id"])
        key = (row["phage_family"], int(row["host_cv_group"]))
        count = family_cvg_positive_counts.get(fold_id, {}).get(key, 0)
        return count >= WITHIN_FAMILY_MIN_POSITIVES, count == 0

    wc = out.apply(within_cross_family, axis=1, result_type="expand")
    wc.columns = ["stratum_within_family", "stratum_cross_family"]
    out = pd.concat([out, wc], axis=1)

    # Phylogroup-orphan: for each holdout bacterium, count training-set bacteria sharing phylogroup
    out["__train_pgrp_siblings"] = 0
    bacteria_to_pgrp = {b: host_phylogroup.get(b, "UNKNOWN") for b in out["bacteria"].astype(str).unique()}
    for fold_id, holdout_bacteria in fold_to_bacteria.items():
        # Training bacteria = all bacteria minus holdout bacteria for this fold.
        all_bacteria = set(out["bacteria"].astype(str).unique())
        train_bacteria = all_bacteria - holdout_bacteria
        pgrp_counts: dict[str, int] = {}
        for b in train_bacteria:
            pgrp = bacteria_to_pgrp.get(b, "UNKNOWN")
            pgrp_counts[pgrp] = pgrp_counts.get(pgrp, 0) + 1
        fold_mask = out["fold_id"] == fold_id
        out.loc[fold_mask, "__train_pgrp_siblings"] = (
            out.loc[fold_mask, "host_phylogroup"].map(pgrp_counts).fillna(0).astype(int)
        )
    out["stratum_phylogroup_orphan"] = out["__train_pgrp_siblings"] <= PHYLOGROUP_ORPHAN_MAX_SIBLINGS
    out = out.drop(columns=["__train_pgrp_siblings"])

    return out


def prepare_rows_for_bootstrap(df: pd.DataFrame) -> list[dict[str, object]]:
    """The bootstrap function expects dict rows with predicted_probability and label_binary/mlc_score."""
    required = {"bacteria", "phage", "predicted_probability", "mlc_score"}
    missing = required - set(df.columns)
    if missing:
        raise ValueError(f"bootstrap input missing: {missing}")
    rows = []
    for _, r in df.iterrows():
        rows.append(
            {
                "bacteria": str(r["bacteria"]),
                "phage": str(r["phage"]),
                "predicted_probability": float(r["predicted_probability"]),
                "mlc_score": float(r["mlc_score"]) if pd.notna(r["mlc_score"]) else 0.0,
                "label_binary": int((r["mlc_score"] or 0) > 0) if pd.notna(r["mlc_score"]) else 0,
            }
        )
    return rows


def bootstrap_per_stratum(
    arm_predictions: pd.DataFrame,
    arm_id: str,
    stratum_cols: Iterable[str],
) -> list[dict[str, object]]:
    """For each stratum (aggregate + one per stratum flag TRUE), compute bootstrap CIs."""
    results: list[dict[str, object]] = []

    def record(stratum: str, subset: pd.DataFrame) -> None:
        if len(subset) < 10 or subset["mlc_score"].fillna(0).max() == 0:
            LOGGER.info(
                "  arm=%s stratum=%s: skipped (n=%d, too few rows or no positives)", arm_id, stratum, len(subset)
            )
            return
        rows = prepare_rows_for_bootstrap(subset)
        cis = bootstrap_spandex_cis(
            rows, bootstrap_samples=BOOTSTRAP_SAMPLES, bootstrap_random_state=BOOTSTRAP_RANDOM_STATE
        )
        for metric, ci in cis.items():
            results.append(
                {
                    "arm_id": arm_id,
                    "stratum": stratum,
                    "n_pairs": len(subset),
                    "n_positive": int((subset["mlc_score"] > 0).sum()),
                    "metric": metric,
                    "point_estimate": ci.point_estimate,
                    "ci_low": ci.ci_low,
                    "ci_high": ci.ci_high,
                }
            )

    record("aggregate", arm_predictions)
    for col in stratum_cols:
        label = col.replace("stratum_", "")
        record(label, arm_predictions[arm_predictions[col]])
    return results


def build_arm_table(
    arm_ids: Iterable[str],
    stratum_cols: Iterable[str],
    source_rows: dict[str, pd.DataFrame],
    metadata: dict[str, object],
) -> tuple[pd.DataFrame, pd.DataFrame]:
    """Returns (stratified_metrics_df, all_predictions_df).

    all_predictions_df has a `source_arm` column + all stratum columns, for reproducibility.
    """
    all_preds: list[pd.DataFrame] = []
    metric_rows: list[dict[str, object]] = []
    for arm_id in arm_ids:
        if arm_id not in source_rows:
            LOGGER.warning("Arm %s not loaded — skipping", arm_id)
            continue
        df = source_rows[arm_id]
        enriched = attach_stratum_labels(df, **metadata)
        enriched.insert(0, "source_arm", arm_id)
        all_preds.append(enriched)
        LOGGER.info(
            "arm=%s: n=%d pairs, within_family=%d, cross_family=%d, narrow_host_phage=%d, phylogroup_orphan=%d",
            arm_id,
            len(enriched),
            enriched["stratum_within_family"].sum(),
            enriched["stratum_cross_family"].sum(),
            enriched["stratum_narrow_host_phage"].sum(),
            enriched["stratum_phylogroup_orphan"].sum(),
        )
        metric_rows.extend(bootstrap_per_stratum(enriched, arm_id, stratum_cols))
    return pd.DataFrame(metric_rows), pd.concat(all_preds, ignore_index=True)


def pivot_for_notebook(metrics: pd.DataFrame) -> str:
    """Render a side-by-side markdown table per metric per arm across strata.

    Manually formats a GFM pipe table (no tabulate dep).
    """
    out_lines = []
    metrics = metrics.copy()
    metrics["value"] = metrics.apply(
        lambda r: f"{r['point_estimate']:.4f} [{r['ci_low']:.3f}, {r['ci_high']:.3f}]"
        if pd.notna(r["ci_low"])
        else f"{r['point_estimate']:.4f}",
        axis=1,
    )
    col_order = ["aggregate", "within_family", "cross_family", "narrow_host_phage", "phylogroup_orphan"]
    for metric_name in ["holdout_ndcg", "holdout_map", "holdout_roc_auc", "holdout_brier_score"]:
        sub = metrics[metrics["metric"] == metric_name]
        if sub.empty:
            continue
        pivot = sub.pivot_table(index="arm_id", columns="stratum", values="value", aggfunc="first")
        cols = [c for c in col_order if c in pivot.columns]
        pivot = pivot[cols]
        out_lines.append(f"\n### {metric_name}\n")
        header = "| arm | " + " | ".join(cols) + " |"
        sep = "|" + "|".join(["---"] * (len(cols) + 1)) + "|"
        out_lines.append(header)
        out_lines.append(sep)
        for arm in pivot.index:
            row_vals = [str(pivot.loc[arm, c]) if pd.notna(pivot.loc[arm, c]) else "—" for c in cols]
            out_lines.append(f"| {arm} | " + " | ".join(row_vals) + " |")
    return "\n".join(out_lines)


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--arms", type=str, default=",".join(PREDICTION_SOURCES.keys()))
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    args.output_dir.mkdir(parents=True, exist_ok=True)

    arm_ids = [a.strip() for a in args.arms.split(",") if a.strip()]
    unknown = set(arm_ids) - set(PREDICTION_SOURCES)
    if unknown:
        raise ValueError(f"Unknown arm ids: {sorted(unknown)}")

    # Load metadata once
    phage_family = load_phage_family_map()
    host_phylogroup = load_host_phylogroup_map()
    host_cv_group = load_host_cv_group_map()
    long_matrix = load_interaction_matrix()
    phage_lysis = compute_phage_lysis_rates(long_matrix)
    LOGGER.info(
        "Phage lysis rates: min=%.3f, median=%.3f, max=%.3f, narrow(<30%%)=%d/%d",
        min(phage_lysis.values()),
        float(np.median(list(phage_lysis.values()))),
        max(phage_lysis.values()),
        sum(1 for v in phage_lysis.values() if v < NARROW_HOST_LYSIS_THRESHOLD),
        len(phage_lysis),
    )

    metadata = {
        "phage_family": phage_family,
        "host_phylogroup": host_phylogroup,
        "host_cv_group": host_cv_group,
        "phage_lysis": phage_lysis,
        "long_matrix": long_matrix,
    }

    stratum_cols = [
        "stratum_within_family",
        "stratum_cross_family",
        "stratum_narrow_host_phage",
        "stratum_phylogroup_orphan",
    ]

    source_rows: dict[str, pd.DataFrame] = {}
    for arm_id in arm_ids:
        src = PREDICTION_SOURCES[arm_id]
        path = src["path"]
        if not path.exists():
            LOGGER.warning("Missing prediction file for %s: %s — skipping", arm_id, path)
            continue
        df = pd.read_csv(path)
        # Normalize column names across sources (SX01 uses predicted_probability directly; some use mlc_score as a renamed field)
        if "mlc_score" not in df.columns:
            df["mlc_score"] = 0.0
        if "fold_id" not in df.columns:
            raise ValueError(f"{path}: missing fold_id column")
        source_rows[arm_id] = df
        LOGGER.info("Loaded %s: %d rows from %s", arm_id, len(df), path)

    metrics_df, all_preds_df = build_arm_table(arm_ids, stratum_cols, source_rows, metadata)
    metrics_csv = args.output_dir / "stratified_metrics.csv"
    all_preds_csv = args.output_dir / "all_predictions.csv"
    metrics_df.to_csv(metrics_csv, index=False)
    all_preds_df.to_csv(all_preds_csv, index=False)
    LOGGER.info("Wrote %d metric rows to %s", len(metrics_df), metrics_csv)
    LOGGER.info("Wrote %d prediction rows to %s", len(all_preds_df), all_preds_csv)

    md_table = pivot_for_notebook(metrics_df)
    notebook_md = args.output_dir / "notebook_table.md"
    notebook_md.write_text(md_table, encoding="utf-8")
    LOGGER.info("Wrote side-by-side table to %s", notebook_md)

    # Also dump a JSON summary for machine consumption
    summary = {
        "arms": arm_ids,
        "strata": ["aggregate"] + [c.replace("stratum_", "") for c in stratum_cols],
        "metrics_csv": str(metrics_csv),
        "predictions_csv": str(all_preds_csv),
    }
    (args.output_dir / "summary.json").write_text(json.dumps(summary, indent=2), encoding="utf-8")


if __name__ == "__main__":
    main()
