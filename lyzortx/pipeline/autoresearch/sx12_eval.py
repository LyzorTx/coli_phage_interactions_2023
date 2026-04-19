#!/usr/bin/env python3
"""SX12: Evaluate Moriniere 5-mer phage-side features as direct model features.

Adds the `phage_moriniere_kmer` slot (815 receptor-predictive k-mers from Moriniere
Dataset S6) to the existing SX10 configuration and re-runs the SX01 + SX03 evaluation
protocols. Compares to the SX10 baseline. Acceptance gate: within-panel AUC OR SX03 Arm C
AUC improves by >=2 pp.

Usage:
    python -m lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot  # one-time slot build
    python -m lyzortx.pipeline.autoresearch.sx12_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import json
import logging
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any

import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot import FEATURE_PREFIX
from lyzortx.pipeline.autoresearch.candidate_replay import (
    build_st03_training_frame,
    load_module_from_path,
    load_st03_holdout_frame,
)
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.spandex_metrics import evaluate_holdout_rows
from lyzortx.pipeline.autoresearch.sx01_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
    bootstrap_spandex_cis,
    enrich_rows_with_mlc,
    load_mlc_scores,
    train_and_predict_fold,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx12_eval")
DEFAULT_KMER_SLOT_PATH = Path(".scratch/moriniere_kmer/features.csv")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")

MORINIERE_SLOT_NAME = "phage_moriniere_kmer"
EXTENDED_PHAGE_SLOTS = ["phage_projection", "phage_stats", MORINIERE_SLOT_NAME]


def attach_moriniere_slot(context: Any, features_path: Path) -> None:
    """Add the Moriniere k-mer slot to context.slot_artifacts in place.

    Reuses the SlotArtifact shape from `lyzortx.autoresearch.train`. Columns are the k-mer
    feature columns; entity key is `phage`.
    """
    from lyzortx.autoresearch.train import SlotArtifact

    if not features_path.exists():
        raise FileNotFoundError(
            f"Moriniere k-mer slot CSV not found at {features_path}; run "
            "`python -m lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot` first"
        )
    frame = pd.read_csv(features_path)
    feature_cols = tuple(c for c in frame.columns if c.startswith(FEATURE_PREFIX))
    if "phage" not in frame.columns or not feature_cols:
        raise ValueError(f"Moriniere slot CSV {features_path} missing phage column or k-mer features")
    artifact = SlotArtifact(
        slot_name=MORINIERE_SLOT_NAME,
        entity_key="phage",
        feature_columns=feature_cols,
        frame=frame,
    )
    # dict assignment is fine; slot_artifacts is mutable
    context.slot_artifacts[MORINIERE_SLOT_NAME] = artifact
    LOGGER.info(
        "Attached Moriniere k-mer slot: %d phages x %d k-mer features",
        len(frame),
        len(feature_cols),
    )


def run_within_panel_eval(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_frame: pd.DataFrame,
    mlc_lookup: dict[tuple[str, str], float],
    device_type: str,
    output_dir: Path,
) -> None:
    """Mirrors SX01's 10-fold CV but injects the Moriniere slot into phage features."""
    fold_assignments = assign_bacteria_folds(bacteria_to_cv_group_map(clean_frame))
    all_predictions: list[dict[str, object]] = []

    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        training_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_fold = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_fold = clean_frame[clean_frame["bacteria"].isin(training_bacteria)].copy()
        LOGGER.info(
            "=== SX12 within-panel Fold %d: %d train bacteria (%d pairs), %d holdout bacteria (%d pairs) ===",
            fold_id,
            len(training_bacteria),
            len(training_fold),
            len(holdout_bacteria),
            len(holdout_fold),
        )

        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_and_predict_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=training_fold,
                holdout_frame=holdout_fold,
                seed=seed,
                device_type=device_type,
                phage_slots=EXTENDED_PHAGE_SLOTS,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_rows.extend(rows)

        df = pd.DataFrame(fold_rows)
        agg = (
            df.groupby(["fold_id", "pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
                "predicted_probability"
            ]
            .mean()
            .sort_values(["bacteria", "phage"])
        )
        enriched = enrich_rows_with_mlc(agg.to_dict(orient="records"), mlc_lookup)
        fold_metrics = evaluate_holdout_rows(enriched)
        LOGGER.info(
            "  SX12 Fold %d: nDCG=%.4f, mAP=%.4f, AUC=%.4f, Brier=%.4f",
            fold_id,
            fold_metrics.get("holdout_ndcg") or 0,
            fold_metrics.get("holdout_map") or 0,
            fold_metrics.get("holdout_roc_auc") or 0,
            fold_metrics.get("holdout_brier_score") or 0,
        )
        all_predictions.extend(enriched)

    LOGGER.info("=== SX12 within-panel bootstrap (%d predictions) ===", len(all_predictions))
    bootstrap_results = bootstrap_spandex_cis(
        all_predictions,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )
    for metric, ci in bootstrap_results.items():
        LOGGER.info("  %s: %.4f [%.4f, %.4f]", metric, ci.point_estimate or 0, ci.ci_low or 0, ci.ci_high or 0)

    output_dir.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(all_predictions).to_csv(output_dir / "within_panel_predictions.csv", index=False)
    with open(output_dir / "within_panel_bootstrap.json", "w", encoding="utf-8") as f:
        json.dump(
            {
                metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
                for metric, ci in bootstrap_results.items()
            },
            f,
            indent=2,
        )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    parser.add_argument("--kmer-slot-path", type=Path, default=DEFAULT_KMER_SLOT_PATH)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    LOGGER.info("SX12 eval starting at %s", datetime.now(timezone.utc).isoformat())

    candidate_module = load_module_from_path("sx12_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)
    attach_moriniere_slot(context, args.kmer_slot_path)

    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    ambiguous = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean_frame), clean_frame["bacteria"].nunique())

    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    run_within_panel_eval(
        candidate_module=candidate_module,
        context=context,
        clean_frame=clean_frame,
        mlc_lookup=mlc_lookup,
        device_type=args.device_type,
        output_dir=args.output_dir,
    )

    LOGGER.info("SX12 within-panel eval complete. For SX03 Arm C re-evaluation, use sx12_eval_arm_c flag (TODO).")


if __name__ == "__main__":
    main()
