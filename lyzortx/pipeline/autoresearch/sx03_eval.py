#!/usr/bin/env python3
"""SX03: BASEL data integration + cross-source evaluation.

Compares three training configurations:
  Arm A: Our clean data only (SX01 baseline replication)
  Arm B: Our clean data + BASEL training pairs
  Arm C: Generalization — train on all our data, predict BASEL × ECOR pairs

Pre-flight: check BASEL phage feature-space overlap with Guelin.

Usage:
    python -m lyzortx.pipeline.autoresearch.sx03_eval --device-type cpu
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
from lyzortx.pipeline.autoresearch.candidate_replay import (
    build_st03_training_frame,
    load_module_from_path,
    load_st03_holdout_frame,
)
from lyzortx.pipeline.autoresearch.gt09_clean_label_eval import identify_ambiguous_pairs
from lyzortx.pipeline.autoresearch.spandex_metrics import (
    evaluate_holdout_rows,
)
from lyzortx.pipeline.autoresearch.sx01_eval import (
    BOOTSTRAP_RANDOM_STATE,
    BOOTSTRAP_SAMPLES,
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
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/sx03_eval")
RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")
BASEL_INTERACTIONS_PATH = Path(".scratch/genophi_data/BASEL_ECOR_interaction_matrix.csv")
EXTENDED_SLOTS_DIR = Path(".scratch/basel/feature_slots")

N_FOLDS = 10


def patch_context_with_extended_slots(context: Any) -> None:
    """Replace phage slot artifacts with extended versions (Guelin + BASEL)."""
    for slot_name in ("phage_projection", "phage_stats", "phage_rbp_struct"):
        extended_csv = EXTENDED_SLOTS_DIR / slot_name / "features.csv"
        if not extended_csv.exists():
            LOGGER.warning("Extended slot %s not found at %s — using original", slot_name, extended_csv)
            continue
        if slot_name not in context.slot_artifacts:
            continue
        art = context.slot_artifacts[slot_name]
        extended_df = pd.read_csv(extended_csv)
        original_count = len(art.frame)
        object.__setattr__(art, "frame", extended_df)
        LOGGER.info(
            "Patched %s: %d → %d phages",
            slot_name,
            original_count,
            len(extended_df),
        )


def load_basel_interactions() -> pd.DataFrame:
    """Load BASEL interaction data and harmonize with our pair format."""
    basel = pd.read_csv(BASEL_INTERACTIONS_PATH)
    # Harmonize bacteria names: ECOR12 → ECOR-12
    basel["bacteria"] = basel["strain"].apply(lambda s: f"ECOR-{s[4:]}" if s.startswith("ECOR") else s)
    basel["label_any_lysis"] = basel["interaction"].astype(str)
    basel["pair_id"] = basel["bacteria"] + "__" + basel["phage"]
    basel["source"] = "basel"
    LOGGER.info(
        "BASEL interactions: %d pairs (%d positive), %d phages × %d bacteria",
        len(basel),
        (basel["interaction"] == 1).sum(),
        basel["phage"].nunique(),
        basel["bacteria"].nunique(),
    )
    return basel


def run_preflight_overlap(
    candidate_module: ModuleType,
    context: Any,
    output_dir: Path,
) -> dict[str, object]:
    """Pre-flight: check BASEL phage feature-space overlap with Guelin."""
    from sklearn.neighbors import NearestNeighbors

    phage_slots = ["phage_projection", "phage_stats"]
    phage_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=phage_slots, entity_key="phage"
    )
    phage_typed, _, _ = candidate_module.type_entity_features(phage_table, "phage")

    # Separate Guelin and BASEL phages.
    guelin_mask = ~phage_typed["phage"].str.startswith("Bas")
    guelin_features = phage_typed[guelin_mask].drop(columns=["phage"]).fillna(0)
    basel_features = phage_typed[~guelin_mask].drop(columns=["phage"]).fillna(0)

    if len(basel_features) == 0:
        LOGGER.warning("No BASEL phages found in feature table — pre-flight skipped")
        return {"status": "skipped", "reason": "no BASEL phages in feature table"}

    # Fit NN on Guelin, query BASEL.
    nn = NearestNeighbors(n_neighbors=1, metric="cosine")
    nn.fit(guelin_features.values)
    distances, _ = nn.kneighbors(basel_features.values)
    mean_dist = float(distances.mean())
    overlap_fraction = float((distances < 0.1).mean())  # fraction within cosine 0.1

    result = {
        "status": "pass" if overlap_fraction < 0.9 else "flag",
        "mean_cosine_distance": round(mean_dist, 4),
        "overlap_fraction_below_0.1": round(overlap_fraction, 4),
        "n_guelin": len(guelin_features),
        "n_basel": len(basel_features),
    }
    LOGGER.info(
        "Pre-flight overlap: mean cosine distance=%.4f, overlap fraction=%.2f%% (threshold 90%%)",
        mean_dist,
        overlap_fraction * 100,
    )

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "preflight_overlap.json", "w", encoding="utf-8") as f:
        json.dump(result, f, indent=2)
    return result


def run_arm_a_baseline(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_frame: pd.DataFrame,
    mlc_lookup: dict[tuple[str, str], float],
    device_type: str,
) -> list[dict[str, object]]:
    """Arm A: k-fold CV on our clean data only (SX01 replication)."""
    LOGGER.info("=== Arm A: Our clean data only (SX01 baseline replication) ===")
    fold_assignments = assign_bacteria_folds(bacteria_to_cv_group_map(clean_frame))

    all_predictions: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_frame = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_frame = clean_frame[clean_frame["bacteria"].isin(train_bacteria)].copy()

        LOGGER.info("Arm A Fold %d: %d train, %d holdout", fold_id, len(training_frame), len(holdout_frame))
        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_and_predict_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=training_frame,
                holdout_frame=holdout_frame,
                seed=seed,
                device_type=device_type,
            )
            fold_rows.extend(rows)

        # Aggregate seeds.
        df = pd.DataFrame(fold_rows)
        agg = df.groupby(["pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
            "predicted_probability"
        ].mean()
        enriched = enrich_rows_with_mlc(agg.to_dict(orient="records"), mlc_lookup)
        for r in enriched:
            r["arm_id"] = "arm_a_baseline"
        all_predictions.extend(enriched)

    metrics = evaluate_holdout_rows(all_predictions)
    LOGGER.info(
        "Arm A: nDCG=%.4f, mAP=%.4f, AUC=%.4f",
        metrics.get("holdout_ndcg") or 0,
        metrics.get("holdout_map") or 0,
        metrics.get("holdout_roc_auc") or 0,
    )
    return all_predictions


def run_arm_b_pooled(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_frame: pd.DataFrame,
    basel_frame: pd.DataFrame,
    mlc_lookup: dict[tuple[str, str], float],
    device_type: str,
    deposcope_dir: Path | None = None,
) -> list[dict[str, object]]:
    """Arm B: k-fold CV with BASEL added to training (same folds as Arm A)."""
    LOGGER.info("=== Arm B: Our clean data + BASEL training ===")
    all_bacteria = sorted(clean_frame["bacteria"].unique())
    fold_assignments = assign_bacteria_folds(bacteria_to_cv_group_map(clean_frame))

    # BASEL bacteria that are also in our panel.
    basel_bacteria = set(basel_frame["bacteria"].unique())
    shared_bacteria = basel_bacteria & set(all_bacteria)
    LOGGER.info("Shared bacteria (our panel ∩ BASEL): %d", len(shared_bacteria))

    all_predictions: list[dict[str, object]] = []
    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        holdout_frame = clean_frame[clean_frame["bacteria"].isin(holdout_bacteria)].copy()
        training_frame = clean_frame[clean_frame["bacteria"].isin(train_bacteria)].copy()

        # Add BASEL pairs where bacteria is in training fold (not holdout).
        basel_train = basel_frame[
            basel_frame["bacteria"].isin(train_bacteria) & basel_frame["bacteria"].isin(shared_bacteria)
        ].copy()
        if len(basel_train) > 0:
            pooled_training = pd.concat([training_frame, basel_train], ignore_index=True)
        else:
            pooled_training = training_frame

        LOGGER.info(
            "Arm B Fold %d: %d train (%d BASEL), %d holdout",
            fold_id,
            len(pooled_training),
            len(basel_train),
            len(holdout_frame),
        )

        fold_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows = train_and_predict_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=pooled_training,
                holdout_frame=holdout_frame,
                seed=seed,
                device_type=device_type,
                deposcope_dir=deposcope_dir,
            )
            fold_rows.extend(rows)

        df = pd.DataFrame(fold_rows)
        agg = df.groupby(["pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
            "predicted_probability"
        ].mean()
        enriched = enrich_rows_with_mlc(agg.to_dict(orient="records"), mlc_lookup)
        for r in enriched:
            r["arm_id"] = "arm_b_pooled"
        all_predictions.extend(enriched)

    metrics = evaluate_holdout_rows(all_predictions)
    LOGGER.info(
        "Arm B: nDCG=%.4f, mAP=%.4f, AUC=%.4f",
        metrics.get("holdout_ndcg") or 0,
        metrics.get("holdout_map") or 0,
        metrics.get("holdout_roc_auc") or 0,
    )
    return all_predictions


def run_arm_c_generalization(
    *,
    candidate_module: ModuleType,
    context: Any,
    clean_frame: pd.DataFrame,
    basel_frame: pd.DataFrame,
    device_type: str,
    deposcope_dir: Path | None = None,
) -> list[dict[str, object]]:
    """Arm C: Train on all our data, predict BASEL × ECOR pairs."""
    LOGGER.info("=== Arm C: Generalization — predict BASEL phage interactions ===")

    # Train on ALL our clean data (no holdout).
    training_frame = clean_frame.copy()

    # Build holdout from BASEL pairs with shared bacteria.
    shared_bacteria = set(clean_frame["bacteria"].unique()) & set(basel_frame["bacteria"].unique())
    holdout_frame = basel_frame[basel_frame["bacteria"].isin(shared_bacteria)].copy()

    LOGGER.info(
        "Arm C: %d train pairs, %d BASEL holdout pairs (%d bacteria)",
        len(training_frame),
        len(holdout_frame),
        len(shared_bacteria),
    )

    all_rows: list[dict[str, object]] = []
    for seed in SEEDS:
        rows = train_and_predict_fold(
            candidate_module=candidate_module,
            context=context,
            training_frame=training_frame,
            holdout_frame=holdout_frame,
            seed=seed,
            device_type=device_type,
            deposcope_dir=deposcope_dir,
        )
        all_rows.extend(rows)

    # Aggregate seeds.
    df = pd.DataFrame(all_rows)
    agg = df.groupby(["pair_id", "bacteria", "phage", "label_hard_any_lysis"], as_index=False)[
        "predicted_probability"
    ].mean()

    # Enrich with BASEL-specific relevance (binary: 0 or 1).
    enriched = []
    for row in agg.to_dict(orient="records"):
        label = int(row["label_hard_any_lysis"])
        enriched.append(
            {
                **row,
                "arm_id": "arm_c_generalization",
                "mlc_score": float(label),  # BASEL: binary 0/1 mapped to relevance
                "label_binary": label,
            }
        )

    metrics = evaluate_holdout_rows(enriched)
    LOGGER.info(
        "Arm C: nDCG=%.4f, mAP=%.4f, AUC=%.4f",
        metrics.get("holdout_ndcg") or 0,
        metrics.get("holdout_map") or 0,
        metrics.get("holdout_roc_auc") or 0,
    )
    return enriched


def run_ranking_displacement_diagnostic(
    arm_a: list[dict[str, object]],
    arm_b: list[dict[str, object]],
    arm_c: list[dict[str, object]],
    mlc_lookup: dict[tuple[str, str], float],
    output_dir: Path,
) -> None:
    """Diagnostic: do BASEL phages steal top-3 ranking slots from Guelin phages?

    For ECOR holdout bacteria with both Guelin and BASEL ground truth, compare
    the top-ranked phages in Arm A (96 Guelin only) vs Arm B (96 Guelin + BASEL training,
    evaluated on same Guelin holdout pairs).
    """
    LOGGER.info("=== Diagnostic: ranking displacement analysis ===")

    # Get ECOR bacteria from Arm C (these have BASEL ground truth).
    arm_c_df = pd.DataFrame(arm_c)
    ecor_bacteria = set(arm_c_df["bacteria"].unique())

    # Filter Arm A and B predictions to ECOR bacteria only.
    arm_a_ecor = [r for r in arm_a if str(r["bacteria"]) in ecor_bacteria]
    arm_b_ecor = [r for r in arm_b if str(r["bacteria"]) in ecor_bacteria]

    if not arm_a_ecor or not arm_b_ecor:
        LOGGER.warning("No ECOR bacteria found in holdout — skipping displacement diagnostic")
        return

    LOGGER.info("ECOR bacteria in holdout: %d", len(ecor_bacteria))

    # Compare top-3 phages per bacterium between Arm A and Arm B.
    displacement_count = 0
    total_bacteria = 0
    details = []

    for bacteria in sorted(ecor_bacteria):
        a_rows = sorted(
            [r for r in arm_a_ecor if str(r["bacteria"]) == bacteria],
            key=lambda r: -float(r["predicted_probability"]),
        )
        b_rows = sorted(
            [r for r in arm_b_ecor if str(r["bacteria"]) == bacteria],
            key=lambda r: -float(r["predicted_probability"]),
        )

        if len(a_rows) < 3 or len(b_rows) < 3:
            continue

        a_top3 = {str(r["phage"]) for r in a_rows[:3]}
        b_top3 = {str(r["phage"]) for r in b_rows[:3]}
        displaced = a_top3 - b_top3
        total_bacteria += 1

        if displaced:
            displacement_count += 1
            details.append(
                {
                    "bacteria": bacteria,
                    "a_top3": sorted(a_top3),
                    "b_top3": sorted(b_top3),
                    "displaced": sorted(displaced),
                }
            )

    LOGGER.info(
        "Displacement: %d/%d ECOR bacteria had different top-3 phages between Arm A and B",
        displacement_count,
        total_bacteria,
    )
    for d in details:
        LOGGER.info("  %s: displaced %s", d["bacteria"], d["displaced"])

    output_dir.mkdir(parents=True, exist_ok=True)
    with open(output_dir / "displacement_diagnostic.json", "w", encoding="utf-8") as f:
        json.dump(
            {"displacement_count": displacement_count, "total_bacteria": total_bacteria, "details": details},
            f,
            indent=2,
        )


def run_sx03_eval(
    *,
    candidate_module: ModuleType,
    context: Any,
    device_type: str,
    output_dir: Path,
) -> None:
    start_time = datetime.now(timezone.utc)

    # Load data.
    holdout_frame = load_st03_holdout_frame()
    training_frame = build_st03_training_frame()
    full_frame = pd.concat([training_frame, holdout_frame], ignore_index=True)
    ambiguous_pairs = identify_ambiguous_pairs(RAW_INTERACTIONS_PATH)
    clean_frame = full_frame[~full_frame["pair_id"].isin(ambiguous_pairs)].copy()
    LOGGER.info("Clean frame: %d pairs, %d bacteria", len(clean_frame), clean_frame["bacteria"].nunique())

    # Load BASEL interactions.
    basel_interactions = load_basel_interactions()

    # Harmonize BASEL with our frame format.
    # BASEL needs the same columns as our training frame for train_and_predict_fold.
    # Minimum: pair_id, bacteria, phage, label_any_lysis, training_weight_v3
    basel_frame = basel_interactions[["pair_id", "bacteria", "phage", "label_any_lysis", "source"]].copy()
    basel_frame["training_weight_v3"] = "1.0"

    # Patch cache with extended phage slots (Guelin + BASEL).
    patch_context_with_extended_slots(context)

    # Load MLC scores for graded evaluation.
    mlc_df = load_mlc_scores()
    mlc_lookup = {(r["bacteria"], r["phage"]): r["mlc_score"] for _, r in mlc_df.iterrows()}

    # Combined DepoScope dir for BASEL phage cross-terms.
    combined_depo_dir = Path(".scratch/deposcope_combined")
    if not combined_depo_dir.exists():
        LOGGER.warning("Combined DepoScope dir not found — BASEL phage depo features will be missing from cross-terms")
        combined_depo_dir = None

    # Pre-flight overlap check.
    run_preflight_overlap(candidate_module, context, output_dir)

    # Run arms.
    arm_a = run_arm_a_baseline(
        candidate_module=candidate_module,
        context=context,
        clean_frame=clean_frame,
        mlc_lookup=mlc_lookup,
        device_type=device_type,
    )

    arm_b = run_arm_b_pooled(
        candidate_module=candidate_module,
        context=context,
        clean_frame=clean_frame,
        basel_frame=basel_frame,
        mlc_lookup=mlc_lookup,
        device_type=device_type,
        deposcope_dir=combined_depo_dir,
    )

    arm_c = run_arm_c_generalization(
        candidate_module=candidate_module,
        context=context,
        clean_frame=clean_frame,
        basel_frame=basel_frame,
        device_type=device_type,
        deposcope_dir=combined_depo_dir,
    )

    # Diagnostic: do BASEL phages steal top ranking slots for ECOR holdout bacteria?
    run_ranking_displacement_diagnostic(arm_a, arm_b, arm_c, mlc_lookup, output_dir)

    # Bootstrap CIs per arm.
    LOGGER.info("Computing bootstrap CIs...")
    results = {}
    for arm_id, predictions in [("arm_a_baseline", arm_a), ("arm_b_pooled", arm_b), ("arm_c_generalization", arm_c)]:
        ci = bootstrap_spandex_cis(
            predictions, bootstrap_samples=BOOTSTRAP_SAMPLES, bootstrap_random_state=BOOTSTRAP_RANDOM_STATE
        )
        results[arm_id] = ci
        LOGGER.info("  %s:", arm_id)
        for metric, val in ci.items():
            LOGGER.info("    %s: %.4f [%.4f, %.4f]", metric, val.point_estimate or 0, val.ci_low or 0, val.ci_high or 0)

    # Write outputs.
    output_dir.mkdir(parents=True, exist_ok=True)
    all_preds = arm_a + arm_b + arm_c
    pd.DataFrame(all_preds).to_csv(output_dir / "all_predictions.csv", index=False)

    bootstrap_json = {}
    for arm_id, ci_dict in results.items():
        bootstrap_json[arm_id] = {
            metric: {"point_estimate": ci.point_estimate, "ci_low": ci.ci_low, "ci_high": ci.ci_high}
            for metric, ci in ci_dict.items()
        }
    with open(output_dir / "bootstrap_results.json", "w", encoding="utf-8") as f:
        json.dump(bootstrap_json, f, indent=2)

    # Summary.
    elapsed = (datetime.now(timezone.utc) - start_time).total_seconds()
    LOGGER.info("=" * 60)
    LOGGER.info("SX03 Cross-Source Evaluation Summary (%.0fs)", elapsed)
    LOGGER.info("=" * 60)
    for arm_id, ci_dict in results.items():
        ndcg = ci_dict.get("holdout_ndcg")
        mAP = ci_dict.get("holdout_map")
        auc = ci_dict.get("holdout_roc_auc")
        if ndcg and ndcg.point_estimate is not None:
            LOGGER.info(
                "  %s: nDCG=%.4f [%.3f,%.3f], mAP=%.4f, AUC=%.4f",
                arm_id,
                ndcg.point_estimate,
                ndcg.ci_low or 0,
                ndcg.ci_high or 0,
                (mAP.point_estimate if mAP else 0) or 0,
                (auc.point_estimate if auc else 0) or 0,
            )


def parse_args(argv: list[str] | None = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--device-type", choices=("cpu", "gpu"), default="cpu")
    parser.add_argument("--cache-dir", type=Path, default=DEFAULT_CACHE_DIR)
    parser.add_argument("--candidate-dir", type=Path, default=DEFAULT_CANDIDATE_DIR)
    parser.add_argument("--output-dir", type=Path, default=DEFAULT_OUTPUT_DIR)
    return parser.parse_args(argv)


def main(argv: list[str] | None = None) -> None:
    setup_logging()
    args = parse_args(argv)
    LOGGER.info("SX03 eval starting at %s", datetime.now(timezone.utc).isoformat())

    candidate_module = load_module_from_path("sx03_candidate", args.candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=args.cache_dir, include_host_defense=True)

    run_sx03_eval(
        candidate_module=candidate_module,
        context=context,
        device_type=args.device_type,
        output_dir=args.output_dir,
    )


if __name__ == "__main__":
    main()
