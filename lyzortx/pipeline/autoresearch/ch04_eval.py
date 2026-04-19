#!/usr/bin/env python3
"""CH04: CHISEL baseline — per-row binary training, concentration as feature, AUC+Brier.

Establishes the first canonical CHISEL baseline. Training unit flips from pair-level
`any_lysis` (SPANDEX/CH02/CH03) to per-row binary `score`: every
(bacterium, phage, log_dilution, replicate, X, Y) raw observation with score ∈ {"0", "1"}
is its own training row. Rows with score == "n" are dropped as missing (not negative),
matching the ST01B rule.

Concentration enters the model as a numeric pair-level feature encoded as absolute
log₁₀ pfu/ml (`pair_concentration__log10_pfu_ml`). Guelin's {0, -1, -2, -4} relative
log_dilution steps map to {8.7, 7.7, 6.7, 4.7} absolute log₁₀ pfu/ml (Guelin neat at
5×10⁸ pfu/ml ≈ 10⁸·⁷). This encoding is chosen so CH05's unified Guelin+BASEL panel
can express BASEL's >10⁹ pfu/ml spot test as log₁₀ = 9.0 on the same feature axis
without the relative-log_dilution ambiguity (BASEL is not a dilution step). For
CH04 Guelin-only, the encoding is an affine shift of the prior `log_dilution`
feature — tree splits are threshold-equivalent and metrics are bit-identical to
pre-encoding runs under the same seeds. No binning, no conditioning, no
per-concentration sub-models — LightGBM learns the effect as part of training. All
other feature engineering is identical to the SX10 baseline (host_surface +
host_typing + host_stats + host_defense + phage_projection + phage_stats +
pair_depo_capsule + pair_receptor_omp, RFE-selected), so the CH04-vs-CH02 delta
isolates the frame change.

Per-phage blending (AX02) is dropped. SPANDEX's per-phage models added ~2 pp AUC to
bacteria-axis evaluation but are `per-phage-not-deployable` — they require training
data for each phage and cannot predict for unseen phages, which contradicts the
CHISEL `deployment-goal`. The SPANDEX bacteria-axis gain was also inflated by the
pre-CH02 fold-leakage bug (45 of 48 multi-bacterium cv_groups split across folds
let per-phage models memorize bacterium-level priors). CH05's phage-axis evaluation
measures the cross-phage generalization cost directly and is the honest successor to
"how much does per-phage help?" CHISEL's baseline is all-pairs only.

Evaluation uses AUC + Brier only. No nDCG, mAP, or top-k: ranking is a product-layer
concern and retired from the CHISEL scorecard (see ranking-metrics-retired). Each
held-out (bacterium, phage) pair is scored at its highest-observed concentration
(max log10_pfu_ml = Guelin neat 8.7 or BASEL 9.0) to match the deployment question
"would this phage lyse this bacterium at the maximum testable titer?" — that single
prediction per pair is aggregated with bacterium-level bootstrap CIs.

Usage:
    PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch04_eval --device-type cpu
"""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from types import ModuleType
from typing import Any, Optional, Sequence

import numpy as np
import pandas as pd
from sklearn.metrics import brier_score_loss, roc_auc_score

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import (
    BootstrapMetricCI,
    load_module_from_path,
    safe_round,
    temporary_module_attribute,
)
from lyzortx.pipeline.autoresearch.ch03_row_expansion import (
    RAW_SCORE_NEGATIVE,
    RAW_SCORE_POSITIVE,
    RAW_SCORE_UNINTERPRETABLE,
    load_row_expanded_frame,
)
from lyzortx.pipeline.autoresearch.gt03_eval import apply_rfe
from lyzortx.pipeline.autoresearch.sx01_eval import (
    FOLD_HASH_NAMESPACE,
    FOLD_SALT,
    N_FOLDS,
    SEEDS,
    assign_bacteria_folds,
    bacteria_to_cv_group_map,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_CACHE_DIR = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1")
DEFAULT_CANDIDATE_DIR = Path("lyzortx/autoresearch")
DEFAULT_OUTPUT_DIR = Path("lyzortx/generated_outputs/ch04_chisel_baseline")

CH02_REVALIDATED_METRICS_PATH = Path("lyzortx/generated_outputs/ch02_cv_group_fix/ch02_sx10_revalidated_metrics.json")
CONCENTRATION_FEATURE_COLUMN = "pair_concentration__log10_pfu_ml"
BOOTSTRAP_SAMPLES = 1000
BOOTSTRAP_RANDOM_STATE = 42


def build_clean_row_training_frame(row_frame: pd.DataFrame) -> pd.DataFrame:
    """Drop score=='n' rows, cast score to {0, 1}, attach CH04 label + concentration feature.

    `score == "n"` rows are dropped from training (they're missing observations, not
    negative). The remaining rows get two derived columns: `label_row_binary` (int 0/1)
    for the training target and `pair_concentration__log10_pfu_ml` (float) for the
    concentration feature LightGBM ingests. The row frame is expected to carry a
    pre-computed `log10_pfu_ml` column (added at the row-frame source: CH03 for Guelin,
    CH05 for BASEL) — CH04 does not derive it here, so each source owns its encoding.

    Also logs pair-level bookkeeping: pairs whose every observation is `n` (no
    interpretable lysis observations) lose all their rows in this filter, so the
    surviving pair count is below the input pair count. Log lets callers notice when
    that number shifts across dataset revisions.
    """
    clean = row_frame.loc[row_frame["score"].isin([RAW_SCORE_POSITIVE, RAW_SCORE_NEGATIVE])].copy()
    dropped_rows = len(row_frame) - len(clean)
    input_pairs = row_frame["pair_id"].nunique() if "pair_id" in row_frame.columns else None
    output_pairs = clean["pair_id"].nunique() if "pair_id" in clean.columns else None
    pair_drop_detail = ""
    if input_pairs is not None and output_pairs is not None:
        pairs_dropped = input_pairs - output_pairs
        pair_drop_detail = f"; {pairs_dropped} pairs lost all interpretable rows ({input_pairs} → {output_pairs} pairs)"
    LOGGER.info(
        "Dropped %d score='%s' rows (missing observations); %d interpretable rows remain%s",
        dropped_rows,
        RAW_SCORE_UNINTERPRETABLE,
        len(clean),
        pair_drop_detail,
    )
    if "log10_pfu_ml" not in clean.columns:
        raise ValueError(
            "Row frame missing required 'log10_pfu_ml' column. Guelin rows get it from "
            "`ch03_row_expansion.load_row_expanded_frame`; BASEL rows get it from "
            "`ch05_eval.load_basel_as_row_frame`."
        )
    clean["label_row_binary"] = clean["score"].astype(int)
    clean[CONCENTRATION_FEATURE_COLUMN] = clean["log10_pfu_ml"].astype(float)
    return clean


def train_and_predict_per_row_fold(
    *,
    candidate_module: ModuleType,
    context: Any,
    training_frame: pd.DataFrame,
    holdout_frame: pd.DataFrame,
    seed: int,
    device_type: str,
) -> tuple[list[dict[str, object]], pd.DataFrame]:
    """Row-level all-pairs trainer for CH04.

    Same host/phage slot bundle as SX10, same pairwise cross-terms, same RFE. Differs
    from SPANDEX/CH02 on three axes: (1) `label_row_binary` replaces `label_any_lysis`
    as the training target, (2) `pair_concentration__log10_pfu_ml` is added as a
    feature column before RFE, and (3) every row is a training example rather than
    one row per pair. Per-phage blending (AX02) is intentionally omitted — see module
    docstring for rationale.
    """
    from lyzortx.pipeline.autoresearch.derive_pairwise_depo_capsule_features import (
        compute_pairwise_depo_capsule_features,
    )
    from lyzortx.pipeline.autoresearch.derive_pairwise_receptor_omp_features import (
        compute_pairwise_receptor_omp_features,
    )

    host_slots = ["host_surface", "host_typing", "host_stats", "host_defense"]
    phage_slots = ["phage_projection", "phage_stats"]

    host_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=host_slots, entity_key="bacteria"
    )
    phage_table = candidate_module.build_entity_feature_table(
        context.slot_artifacts, slot_names=phage_slots, entity_key="phage"
    )
    host_typed, _, host_categorical = candidate_module.type_entity_features(host_table, "bacteria")
    phage_typed, _, phage_categorical = candidate_module.type_entity_features(phage_table, "phage")

    train_design = candidate_module.build_raw_pair_design_matrix(
        training_frame, host_features=host_typed, phage_features=phage_typed
    )
    holdout_design = candidate_module.build_raw_pair_design_matrix(
        holdout_frame, host_features=host_typed, phage_features=phage_typed
    )

    compute_pairwise_depo_capsule_features(train_design)
    compute_pairwise_depo_capsule_features(holdout_design)
    compute_pairwise_receptor_omp_features(train_design)
    compute_pairwise_receptor_omp_features(holdout_design)

    # Propagate row-level columns onto the merged design (many_to_one merges may drop
    # row-level columns in some candidate implementations; be explicit).
    train_design[CONCENTRATION_FEATURE_COLUMN] = training_frame[CONCENTRATION_FEATURE_COLUMN].to_numpy()
    holdout_design[CONCENTRATION_FEATURE_COLUMN] = holdout_frame[CONCENTRATION_FEATURE_COLUMN].to_numpy()
    train_design["label_row_binary"] = training_frame["label_row_binary"].to_numpy()
    holdout_design["label_row_binary"] = holdout_frame["label_row_binary"].to_numpy()
    train_design["log10_pfu_ml"] = training_frame["log10_pfu_ml"].to_numpy()
    holdout_design["log10_pfu_ml"] = holdout_frame["log10_pfu_ml"].to_numpy()
    train_design["log_dilution"] = training_frame["log_dilution"].to_numpy()
    holdout_design["log_dilution"] = holdout_frame["log_dilution"].to_numpy()
    train_design["replicate"] = training_frame["replicate"].to_numpy()
    holdout_design["replicate"] = holdout_frame["replicate"].to_numpy()

    prefixes = tuple(f"{s}__" for s in host_slots + phage_slots) + (
        "pair_depo_capsule__",
        "pair_receptor_omp__",
        "pair_concentration__",
    )
    feature_columns = [col for col in train_design.columns if col.startswith(prefixes)]
    categorical_columns = [col for col in (host_categorical + phage_categorical) if col in feature_columns]

    y_train = train_design["label_row_binary"].astype(int).to_numpy(dtype=int)
    rfe_features = apply_rfe(train_design, feature_columns, categorical_columns, y_train, seed=42)
    rfe_categorical = [c for c in categorical_columns if c in rfe_features]

    LOGGER.info(
        "Fold training (per-row): %d features after RFE (from %d), %d train rows, %d holdout rows, "
        "concentration_in_rfe=%s",
        len(rfe_features),
        len(feature_columns),
        len(train_design),
        len(holdout_design),
        CONCENTRATION_FEATURE_COLUMN in rfe_features,
    )

    with temporary_module_attribute(candidate_module, "PAIR_SCORER_RANDOM_STATE", seed):
        estimator = candidate_module.build_pair_scorer(device_type=device_type)
    estimator.fit(
        train_design[rfe_features],
        y_train,
        categorical_feature=rfe_categorical,
    )
    predictions = estimator.predict_proba(holdout_design[rfe_features])[:, 1]

    rows = []
    for (_, row), probability in zip(
        holdout_design[
            ["pair_id", "bacteria", "phage", "label_row_binary", "log_dilution", "log10_pfu_ml", "replicate"]
        ].iterrows(),
        predictions,
    ):
        rows.append(
            {
                "arm_id": "chisel_baseline",
                "seed": seed,
                "pair_id": str(row["pair_id"]),
                "bacteria": str(row["bacteria"]),
                "phage": str(row["phage"]),
                "log_dilution": int(row["log_dilution"]),
                "log10_pfu_ml": float(row["log10_pfu_ml"]),
                "replicate": int(row["replicate"]),
                "label_row_binary": int(row["label_row_binary"]),
                "predicted_probability": safe_round(float(probability)),
            }
        )

    feature_importance = pd.DataFrame(
        {
            "feature": rfe_features,
            "importance": estimator.feature_importances_,
        }
    )
    feature_importance["is_concentration_feature"] = feature_importance["feature"] == CONCENTRATION_FEATURE_COLUMN
    return rows, feature_importance


def select_pair_max_concentration_rows(per_row_predictions: pd.DataFrame) -> pd.DataFrame:
    """Keep one prediction row per held-out (bacterium, phage) pair, at max log10_pfu_ml.

    "Highest observed concentration" per pair is picked by max `log10_pfu_ml` — the
    absolute-titer feature the model actually sees. For Guelin-only inputs this is
    equivalent to max log_dilution (values {8.7, 7.7, 6.7, 4.7} order-preserve {0, -1,
    -2, -4}). For the unified Guelin+BASEL frame, BASEL's single row per pair always
    wins at log10_pfu_ml = 9.0, so the aggregation trivially surfaces it. Replicate
    rows at that top concentration are collapsed by averaging the predicted
    probability (replicates have identical feature vectors, so their predictions are
    identical anyway — the mean is a no-op in practice but it documents intent) and
    by max-aggregating the binary label (any positive replicate makes the pair
    positive). This "any replicate positive" rule matches the deployment semantics
    of the test strip read-out: a pair is declared positive if lysis is observed in
    at least one replicate at the top titer. CH04's core training still operates
    per-row.
    """
    ranked = per_row_predictions.sort_values(["pair_id", "log10_pfu_ml", "replicate"], ascending=[True, False, True])
    max_conc = ranked.groupby("pair_id", as_index=False)["log10_pfu_ml"].max()
    top_rows = ranked.merge(max_conc, on=["pair_id", "log10_pfu_ml"], how="inner")
    aggregated = top_rows.groupby(["pair_id", "bacteria", "phage", "log10_pfu_ml"], as_index=False).agg(
        predicted_probability=("predicted_probability", "mean"),
        label_row_binary=("label_row_binary", "max"),
        n_replicates_at_max=("replicate", "count"),
    )
    return aggregated


def compute_aggregate_auc_brier(pair_rows: Sequence[dict[str, object]]) -> dict[str, float]:
    labels = np.array([int(r["label_row_binary"]) for r in pair_rows])
    preds = np.array([float(r["predicted_probability"]) for r in pair_rows])
    return {
        "auc": float(roc_auc_score(labels, preds)),
        "brier": float(brier_score_loss(labels, preds)),
    }


def bootstrap_auc_brier_by_bacterium(
    pair_rows: Sequence[dict[str, object]],
    *,
    bootstrap_samples: int,
    bootstrap_random_state: int,
) -> dict[str, BootstrapMetricCI]:
    by_bacterium: dict[str, list[dict[str, object]]] = defaultdict(list)
    for row in pair_rows:
        by_bacterium[str(row["bacteria"])].append(row)
    bacteria_ids = tuple(sorted(by_bacterium.keys()))
    n_bacteria = len(bacteria_ids)
    rng = np.random.default_rng(bootstrap_random_state)

    aucs: list[float] = []
    briers: list[float] = []
    progress_interval = max(1, bootstrap_samples // 5)
    for i in range(bootstrap_samples):
        if i == 0 or (i + 1) % progress_interval == 0 or i + 1 == bootstrap_samples:
            LOGGER.info("Bootstrap progress: %d/%d", i + 1, bootstrap_samples)
        idx = rng.integers(0, n_bacteria, size=n_bacteria)
        sampled: list[dict[str, object]] = []
        for j in idx.tolist():
            sampled.extend(by_bacterium[bacteria_ids[j]])
        labels = np.array([int(r["label_row_binary"]) for r in sampled])
        preds = np.array([float(r["predicted_probability"]) for r in sampled])
        if len(np.unique(labels)) < 2:
            continue  # degenerate resample — skip AUC
        aucs.append(float(roc_auc_score(labels, preds)))
        briers.append(float(brier_score_loss(labels, preds)))

    point = compute_aggregate_auc_brier(pair_rows)

    def _ci(values: Sequence[float]) -> tuple[Optional[float], Optional[float], int]:
        if not values:
            return None, None, 0
        arr = np.asarray(values, dtype=float)
        low, high = np.quantile(arr, [0.025, 0.975])
        return safe_round(float(low)), safe_round(float(high)), len(values)

    auc_low, auc_high, auc_used = _ci(aucs)
    brier_low, brier_high, brier_used = _ci(briers)
    return {
        "holdout_roc_auc": BootstrapMetricCI(
            point_estimate=point["auc"],
            ci_low=auc_low,
            ci_high=auc_high,
            bootstrap_samples_requested=bootstrap_samples,
            bootstrap_samples_used=auc_used,
        ),
        "holdout_brier_score": BootstrapMetricCI(
            point_estimate=point["brier"],
            ci_low=brier_low,
            ci_high=brier_high,
            bootstrap_samples_requested=bootstrap_samples,
            bootstrap_samples_used=brier_used,
        ),
    }


def run_ch04_eval(
    *,
    device_type: str,
    output_dir: Path,
    cache_dir: Path,
    candidate_dir: Path,
) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    start_time = datetime.now(timezone.utc)
    LOGGER.info("CH04 evaluation starting at %s", start_time.isoformat())

    row_frame = load_row_expanded_frame()
    clean_rows = build_clean_row_training_frame(row_frame)

    training = clean_rows[
        (clean_rows["split_holdout"] == "train_non_holdout") & (clean_rows["is_hard_trainable"] == "1")
    ].copy()
    holdout = clean_rows[
        (clean_rows["split_holdout"] == "holdout_test") & (clean_rows["is_hard_trainable"] == "1")
    ].copy()
    full_frame = pd.concat([training, holdout], ignore_index=True)
    LOGGER.info(
        "CH04 clean row frame: %d total rows, %d bacteria, %d pairs",
        len(full_frame),
        full_frame["bacteria"].nunique(),
        full_frame["pair_id"].nunique(),
    )

    # Reuse CH02 fold assignment on the clean row frame.
    mapping = bacteria_to_cv_group_map(full_frame)
    fold_assignments = assign_bacteria_folds(mapping)
    LOGGER.info(
        "k-fold CV: %d bacteria / %d folds, fold hash=%s",
        len(mapping),
        N_FOLDS,
        hashlib.sha256(f"{FOLD_SALT}:{FOLD_HASH_NAMESPACE}".encode()).hexdigest()[:8],
    )

    candidate_module = load_module_from_path("ch04_candidate", candidate_dir / "train.py")
    context = candidate_module.load_and_validate_cache(cache_dir=cache_dir, include_host_defense=True)

    all_per_row_predictions: list[dict[str, object]] = []
    all_feature_importance: list[pd.DataFrame] = []

    for fold_id in range(N_FOLDS):
        holdout_bacteria = {b for b, f in fold_assignments.items() if f == fold_id}
        train_bacteria = {b for b, f in fold_assignments.items() if f != fold_id}
        fold_train = full_frame[full_frame["bacteria"].isin(train_bacteria)].copy()
        fold_holdout = full_frame[full_frame["bacteria"].isin(holdout_bacteria)].copy()
        LOGGER.info(
            "=== Fold %d: %d train bacteria (%d rows), %d holdout bacteria (%d rows) ===",
            fold_id,
            len(train_bacteria),
            len(fold_train),
            len(holdout_bacteria),
            len(fold_holdout),
        )

        fold_seed_rows: list[dict[str, object]] = []
        for seed in SEEDS:
            rows, fi = train_and_predict_per_row_fold(
                candidate_module=candidate_module,
                context=context,
                training_frame=fold_train,
                holdout_frame=fold_holdout,
                seed=seed,
                device_type=device_type,
            )
            for r in rows:
                r["fold_id"] = fold_id
            fold_seed_rows.extend(rows)
            fi["fold_id"] = fold_id
            fi["seed"] = seed
            all_feature_importance.append(fi)

        per_row_df = pd.DataFrame(fold_seed_rows)
        aggregated = (
            per_row_df.groupby(
                [
                    "fold_id",
                    "pair_id",
                    "bacteria",
                    "phage",
                    "log_dilution",
                    "log10_pfu_ml",
                    "replicate",
                    "label_row_binary",
                ],
                as_index=False,
            )["predicted_probability"]
            .mean()
            .sort_values(["bacteria", "phage", "log10_pfu_ml", "replicate"])
        )
        all_per_row_predictions.extend(aggregated.to_dict(orient="records"))

        fold_pair_predictions = select_pair_max_concentration_rows(aggregated)
        fold_pair_rows = fold_pair_predictions.to_dict(orient="records")
        fold_point = compute_aggregate_auc_brier(fold_pair_rows)
        LOGGER.info(
            "Fold %d metrics: AUC=%.4f, Brier=%.4f, n_pairs=%d, pos_rate=%.3f",
            fold_id,
            fold_point["auc"],
            fold_point["brier"],
            len(fold_pair_rows),
            fold_pair_predictions["label_row_binary"].mean(),
        )

    per_row_df = pd.DataFrame(all_per_row_predictions)
    per_row_df.to_csv(output_dir / "ch04_per_row_predictions.csv", index=False)

    pair_predictions = select_pair_max_concentration_rows(per_row_df)
    pair_predictions.to_csv(output_dir / "ch04_predictions.csv", index=False)
    LOGGER.info(
        "Pair-level (max-conc) predictions: %d pairs, %d bacteria",
        len(pair_predictions),
        pair_predictions["bacteria"].nunique(),
    )

    pair_rows = pair_predictions.to_dict(orient="records")
    point = compute_aggregate_auc_brier(pair_rows)
    LOGGER.info("Aggregate AUC=%.4f, Brier=%.4f (n=%d pairs)", point["auc"], point["brier"], len(pair_rows))
    ci_results = bootstrap_auc_brier_by_bacterium(
        pair_rows,
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )
    for name, ci in ci_results.items():
        LOGGER.info("  %s: %.4f [%.4f, %.4f]", name, ci.point_estimate, ci.ci_low or 0, ci.ci_high or 0)

    feature_importance_df = pd.concat(all_feature_importance, ignore_index=True)
    fi_agg = (
        feature_importance_df.groupby(["feature", "is_concentration_feature"], as_index=False)
        .agg(mean_importance=("importance", "mean"), n_folds_selected=("importance", "count"))
        .sort_values("mean_importance", ascending=False)
    )
    fi_agg.to_csv(output_dir / "ch04_feature_importance.csv", index=False)
    concentration_rows = fi_agg[fi_agg["is_concentration_feature"]]
    concentration_importance = (
        float(concentration_rows["mean_importance"].iloc[0]) if not concentration_rows.empty else 0.0
    )
    LOGGER.info(
        "Concentration feature importance (mean across %d folds × seeds): %.2f",
        int(concentration_rows["n_folds_selected"].iloc[0]) if not concentration_rows.empty else 0,
        concentration_importance,
    )

    aggregate = {
        "task_id": "CH04",
        "scorecard": "AUC + Brier (nDCG/mAP/top-k retired; see ranking-metrics-retired)",
        "training_unit": "per-row (bacterium, phage, log_dilution, replicate, X, Y) with score ∈ {0, 1}",
        "evaluation_unit": "per-pair at max observed log10_pfu_ml (one prediction per held-out pair)",
        "n_bacteria_total": int(pair_predictions["bacteria"].nunique()),
        "n_pairs_evaluated": int(len(pair_predictions)),
        "n_training_rows_dropped_as_n": int(len(row_frame) - len(clean_rows)),
        "concentration_feature_column": CONCENTRATION_FEATURE_COLUMN,
        "concentration_mean_importance": safe_round(concentration_importance),
        "chisel_baseline": {
            "holdout_roc_auc": {
                "point_estimate": safe_round(ci_results["holdout_roc_auc"].point_estimate),
                "ci_low": ci_results["holdout_roc_auc"].ci_low,
                "ci_high": ci_results["holdout_roc_auc"].ci_high,
            },
            "holdout_brier_score": {
                "point_estimate": safe_round(ci_results["holdout_brier_score"].point_estimate),
                "ci_low": ci_results["holdout_brier_score"].ci_low,
                "ci_high": ci_results["holdout_brier_score"].ci_high,
            },
        },
        "comparison_ch02_revalidated": _compare_against_ch02(ci_results),
        "elapsed_seconds": round((datetime.now(timezone.utc) - start_time).total_seconds(), 1),
    }
    aggregate_path = output_dir / "ch04_aggregate_metrics.json"
    with open(aggregate_path, "w", encoding="utf-8") as f:
        json.dump(aggregate, f, indent=2)
    LOGGER.info("Wrote aggregate metrics: %s", aggregate_path)
    return aggregate


def _compare_against_ch02(ci_results: dict[str, BootstrapMetricCI]) -> dict[str, object]:
    if not CH02_REVALIDATED_METRICS_PATH.exists():
        return {"note": f"CH02 baseline not found at {CH02_REVALIDATED_METRICS_PATH}"}
    with open(CH02_REVALIDATED_METRICS_PATH, encoding="utf-8") as f:
        ch02 = json.load(f)
    ch02_auc = ch02["chisel_fixed_folds"]["holdout_roc_auc"]["point_estimate"]
    ch02_brier = ch02["chisel_fixed_folds"]["holdout_brier_score"]["point_estimate"]
    return {
        "ch02_auc": ch02_auc,
        "ch02_brier": ch02_brier,
        "delta_auc": safe_round(ci_results["holdout_roc_auc"].point_estimate - ch02_auc),
        "delta_brier": safe_round(ci_results["holdout_brier_score"].point_estimate - ch02_brier),
    }


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
    run_ch04_eval(
        device_type=args.device_type,
        output_dir=args.output_dir,
        cache_dir=args.cache_dir,
        candidate_dir=args.candidate_dir,
    )


if __name__ == "__main__":
    main()
