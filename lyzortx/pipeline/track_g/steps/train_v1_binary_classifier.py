#!/usr/bin/env python3
"""TG01: Train tuned LightGBM and logistic-regression models on the v1 expanded feature set."""

from __future__ import annotations

import argparse
import hashlib
import json
import logging
import warnings
from dataclasses import dataclass
from pathlib import Path
from typing import Any, Callable, Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from sklearn.feature_extraction import DictVectorizer
from sklearn.linear_model import LogisticRegression
from sklearn.metrics import average_precision_score, brier_score_loss, log_loss, roc_auc_score
from sklearn.model_selection import ParameterGrid

from lyzortx.pipeline.steel_thread_v0.io.write_outputs import ensure_directory, write_csv, write_json
from lyzortx.pipeline.steel_thread_v0.steps._io_helpers import parse_float, read_csv_rows, safe_round
from lyzortx.pipeline.steel_thread_v0.steps.st04_train_baselines import (
    CATEGORICAL_FEATURE_COLUMNS as V0_CATEGORICAL_FEATURE_COLUMNS,
)
from lyzortx.pipeline.steel_thread_v0.steps.st04_train_baselines import (
    NUMERIC_FEATURE_COLUMNS as V0_NUMERIC_FEATURE_COLUMNS,
)
from lyzortx.pipeline.track_c.steps.build_v1_host_feature_pair_table import EXTENDED_CATEGORICAL_COLUMNS
from lyzortx.pipeline.track_c.steps import build_v1_host_feature_pair_table
from lyzortx.pipeline.track_d import run_track_d
from lyzortx.pipeline.steel_thread_v0.steps import (
    st01_label_policy,
    st01b_confidence_tiers,
    st02_build_pair_table,
    st03_build_splits,
)

logger = logging.getLogger(__name__)

IDENTIFIER_COLUMNS: Tuple[str, ...] = ("pair_id", "bacteria", "phage")
TRACK_E_REQUIRED_BLOCKS: Tuple[Tuple[str, Path], ...] = (
    (
        "rbp_receptor_compatibility",
        Path(
            "lyzortx/generated_outputs/track_e/rbp_receptor_compatibility_feature_block/"
            "rbp_receptor_compatibility_features_v1.csv"
        ),
    ),
    (
        "isolation_host_distance",
        Path(
            "lyzortx/generated_outputs/track_e/isolation_host_distance_feature_block/"
            "isolation_host_distance_features_v1.csv"
        ),
    ),
)
LIGHTGBM_PARAMETER_GRID: Tuple[Dict[str, object], ...] = tuple(
    ParameterGrid(
        {
            "n_estimators": [150, 300],
            "learning_rate": [0.03, 0.05],
            "num_leaves": [15, 31],
            "min_child_samples": [10, 25],
        }
    )
)
LOGISTIC_PARAMETER_GRID: Tuple[Dict[str, object], ...] = tuple(ParameterGrid({"C": [0.1, 0.3, 1.0, 3.0, 10.0]}))


@dataclass(frozen=True)
class FeatureSpace:
    categorical_columns: Tuple[str, ...]
    numeric_columns: Tuple[str, ...]
    track_c_additional_columns: Tuple[str, ...]
    track_d_columns: Tuple[str, ...]
    track_e_columns: Tuple[str, ...]


@dataclass(frozen=True)
class FoldDataset:
    fold_id: int
    train_rows: Tuple[Dict[str, object], ...]
    valid_rows: Tuple[Dict[str, object], ...]
    y_train: Tuple[int, ...]
    y_valid: Tuple[int, ...]
    X_train: Any
    X_valid: Any
    sample_weights: Tuple[float, ...]


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--st02-pair-table-path",
        type=Path,
        default=Path("lyzortx/generated_outputs/steel_thread_v0/intermediate/st02_pair_table.csv"),
        help="Input ST0.2 pair table path.",
    )
    parser.add_argument(
        "--st03-split-assignments-path",
        type=Path,
        default=Path("lyzortx/generated_outputs/steel_thread_v0/intermediate/st03_split_assignments.csv"),
        help="Input ST0.3 split assignments path.",
    )
    parser.add_argument(
        "--track-c-pair-table-path",
        type=Path,
        default=Path("lyzortx/generated_outputs/track_c/v1_host_feature_pair_table/pair_table_v1.csv"),
        help="Input Track C v1 pair table path.",
    )
    parser.add_argument(
        "--track-d-genome-kmer-path",
        type=Path,
        default=Path("lyzortx/generated_outputs/track_d/phage_genome_kmer_features/phage_genome_kmer_features.csv"),
        help="Input Track D genome k-mer feature CSV.",
    )
    parser.add_argument(
        "--track-d-distance-path",
        type=Path,
        default=Path(
            "lyzortx/generated_outputs/track_d/phage_distance_embedding/phage_distance_embedding_features.csv"
        ),
        help="Input Track D phage-distance feature CSV.",
    )
    parser.add_argument(
        "--track-e-rbp-compatibility-path",
        type=Path,
        default=TRACK_E_REQUIRED_BLOCKS[0][1],
        help="Input Track E RBP-receptor compatibility feature CSV.",
    )
    parser.add_argument(
        "--track-e-isolation-distance-path",
        type=Path,
        default=TRACK_E_REQUIRED_BLOCKS[1][1],
        help="Input Track E isolation-host distance feature CSV.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lyzortx/generated_outputs/track_g/tg01_v1_binary_classifier"),
        help="Directory for generated TG01 artifacts.",
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=42,
        help="Base random seed for tuned models.",
    )
    parser.add_argument(
        "--skip-prerequisites",
        action="store_true",
        help="Assume prerequisite Track C/D/E outputs already exist instead of generating missing artifacts.",
    )
    return parser.parse_args(argv)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(65536), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _deduplicate_preserving_order(values: Iterable[str]) -> Tuple[str, ...]:
    out: List[str] = []
    seen: set[str] = set()
    for value in values:
        if value in seen:
            continue
        seen.add(value)
        out.append(value)
    return tuple(out)


def _build_feature_dict(
    row: Mapping[str, object],
    *,
    categorical_columns: Sequence[str],
    numeric_columns: Sequence[str],
) -> Dict[str, object]:
    features: Dict[str, object] = {}
    for column in categorical_columns:
        value = row.get(column, "")
        if value not in {"", None}:
            features[column] = str(value)
    for column in numeric_columns:
        raw = row.get(column, "")
        if raw in {"", None}:
            continue
        if isinstance(raw, (int, float)):
            features[column] = float(raw)
            continue
        parsed = parse_float(str(raw))
        if parsed is not None:
            features[column] = parsed
    return features


def deduplicate_preserving_order(values: Iterable[str]) -> Tuple[str, ...]:
    """Return `values` with duplicates removed while preserving the first occurrence."""
    return _deduplicate_preserving_order(values)


def build_feature_dict(
    row: Mapping[str, object],
    *,
    categorical_columns: Sequence[str],
    numeric_columns: Sequence[str],
) -> Dict[str, object]:
    """Public wrapper for the feature-vector dictionary builder used by downstream steps."""
    return _build_feature_dict(
        row,
        categorical_columns=categorical_columns,
        numeric_columns=numeric_columns,
    )


def compute_binary_metrics(y_true: Sequence[int], y_prob: Sequence[float]) -> Dict[str, Optional[float]]:
    if not y_true:
        raise ValueError("No labels available for metric computation.")

    metrics: Dict[str, Optional[float]] = {
        "n": float(len(y_true)),
        "positive_rate": safe_round(sum(y_true) / len(y_true)),
        "average_precision": None,
        "roc_auc": None,
        "brier_score": safe_round(brier_score_loss(y_true, y_prob)),
        "log_loss": safe_round(log_loss(y_true, y_prob, labels=[0, 1])),
    }
    if len(set(y_true)) >= 2:
        metrics["average_precision"] = safe_round(average_precision_score(y_true, y_prob))
        metrics["roc_auc"] = safe_round(roc_auc_score(y_true, y_prob))
    return metrics


def _predict_probabilities(estimator: Any, X: Any) -> List[float]:
    with warnings.catch_warnings():
        warnings.filterwarnings(
            "ignore",
            message="X does not have valid feature names, but LGBMClassifier was fitted with feature names",
            category=UserWarning,
        )
        return [float(value) for value in estimator.predict_proba(X)[:, 1]]


def predict_probabilities(estimator: Any, X: Any) -> List[float]:
    """Public wrapper for probability scoring used by TL05 and SHAP explanation steps."""
    return _predict_probabilities(estimator, X)


def compute_top3_hit_rate(rows: Sequence[Mapping[str, object]], *, probability_key: str) -> Dict[str, object]:
    rows_by_bacteria: Dict[str, List[Mapping[str, object]]] = {}
    for row in rows:
        rows_by_bacteria.setdefault(str(row["bacteria"]), []).append(row)

    total = 0
    hits = 0
    susceptible_total = 0
    susceptible_hits = 0
    for bacteria_rows in rows_by_bacteria.values():
        ranked = sorted(bacteria_rows, key=lambda row: (-float(row[probability_key]), str(row["phage"])))
        top3 = ranked[:3]
        is_hit = any(int(row["label_hard_any_lysis"]) == 1 for row in top3)
        total += 1
        hits += 1 if is_hit else 0
        is_susceptible = any(int(row["label_hard_any_lysis"]) == 1 for row in ranked)
        if is_susceptible:
            susceptible_total += 1
            susceptible_hits += 1 if is_hit else 0

    return {
        "strain_count": total,
        "hit_count": hits,
        "top3_hit_rate_all_strains": safe_round(hits / total if total else 0.0),
        "susceptible_strain_count": susceptible_total,
        "susceptible_hit_count": susceptible_hits,
        "top3_hit_rate_susceptible_only": safe_round(
            susceptible_hits / susceptible_total if susceptible_total else 0.0
        ),
    }


def select_best_candidate(results: Sequence[Mapping[str, object]]) -> Dict[str, object]:
    if not results:
        raise ValueError("No CV candidate results available.")

    def _sort_key(result: Mapping[str, object]) -> Tuple[float, float, float, str]:
        summary = result["summary"]
        roc_auc = float(summary["mean_roc_auc"] or float("-inf"))
        top3 = float(summary["mean_top3_hit_rate_all_strains"] or float("-inf"))
        brier = float(summary["mean_brier_score"] or float("inf"))
        params = json.dumps(result["params"], sort_keys=True)
        return (roc_auc, top3, -brier, params)

    return dict(max(results, key=_sort_key))


def prepare_fold_datasets(
    rows: Sequence[Mapping[str, object]],
    feature_space: FeatureSpace,
) -> List[FoldDataset]:
    training_rows = [
        dict(row)
        for row in rows
        if row["split_holdout"] == "train_non_holdout" and str(row["is_hard_trainable"]) == "1"
    ]
    if not training_rows:
        raise ValueError("No non-holdout hard-trainable rows available for TG01.")

    fold_ids = sorted(
        {int(str(row["split_cv5_fold"])) for row in training_rows if int(str(row["split_cv5_fold"])) >= 0}
    )
    if len(fold_ids) < 2:
        raise ValueError("TG01 requires at least two non-holdout folds.")

    datasets: List[FoldDataset] = []
    for fold_id in fold_ids:
        train_rows = tuple(row for row in training_rows if int(str(row["split_cv5_fold"])) != fold_id)
        valid_rows = tuple(row for row in training_rows if int(str(row["split_cv5_fold"])) == fold_id)
        if not train_rows or not valid_rows:
            raise ValueError(f"Fold {fold_id} is missing train or validation rows.")
        y_train = tuple(int(str(row["label_hard_any_lysis"])) for row in train_rows)
        y_valid = tuple(int(str(row["label_hard_any_lysis"])) for row in valid_rows)
        w_train = tuple(float(row.get("training_weight_v3", 1.0) or 1.0) for row in train_rows)
        vectorizer = DictVectorizer(sparse=True, sort=True)
        X_train = vectorizer.fit_transform(
            [
                _build_feature_dict(
                    row,
                    categorical_columns=feature_space.categorical_columns,
                    numeric_columns=feature_space.numeric_columns,
                )
                for row in train_rows
            ]
        )
        X_valid = vectorizer.transform(
            [
                _build_feature_dict(
                    row,
                    categorical_columns=feature_space.categorical_columns,
                    numeric_columns=feature_space.numeric_columns,
                )
                for row in valid_rows
            ]
        )
        datasets.append(
            FoldDataset(
                fold_id=fold_id,
                train_rows=train_rows,
                valid_rows=valid_rows,
                y_train=y_train,
                y_valid=y_valid,
                X_train=X_train,
                X_valid=X_valid,
                sample_weights=w_train,
            )
        )
    return datasets


def evaluate_candidate_grid(
    fold_datasets: Sequence[FoldDataset],
    *,
    candidate_params: Sequence[Mapping[str, object]],
    estimator_factory: Callable[[Mapping[str, object], int], Any],
    model_label: str,
) -> List[Dict[str, object]]:
    results: List[Dict[str, object]] = []
    for params in candidate_params:
        fold_metrics: List[Dict[str, object]] = []
        for dataset in fold_datasets:
            estimator = estimator_factory(params, dataset.fold_id)
            estimator.fit(dataset.X_train, dataset.y_train, sample_weight=dataset.sample_weights)
            probabilities = _predict_probabilities(estimator, dataset.X_valid)
            scored_rows = []
            for row, probability in zip(dataset.valid_rows, probabilities):
                scored = dict(row)
                scored["predicted_probability"] = probability
                scored_rows.append(scored)

            binary_metrics = compute_binary_metrics(dataset.y_valid, probabilities)
            top3_metrics = compute_top3_hit_rate(scored_rows, probability_key="predicted_probability")
            fold_metrics.append(
                {
                    "fold_id": dataset.fold_id,
                    "train_rows": len(dataset.train_rows),
                    "validation_rows": len(dataset.valid_rows),
                    "binary_metrics": binary_metrics,
                    "top3_metrics": top3_metrics,
                }
            )

        summary = summarize_fold_metrics(fold_metrics)
        results.append(
            {
                "model_label": model_label,
                "params": dict(params),
                "fold_metrics": fold_metrics,
                "summary": summary,
            }
        )
    return results


def summarize_fold_metrics(fold_metrics: Sequence[Mapping[str, object]]) -> Dict[str, Optional[float]]:
    if not fold_metrics:
        raise ValueError("No fold metrics to summarize.")

    def _mean(path: Tuple[str, ...]) -> Optional[float]:
        values: List[float] = []
        for metrics in fold_metrics:
            current: Any = metrics
            for key in path:
                current = current[key]
            if current is None:
                continue
            values.append(float(current))
        if not values:
            return None
        return safe_round(sum(values) / len(values))

    return {
        "mean_average_precision": _mean(("binary_metrics", "average_precision")),
        "mean_roc_auc": _mean(("binary_metrics", "roc_auc")),
        "mean_brier_score": _mean(("binary_metrics", "brier_score")),
        "mean_log_loss": _mean(("binary_metrics", "log_loss")),
        "mean_top3_hit_rate_all_strains": _mean(("top3_metrics", "top3_hit_rate_all_strains")),
        "mean_top3_hit_rate_susceptible_only": _mean(("top3_metrics", "top3_hit_rate_susceptible_only")),
    }


def score_rows_with_cv_predictions(
    fold_datasets: Sequence[FoldDataset],
    *,
    estimator_factory: Callable[[Mapping[str, object], int], Any],
    best_params: Mapping[str, object],
    probability_column: str,
) -> List[Dict[str, object]]:
    scored_rows: List[Dict[str, object]] = []
    for dataset in fold_datasets:
        estimator = estimator_factory(best_params, dataset.fold_id)
        estimator.fit(dataset.X_train, dataset.y_train, sample_weight=dataset.sample_weights)
        probabilities = _predict_probabilities(estimator, dataset.X_valid)
        for row, probability in zip(dataset.valid_rows, probabilities):
            scored = dict(row)
            scored[probability_column] = probability
            scored["prediction_context"] = "non_holdout_oof"
            scored_rows.append(scored)
    scored_rows.sort(key=lambda row: (str(row["bacteria"]), str(row["phage"])))
    return scored_rows


def build_feature_space(
    st02_rows: Sequence[Mapping[str, str]],
    track_c_pair_rows: Sequence[Mapping[str, str]],
    track_d_feature_columns: Sequence[str],
    track_e_feature_columns: Sequence[str],
) -> FeatureSpace:
    st02_columns = tuple(st02_rows[0].keys())
    track_c_new_columns = tuple(column for column in track_c_pair_rows[0].keys() if column not in st02_columns)
    track_c_categorical = tuple(column for column in track_c_new_columns if column in set(EXTENDED_CATEGORICAL_COLUMNS))
    track_c_numeric = tuple(column for column in track_c_new_columns if column not in set(track_c_categorical))
    categorical_columns = _deduplicate_preserving_order(
        list(V0_CATEGORICAL_FEATURE_COLUMNS) + list(track_c_categorical)
    )
    numeric_columns = _deduplicate_preserving_order(
        list(V0_NUMERIC_FEATURE_COLUMNS)
        + list(track_c_numeric)
        + list(track_d_feature_columns)
        + list(track_e_feature_columns)
    )
    return FeatureSpace(
        categorical_columns=categorical_columns,
        numeric_columns=numeric_columns,
        track_c_additional_columns=track_c_new_columns,
        track_d_columns=tuple(track_d_feature_columns),
        track_e_columns=tuple(track_e_feature_columns),
    )


def merge_expanded_feature_rows(
    track_c_pair_rows: Sequence[Mapping[str, str]],
    split_rows: Sequence[Mapping[str, str]],
    phage_feature_blocks: Sequence[Sequence[Mapping[str, str]]],
    pair_feature_blocks: Sequence[Sequence[Mapping[str, str]]],
    *,
    allow_missing_pair_features: bool = False,
) -> List[Dict[str, object]]:
    split_by_pair = {row["pair_id"]: row for row in split_rows}
    phage_indexes: List[Tuple[Dict[str, Dict[str, str]], Tuple[str, ...]]] = []
    for rows in phage_feature_blocks:
        if not rows:
            raise ValueError("Encountered empty phage feature block.")
        columns = tuple(column for column in rows[0].keys() if column != "phage")
        phage_indexes.append(({row["phage"]: dict(row) for row in rows}, columns))

    pair_indexes: List[Tuple[Dict[str, Dict[str, str]], Tuple[str, ...]]] = []
    for rows in pair_feature_blocks:
        if not rows:
            raise ValueError("Encountered empty pair feature block.")
        columns = tuple(column for column in rows[0].keys() if column not in IDENTIFIER_COLUMNS)
        pair_indexes.append(({row["pair_id"]: dict(row) for row in rows}, columns))

    output_rows: List[Dict[str, object]] = []
    for row in track_c_pair_rows:
        pair_id = row["pair_id"]
        merged = dict(row)
        split_row = split_by_pair.get(pair_id)
        if split_row is None:
            raise KeyError(f"Missing ST0.3 split assignment for pair_id {pair_id}")
        merged.update(split_row)

        for phage_index, columns in phage_indexes:
            phage_row = phage_index.get(row["phage"])
            if phage_row is None:
                raise KeyError(f"Missing phage-level feature row for phage {row['phage']}")
            for column in columns:
                merged[column] = phage_row[column]

        for pair_index, columns in pair_indexes:
            pair_row = pair_index.get(pair_id)
            if pair_row is None:
                if not allow_missing_pair_features or str(merged["split_holdout"]) != "holdout_test":
                    raise KeyError(f"Missing pair-level feature row for pair_id {pair_id}")
                for column in columns:
                    merged[column] = 0.0
                continue
            for column in columns:
                merged[column] = pair_row[column]

        output_rows.append(merged)
    output_rows.sort(key=lambda row: (str(row["bacteria"]), str(row["phage"])))
    return output_rows


def ensure_prerequisite_outputs(args: argparse.Namespace) -> None:
    if args.skip_prerequisites:
        return

    st01_output_path = Path("lyzortx/generated_outputs/steel_thread_v0/intermediate/st01_pair_label_audit.csv")
    st01b_output_path = Path("lyzortx/generated_outputs/steel_thread_v0/intermediate/st01b_pair_confidence_audit.csv")
    if not st01_output_path.exists():
        st01_label_policy.main([])
    if not st01b_output_path.exists():
        st01b_confidence_tiers.main([])
    if not args.st02_pair_table_path.exists():
        st02_build_pair_table.main([])
    if not args.st03_split_assignments_path.exists():
        st03_build_splits.main([])
    if not args.track_c_pair_table_path.exists():
        build_v1_host_feature_pair_table.main([])
    if not args.track_d_genome_kmer_path.exists() or not args.track_d_distance_path.exists():
        run_track_d.main(["--step", "all"])
    if not args.track_e_rbp_compatibility_path.exists() or not args.track_e_isolation_distance_path.exists():
        raise FileNotFoundError(
            "Track E feature artifacts missing. Track E is a dead-ended track and cannot be rebuilt automatically."
        )


def make_lightgbm_estimator(params: Mapping[str, object], seed_offset: int, *, base_random_state: int) -> Any:
    from lightgbm import LGBMClassifier

    return LGBMClassifier(
        objective="binary",
        class_weight="balanced",
        random_state=base_random_state + seed_offset,
        deterministic=True,
        verbosity=-1,
        force_col_wise=True,
        colsample_bytree=0.8,
        subsample=0.8,
        reg_lambda=0.0,
        **params,
    )


def make_logistic_estimator(
    params: Mapping[str, object],
    seed_offset: int,
    *,
    base_random_state: int,
) -> LogisticRegression:
    return LogisticRegression(
        solver="liblinear",
        class_weight="balanced",
        max_iter=3000,
        random_state=base_random_state + seed_offset,
        **params,
    )


def fit_final_estimator(
    rows: Sequence[Mapping[str, object]],
    feature_space: FeatureSpace,
    *,
    estimator_factory: Callable[[Mapping[str, object], int], Any],
    params: Mapping[str, object],
    sample_weight_key: Optional[str] = None,
) -> Tuple[Any, DictVectorizer, List[Dict[str, object]], List[Dict[str, object]], List[float]]:
    train_rows = [
        dict(row)
        for row in rows
        if row["split_holdout"] == "train_non_holdout" and str(row["is_hard_trainable"]) == "1"
    ]
    eval_rows = [
        dict(row) for row in rows if row["split_holdout"] == "holdout_test" and str(row["is_hard_trainable"]) == "1"
    ]
    if not train_rows:
        raise ValueError("No non-holdout hard-trainable rows available for final training.")
    if not eval_rows:
        raise ValueError("No holdout hard-trainable rows available for final evaluation.")

    vectorizer = DictVectorizer(sparse=True, sort=True)
    X_train = vectorizer.fit_transform(
        [
            _build_feature_dict(
                row,
                categorical_columns=feature_space.categorical_columns,
                numeric_columns=feature_space.numeric_columns,
            )
            for row in train_rows
        ]
    )
    y_train = [int(str(row["label_hard_any_lysis"])) for row in train_rows]
    sample_weights = None
    if sample_weight_key is not None:
        sample_weights = [float(row.get(sample_weight_key, 1.0) or 1.0) for row in train_rows]
    X_eval = vectorizer.transform(
        [
            _build_feature_dict(
                row,
                categorical_columns=feature_space.categorical_columns,
                numeric_columns=feature_space.numeric_columns,
            )
            for row in eval_rows
        ]
    )
    estimator = estimator_factory(params, 0)
    if sample_weights is None:
        estimator.fit(X_train, y_train)
    else:
        estimator.fit(X_train, y_train, sample_weight=sample_weights)
    probabilities = _predict_probabilities(estimator, X_eval)
    return estimator, vectorizer, train_rows, eval_rows, probabilities


def build_top3_ranking_rows(
    rows: Sequence[Mapping[str, object]],
    *,
    probability_key: str,
    model_label: str,
) -> List[Dict[str, object]]:
    rows_by_bacteria: Dict[str, List[Mapping[str, object]]] = {}
    for row in rows:
        rows_by_bacteria.setdefault(str(row["bacteria"]), []).append(row)

    ranking_rows: List[Dict[str, object]] = []
    for bacteria, bacteria_rows in sorted(rows_by_bacteria.items()):
        ranked = sorted(bacteria_rows, key=lambda row: (-float(row[probability_key]), str(row["phage"])))
        for rank, row in enumerate(ranked[:3], start=1):
            ranking_rows.append(
                {
                    "model_label": model_label,
                    "bacteria": bacteria,
                    "phage": row["phage"],
                    "pair_id": row["pair_id"],
                    "rank": rank,
                    "predicted_probability": safe_round(float(row[probability_key])),
                    "label_hard_any_lysis": row["label_hard_any_lysis"],
                }
            )
    return ranking_rows


def flatten_candidate_rows(model_label: str, results: Sequence[Mapping[str, object]]) -> List[Dict[str, object]]:
    rows: List[Dict[str, object]] = []
    for result in results:
        summary = result["summary"]
        rows.append(
            {
                "model_label": model_label,
                "params_json": json.dumps(result["params"], sort_keys=True),
                "mean_average_precision": summary["mean_average_precision"],
                "mean_roc_auc": summary["mean_roc_auc"],
                "mean_brier_score": summary["mean_brier_score"],
                "mean_log_loss": summary["mean_log_loss"],
                "mean_top3_hit_rate_all_strains": summary["mean_top3_hit_rate_all_strains"],
                "mean_top3_hit_rate_susceptible_only": summary["mean_top3_hit_rate_susceptible_only"],
            }
        )
    return rows


def project_rows_to_fields(rows: Sequence[Mapping[str, object]], fieldnames: Sequence[str]) -> List[Dict[str, object]]:
    return [{fieldname: row.get(fieldname, "") for fieldname in fieldnames} for row in rows]


def main(argv: Optional[Sequence[str]] = None) -> int:
    args = parse_args(argv)
    logger.info("TG01 starting: train v1 binary classifier")
    ensure_directory(args.output_dir)
    ensure_prerequisite_outputs(args)

    st02_rows = read_csv_rows(args.st02_pair_table_path)
    split_rows = read_csv_rows(args.st03_split_assignments_path)
    track_c_pair_rows = read_csv_rows(args.track_c_pair_table_path)
    track_d_genome_rows = read_csv_rows(args.track_d_genome_kmer_path)
    track_d_distance_rows = read_csv_rows(args.track_d_distance_path)
    track_e_rbp_rows = read_csv_rows(args.track_e_rbp_compatibility_path)
    track_e_isolation_rows = read_csv_rows(args.track_e_isolation_distance_path)

    track_d_feature_columns = _deduplicate_preserving_order(
        [column for column in track_d_genome_rows[0].keys() if column != "phage"]
        + [column for column in track_d_distance_rows[0].keys() if column != "phage"]
    )
    track_e_feature_columns = _deduplicate_preserving_order(
        [column for column in track_e_rbp_rows[0].keys() if column not in IDENTIFIER_COLUMNS]
        + [column for column in track_e_isolation_rows[0].keys() if column not in IDENTIFIER_COLUMNS]
    )
    feature_space = build_feature_space(
        st02_rows,
        track_c_pair_rows,
        track_d_feature_columns,
        track_e_feature_columns,
    )
    merged_rows = merge_expanded_feature_rows(
        track_c_pair_rows,
        split_rows,
        phage_feature_blocks=(track_d_genome_rows, track_d_distance_rows),
        pair_feature_blocks=(track_e_rbp_rows, track_e_isolation_rows),
    )
    fold_datasets = prepare_fold_datasets(merged_rows, feature_space)
    lightgbm_factory = lambda params, seed_offset: make_lightgbm_estimator(  # noqa: E731
        params,
        seed_offset,
        base_random_state=args.random_state,
    )
    logreg_factory = lambda params, seed_offset: make_logistic_estimator(  # noqa: E731
        params,
        seed_offset,
        base_random_state=args.random_state,
    )

    lightgbm_results = evaluate_candidate_grid(
        fold_datasets,
        candidate_params=LIGHTGBM_PARAMETER_GRID,
        estimator_factory=lightgbm_factory,
        model_label="lightgbm",
    )
    logreg_results = evaluate_candidate_grid(
        fold_datasets,
        candidate_params=LOGISTIC_PARAMETER_GRID,
        estimator_factory=logreg_factory,
        model_label="logistic_regression",
    )
    best_lightgbm = select_best_candidate(lightgbm_results)
    best_logreg = select_best_candidate(logreg_results)

    cv_scored_lightgbm = score_rows_with_cv_predictions(
        fold_datasets,
        estimator_factory=lightgbm_factory,
        best_params=best_lightgbm["params"],
        probability_column="lightgbm_probability",
    )
    cv_scored_logreg = score_rows_with_cv_predictions(
        fold_datasets,
        estimator_factory=logreg_factory,
        best_params=best_logreg["params"],
        probability_column="logreg_probability",
    )
    cv_scored_by_pair = {row["pair_id"]: dict(row) for row in cv_scored_lightgbm}
    for row in cv_scored_logreg:
        cv_scored_by_pair[row["pair_id"]]["logreg_probability"] = row["logreg_probability"]
    cv_prediction_rows = [cv_scored_by_pair[pair_id] for pair_id in sorted(cv_scored_by_pair)]

    _, _, _, holdout_rows_lightgbm, holdout_lightgbm_prob = fit_final_estimator(
        merged_rows,
        feature_space,
        estimator_factory=lightgbm_factory,
        params=best_lightgbm["params"],
        sample_weight_key="training_weight_v3",
    )
    _, _, _, holdout_rows_logreg, holdout_logreg_prob = fit_final_estimator(
        merged_rows,
        feature_space,
        estimator_factory=logreg_factory,
        params=best_logreg["params"],
        sample_weight_key="training_weight_v3",
    )
    holdout_prediction_rows: List[Dict[str, object]] = []
    for lightgbm_row, logreg_row, lightgbm_prob, logreg_prob in zip(
        holdout_rows_lightgbm,
        holdout_rows_logreg,
        holdout_lightgbm_prob,
        holdout_logreg_prob,
    ):
        if lightgbm_row["pair_id"] != logreg_row["pair_id"]:
            raise ValueError("Holdout row alignment mismatch between LightGBM and logistic regression.")
        scored = dict(lightgbm_row)
        scored["lightgbm_probability"] = lightgbm_prob
        scored["logreg_probability"] = logreg_prob
        scored["prediction_context"] = "holdout_final"
        holdout_prediction_rows.append(scored)

    holdout_y = [int(str(row["label_hard_any_lysis"])) for row in holdout_prediction_rows]
    lightgbm_holdout_metrics = compute_binary_metrics(
        holdout_y,
        [float(row["lightgbm_probability"]) for row in holdout_prediction_rows],
    )
    logreg_holdout_metrics = compute_binary_metrics(
        holdout_y,
        [float(row["logreg_probability"]) for row in holdout_prediction_rows],
    )
    lightgbm_holdout_top3 = compute_top3_hit_rate(holdout_prediction_rows, probability_key="lightgbm_probability")
    logreg_holdout_top3 = compute_top3_hit_rate(holdout_prediction_rows, probability_key="logreg_probability")

    pair_prediction_rows = cv_prediction_rows + holdout_prediction_rows
    pair_prediction_rows.sort(key=lambda row: (str(row["prediction_context"]), str(row["bacteria"]), str(row["phage"])))

    top3_ranking_rows = build_top3_ranking_rows(
        holdout_prediction_rows,
        probability_key="lightgbm_probability",
        model_label="lightgbm",
    ) + build_top3_ranking_rows(
        holdout_prediction_rows,
        probability_key="logreg_probability",
        model_label="logistic_regression",
    )
    top3_ranking_rows.sort(key=lambda row: (str(row["model_label"]), str(row["bacteria"]), int(row["rank"])))

    summary = {
        "task_id": "TG01",
        "feature_space": {
            "categorical_columns": list(feature_space.categorical_columns),
            "numeric_columns": list(feature_space.numeric_columns),
            "track_c_additional_columns": list(feature_space.track_c_additional_columns),
            "track_d_columns": list(feature_space.track_d_columns),
            "track_e_columns": list(feature_space.track_e_columns),
            "categorical_feature_count": len(feature_space.categorical_columns),
            "numeric_feature_count": len(feature_space.numeric_columns),
        },
        "cv_protocol": {
            "group_source": "ST0.2 cv_group via canonical ST0.3 split assignments",
            "fold_ids": sorted(
                {int(str(row["split_cv5_fold"])) for row in merged_rows if int(str(row["split_cv5_fold"])) >= 0}
            ),
            "holdout_split": "ST0.3 holdout_test",
            "training_subset": "split_holdout=train_non_holdout and is_hard_trainable=1",
            "candidate_counts": {
                "lightgbm": len(LIGHTGBM_PARAMETER_GRID),
                "logistic_regression": len(LOGISTIC_PARAMETER_GRID),
            },
        },
        "lightgbm": {
            "best_params": best_lightgbm["params"],
            "cv_summary": best_lightgbm["summary"],
            "holdout_binary_metrics": lightgbm_holdout_metrics,
            "holdout_top3_metrics": lightgbm_holdout_top3,
        },
        "logistic_regression": {
            "best_params": best_logreg["params"],
            "cv_summary": best_logreg["summary"],
            "holdout_binary_metrics": logreg_holdout_metrics,
            "holdout_top3_metrics": logreg_holdout_top3,
        },
        "inputs": {
            "st02_pair_table": {"path": str(args.st02_pair_table_path), "sha256": _sha256(args.st02_pair_table_path)},
            "st03_split_assignments": {
                "path": str(args.st03_split_assignments_path),
                "sha256": _sha256(args.st03_split_assignments_path),
            },
            "track_c_pair_table": {
                "path": str(args.track_c_pair_table_path),
                "sha256": _sha256(args.track_c_pair_table_path),
            },
            "track_d_genome_kmers": {
                "path": str(args.track_d_genome_kmer_path),
                "sha256": _sha256(args.track_d_genome_kmer_path),
            },
            "track_d_distance": {
                "path": str(args.track_d_distance_path),
                "sha256": _sha256(args.track_d_distance_path),
            },
            "track_e_rbp_receptor_compatibility": {
                "path": str(args.track_e_rbp_compatibility_path),
                "sha256": _sha256(args.track_e_rbp_compatibility_path),
            },
            "track_e_isolation_host_distance": {
                "path": str(args.track_e_isolation_distance_path),
                "sha256": _sha256(args.track_e_isolation_distance_path),
            },
        },
    }

    write_json(args.output_dir / "tg01_model_summary.json", summary)
    write_csv(
        args.output_dir / "tg01_cv_candidate_results.csv",
        [
            "model_label",
            "params_json",
            "mean_average_precision",
            "mean_roc_auc",
            "mean_brier_score",
            "mean_log_loss",
            "mean_top3_hit_rate_all_strains",
            "mean_top3_hit_rate_susceptible_only",
        ],
        flatten_candidate_rows("lightgbm", lightgbm_results)
        + flatten_candidate_rows("logistic_regression", logreg_results),
    )
    pair_prediction_fieldnames = [
        "pair_id",
        "bacteria",
        "phage",
        "cv_group",
        "split_holdout",
        "split_cv5_fold",
        "label_hard_any_lysis",
        "prediction_context",
        "lightgbm_probability",
        "logreg_probability",
    ]
    write_csv(
        args.output_dir / "tg01_pair_predictions.csv",
        pair_prediction_fieldnames,
        project_rows_to_fields(pair_prediction_rows, pair_prediction_fieldnames),
    )
    write_csv(
        args.output_dir / "tg01_holdout_top3_rankings.csv",
        [
            "model_label",
            "bacteria",
            "phage",
            "pair_id",
            "rank",
            "predicted_probability",
            "label_hard_any_lysis",
        ],
        top3_ranking_rows,
    )

    logger.info("TG01 completed.")
    logger.info("- LightGBM best CV ROC-AUC: %s", best_lightgbm["summary"]["mean_roc_auc"])
    logger.info("- LightGBM holdout ROC-AUC: %s", lightgbm_holdout_metrics["roc_auc"])
    logger.info("- LightGBM holdout top-3 hit rate: %s", lightgbm_holdout_top3["top3_hit_rate_all_strains"])
    logger.info("- Logistic holdout ROC-AUC: %s", logreg_holdout_metrics["roc_auc"])
    logger.info("- Output directory: %s", args.output_dir)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
