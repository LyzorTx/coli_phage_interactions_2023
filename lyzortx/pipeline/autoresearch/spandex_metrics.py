"""SPANDEX evaluation metrics: nDCG (graded) and mAP (binary) per bacterium.

Replaces top-3 hit rate with ranking metrics that:
- Use graded MLC 0-4 relevance for nDCG
- Support partial ground truth (score only observed pairs per bacterium)
- Handle mixed-source data (different relevance granularity)
"""

from __future__ import annotations

import logging
from typing import Optional, Sequence

import numpy as np
from sklearn.metrics import average_precision_score, brier_score_loss, ndcg_score, roc_auc_score

LOGGER = logging.getLogger(__name__)


def _is_nan(value: object) -> bool:
    """Check if a value is NaN (works for float, np.float64, etc.)."""
    try:
        return value != value  # NaN != NaN is True
    except (TypeError, ValueError):
        return False


def compute_per_bacterium_ndcg(
    rows: Sequence[dict[str, object]],
    *,
    relevance_key: str = "mlc_score",
    probability_key: str = "predicted_probability",
) -> Optional[float]:
    """Compute mean nDCG across bacteria using graded relevance.

    For each bacterium, ranks phages by predicted probability and scores the
    ranking against graded MLC relevance (0-4). Pairs with None relevance
    are excluded (partial ground truth support).

    Returns None if no bacteria have at least one positive.
    """
    by_bacterium: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_bacterium.setdefault(str(row["bacteria"]), []).append(row)

    ndcg_values = []
    for bacteria, bact_rows in by_bacterium.items():
        observed = [r for r in bact_rows if r.get(relevance_key) is not None and not _is_nan(r[relevance_key])]
        if len(observed) < 2:
            continue
        y_rel = np.array([float(r[relevance_key]) for r in observed])
        if y_rel.max() == 0:
            continue  # no positives — nDCG undefined
        y_pred = np.array([float(r[probability_key]) for r in observed])
        ndcg_values.append(float(ndcg_score([y_rel], [y_pred])))

    if not ndcg_values:
        return None
    return float(np.mean(ndcg_values))


def compute_per_bacterium_map(
    rows: Sequence[dict[str, object]],
    *,
    label_key: str = "label_binary",
    probability_key: str = "predicted_probability",
) -> Optional[float]:
    """Compute mean Average Precision across bacteria (binary relevance).

    For each bacterium, computes AP from the ranked list of phages.
    Pairs with None label are excluded (partial ground truth support).

    Returns None if no bacteria have at least one positive.
    """
    by_bacterium: dict[str, list[dict[str, object]]] = {}
    for row in rows:
        by_bacterium.setdefault(str(row["bacteria"]), []).append(row)

    ap_values = []
    for bacteria, bact_rows in by_bacterium.items():
        observed = [r for r in bact_rows if r.get(label_key) is not None and not _is_nan(r[label_key])]
        if len(observed) < 2:
            continue
        y_true = np.array([int(r[label_key]) for r in observed])
        if y_true.sum() == 0 or y_true.sum() == len(y_true):
            continue  # all same class — AP undefined
        y_pred = np.array([float(r[probability_key]) for r in observed])
        ap_values.append(float(average_precision_score(y_true, y_pred)))

    if not ap_values:
        return None
    return float(np.mean(ap_values))


def compute_binary_metrics(
    rows: Sequence[dict[str, object]],
    *,
    label_key: str = "label_binary",
    probability_key: str = "predicted_probability",
) -> dict[str, Optional[float]]:
    """Compute AUC and Brier score (pair-level, not per-bacterium)."""
    observed = [r for r in rows if r.get(label_key) is not None and not _is_nan(r[label_key])]
    if len(observed) < 2:
        return {"roc_auc": None, "brier_score": None}
    y_true = np.array([int(r[label_key]) for r in observed])
    y_pred = np.array([float(r[probability_key]) for r in observed])
    auc = float(roc_auc_score(y_true, y_pred)) if len(np.unique(y_true)) > 1 else None
    brier = float(brier_score_loss(y_true, y_pred))
    return {"roc_auc": auc, "brier_score": brier}


SPANDEX_METRIC_NAMES = ("holdout_ndcg", "holdout_map", "holdout_roc_auc", "holdout_brier_score")


def evaluate_holdout_rows(
    rows: Sequence[dict[str, object]],
) -> dict[str, Optional[float]]:
    """Compute all SPANDEX metrics on a set of holdout rows.

    Each row must have: bacteria, predicted_probability, mlc_score (graded 0-4),
    label_binary (0 or 1). Either can be None for unobserved pairs.
    """
    ndcg = compute_per_bacterium_ndcg(rows, relevance_key="mlc_score", probability_key="predicted_probability")
    mAP = compute_per_bacterium_map(rows, label_key="label_binary", probability_key="predicted_probability")
    binary = compute_binary_metrics(rows, label_key="label_binary", probability_key="predicted_probability")
    return {
        "holdout_ndcg": ndcg,
        "holdout_map": mAP,
        "holdout_roc_auc": binary["roc_auc"],
        "holdout_brier_score": binary["brier_score"],
    }
