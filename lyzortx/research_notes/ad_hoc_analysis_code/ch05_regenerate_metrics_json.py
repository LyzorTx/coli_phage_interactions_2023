#!/usr/bin/env python3
"""CH05 post-hoc: regenerate metrics JSONs from existing predictions CSVs.

The initial CH05 artifacts were written before `_ci_to_dict` was updated to
expose `bootstrap_samples_requested` / `bootstrap_samples_used`. Rather than
rerun CH05 end-to-end (~3 hr), this script reads the existing predictions,
reruns the bootstrap from scratch, and rewrites the three JSONs under HEAD
code so they carry the new fields.

Outputs (in place):
  - ch05_bacteria_axis_metrics.json
  - ch05_phage_axis_metrics.json
  - ch05_combined_summary.json
"""

from __future__ import annotations

import json
import logging
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch.candidate_replay import safe_round
from lyzortx.pipeline.autoresearch.ch04_eval import BOOTSTRAP_RANDOM_STATE, BOOTSTRAP_SAMPLES
from lyzortx.pipeline.autoresearch.ch05_eval import (
    CH04_AGGREGATE_METRICS_PATH,
    DEFAULT_OUTPUT_DIR,
    SOURCE_BASEL,
    SOURCE_GUELIN,
    _bootstrap_by_unit,
    _ci_to_dict,
)

LOGGER = logging.getLogger(__name__)


def main() -> None:
    setup_logging()
    output_dir = Path(DEFAULT_OUTPUT_DIR)

    bacteria_pairs = pd.read_csv(output_dir / "ch05_bacteria_axis_predictions.csv")
    phage_pairs = pd.read_csv(output_dir / "ch05_phage_axis_predictions.csv")

    bacteria_rows = bacteria_pairs.to_dict(orient="records")
    phage_rows = phage_pairs.to_dict(orient="records")

    LOGGER.info("Re-running bacteria-axis bootstrap on %d pairs", len(bacteria_rows))
    bacteria_cis = _bootstrap_by_unit(
        bacteria_rows,
        unit_key="bacteria",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )
    LOGGER.info("Re-running phage-axis bootstrap on %d pairs", len(phage_rows))
    phage_cis_all = _bootstrap_by_unit(
        phage_rows,
        unit_key="phage",
        bootstrap_samples=BOOTSTRAP_SAMPLES,
        bootstrap_random_state=BOOTSTRAP_RANDOM_STATE,
    )

    bacteria_summary = {
        "n_pairs": len(bacteria_pairs),
        "n_bacteria": int(bacteria_pairs["bacteria"].nunique()),
        "n_phages": int(bacteria_pairs["phage"].nunique()),
        "aggregate": {name: _ci_to_dict(ci) for name, ci in bacteria_cis.items()},
    }
    with open(output_dir / "ch05_bacteria_axis_metrics.json", "w", encoding="utf-8") as f:
        json.dump(bacteria_summary, f, indent=2)

    cross_source_rows: list[dict[str, object]] = []
    for source_label in (SOURCE_GUELIN, SOURCE_BASEL):
        subset = [r for r in phage_rows if r["source"] == source_label]
        LOGGER.info("Re-running phage-axis bootstrap on %s subset (%d pairs)", source_label, len(subset))
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
                "auc_point": safe_round(subset_cis["holdout_roc_auc"].point_estimate),
                "auc_low": subset_cis["holdout_roc_auc"].ci_low,
                "auc_high": subset_cis["holdout_roc_auc"].ci_high,
                "brier_point": safe_round(subset_cis["holdout_brier_score"].point_estimate),
                "brier_low": subset_cis["holdout_brier_score"].ci_low,
                "brier_high": subset_cis["holdout_brier_score"].ci_high,
                "bootstrap_samples_requested": subset_cis["holdout_roc_auc"].bootstrap_samples_requested,
                "bootstrap_samples_used": subset_cis["holdout_roc_auc"].bootstrap_samples_used,
            }
        )
    phage_summary = {
        "n_pairs": len(phage_pairs),
        "n_bacteria": int(phage_pairs["bacteria"].nunique()),
        "n_phages": int(phage_pairs["phage"].nunique()),
        "aggregate": {name: _ci_to_dict(ci) for name, ci in phage_cis_all.items()},
        "cross_source": cross_source_rows,
    }
    with open(output_dir / "ch05_phage_axis_metrics.json", "w", encoding="utf-8") as f:
        json.dump(phage_summary, f, indent=2)

    phage_axis_gap = safe_round(
        bacteria_cis["holdout_roc_auc"].point_estimate - phage_cis_all["holdout_roc_auc"].point_estimate
    )
    guelin_auc = cross_source_rows[0]["auc_point"]
    basel_auc = cross_source_rows[1]["auc_point"]
    cross_source_gap = safe_round(abs(float(guelin_auc) - float(basel_auc)))

    ch04_auc = None
    if CH04_AGGREGATE_METRICS_PATH.exists():
        with open(CH04_AGGREGATE_METRICS_PATH, encoding="utf-8") as f:
            ch04 = json.load(f)
        ch04_auc = ch04["chisel_baseline"]["holdout_roc_auc"]["point_estimate"]

    combined = {
        "task_id": "CH05",
        "scorecard": "AUC + Brier (nDCG/mAP/top-k retired)",
        "per_phage_blending": "retired (see per-phage-retired-under-chisel)",
        "bacteria_axis": bacteria_summary,
        "phage_axis": phage_summary,
        "phage_axis_generalization_gap_auc": phage_axis_gap,
        "cross_source_auc_delta": cross_source_gap,
        "ch04_baseline_auc": ch04_auc,
        "regenerated_at": datetime.now(timezone.utc).isoformat(),
        "regenerated_from": "existing predictions CSVs (model not rerun)",
    }
    with open(output_dir / "ch05_combined_summary.json", "w", encoding="utf-8") as f:
        json.dump(combined, f, indent=2)
    LOGGER.info("Regenerated 3 metric JSONs under %s", output_dir)


if __name__ == "__main__":
    main()
