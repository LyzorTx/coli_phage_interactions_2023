#!/usr/bin/env python3
"""AR01: lock the AUTORESEARCH corpus, label policy, and sealed split contract."""

from __future__ import annotations

import argparse
import logging
from collections import Counter, defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Tuple

from lyzortx.log_config import setup_logging
from lyzortx.pipeline.autoresearch import runtime_contract
from lyzortx.pipeline.deployment_paired_features.download_picard_assemblies import (
    DEFAULT_ASSEMBLY_DIR,
    FASTA_SUFFIXES,
    download_picard_assemblies,
    list_fasta_files,
)
from lyzortx.pipeline.steel_thread_v0.io.load_inputs import iter_raw_interactions
from lyzortx.pipeline.steel_thread_v0.io.write_outputs import ensure_directory, write_csv, write_json
from lyzortx.pipeline.steel_thread_v0.steps.st02_build_pair_table import BORDERLINE_NOISE_WEIGHT
from lyzortx.pipeline.track_a.steps.build_track_a_foundation import (
    EXCLUDED_LOG_DILUTIONS,
    LabelPolicyV1,
    compute_label_v1,
)

LOGGER = logging.getLogger(__name__)

DEFAULT_RAW_INTERACTIONS_PATH = Path("data/interactions/raw/raw_interactions.csv")
DEFAULT_PHAGE_FASTA_DIR = Path("data/genomics/phages/FNA")
DEFAULT_OUTPUT_DIR = runtime_contract.DEFAULT_OUTPUT_ROOT

PAIR_TABLE_FILENAME = "ar01_canonical_pair_table_v1.csv"
CONTRACT_MANIFEST_FILENAME = "ar01_split_benchmark_manifest_v1.json"
INPUT_CHECKSUMS_FILENAME = "ar01_input_checksums_v1.json"
LABEL_POLICY_FILENAME = "ar01_label_policy_v1.json"

LABEL_POLICY_ID = "autoresearch_label_policy_v1"
SPLIT_CONTRACT_ID = "autoresearch_bacteria_disjoint_split_v1"
PAIR_TABLE_ID = "autoresearch_pair_table_v1"
INPUT_CONTRACT_ID = "autoresearch_raw_corpus_contract_v1"

TRAIN_SPLIT = runtime_contract.TRAIN_SPLIT
INNER_VAL_SPLIT = runtime_contract.INNER_VAL_SPLIT
HOLDOUT_SPLIT = runtime_contract.HOLDOUT_SPLIT
SPLIT_ORDER = runtime_contract.SPLIT_ORDER

DEFAULT_HOLDOUT_FRACTION = 0.2
DEFAULT_INNER_VAL_FRACTION = 0.2
DEFAULT_SPLIT_SALT = "autoresearch_ar01_bacteria_split_v1"

EXCLUSION_REASON_UNRESOLVED_LABEL = "unresolved_label"
EXCLUSION_REASON_MISSING_HOST_FASTA = "missing_host_fasta"
EXCLUSION_REASON_MISSING_PHAGE_FASTA = "missing_phage_fasta"

COMPARATOR_BENCHMARK = {
    "artifact_id": "track_g_clean_v1_locked_benchmark",
    "benchmark_summary_path": "lyzortx/generated_outputs/track_g/tg02_gbm_calibration/tg02_benchmark_summary.json",
    "feature_lock_path": (
        "lyzortx/generated_outputs/track_g/tg05_feature_subset_sweep/tg05_locked_v1_feature_config.json"
    ),
    "model_summary_path": "lyzortx/generated_outputs/track_g/tg01_v1_binary_classifier/tg01_model_summary.json",
    "evaluation_protocol_id": "steel_thread_v0_st03_split_v1",
    "locked_feature_blocks": ["defense", "phage_genomic"],
    "selection_source": "TG05/TG09 clean v1 lock carried forward as the current production-intent comparator",
}
COMPARATOR_BENCHMARK_PATH_KEYS = ("benchmark_summary_path", "feature_lock_path", "model_summary_path")


def parse_args(argv: Optional[Sequence[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-interactions-path",
        type=Path,
        default=DEFAULT_RAW_INTERACTIONS_PATH,
        help="Semicolon-delimited raw interaction table.",
    )
    parser.add_argument(
        "--host-assembly-dir",
        type=Path,
        default=DEFAULT_ASSEMBLY_DIR,
        help="Directory containing Picard host FASTAs.",
    )
    parser.add_argument(
        "--phage-fasta-dir",
        type=Path,
        default=DEFAULT_PHAGE_FASTA_DIR,
        help="Directory containing phage FASTA files.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=DEFAULT_OUTPUT_DIR,
        help="Directory for AUTORESEARCH generated outputs.",
    )
    parser.add_argument(
        "--holdout-fraction",
        type=float,
        default=DEFAULT_HOLDOUT_FRACTION,
        help="Fraction of bacteria routed to sealed holdout.",
    )
    parser.add_argument(
        "--inner-val-fraction",
        type=float,
        default=DEFAULT_INNER_VAL_FRACTION,
        help="Fraction of bacteria routed to inner validation.",
    )
    parser.add_argument(
        "--split-salt",
        default=DEFAULT_SPLIT_SALT,
        help="Deterministic salt for bacteria split assignment.",
    )
    parser.add_argument(
        "--skip-host-assembly-resolution",
        action="store_true",
        help="Skip download_picard_assemblies() and trust the provided host assembly directory as-is.",
    )
    parser.add_argument(
        "--skip-comparator-lock",
        action="store_true",
        help=(
            "Skip comparator artifact checksum validation.  Use in CI/RunPod where "
            "Track G generated outputs are not available.  The comparator reference "
            "(paths + artifact ID) is still recorded; checksums are deferred to AR09 replay."
        ),
    )
    return parser.parse_args(argv)


def sha256_file(path: Path) -> str:
    return runtime_contract.sha256_file(path)


def normalized_hash_01(text: str) -> float:
    return runtime_contract.normalized_hash_01(text)


def sha256_strings(values: Iterable[str]) -> str:
    return runtime_contract.sha256_strings(values)


def discover_fasta_map(root: Path) -> Dict[str, Path]:
    if not root.exists():
        raise FileNotFoundError(f"FASTA directory not found: {root}")
    fasta_paths = sorted(path for path in root.rglob("*") if path.is_file() and path.suffix.lower() in FASTA_SUFFIXES)
    if not fasta_paths:
        raise ValueError(f"No FASTA files found under {root}")

    out: Dict[str, Path] = {}
    for path in fasta_paths:
        stem = path.stem
        if stem in out:
            raise ValueError(f"Duplicate FASTA stem '{stem}' detected under {root}: {out[stem]} and {path}")
        out[stem] = path
    return out


def build_fasta_checksums(paths_by_id: Mapping[str, Path], *, root: Path) -> Dict[str, object]:
    items = []
    for key, path in sorted(paths_by_id.items()):
        items.append(
            {
                "id": key,
                "relative_path": str(path.relative_to(root)),
                "sha256": sha256_file(path),
            }
        )
    return {
        "root": str(root),
        "file_count": len(items),
        "aggregate_sha256": sha256_strings(f"{item['relative_path']}:{item['sha256']}" for item in items),
        "files": items,
    }


def stable_float_string(value: float) -> str:
    return f"{value:.6f}".rstrip("0").rstrip(".")


def build_dilution_counts(rows: Sequence[Dict[str, str]]) -> Dict[int, Counter[str]]:
    dilution_counts: Dict[int, Counter[str]] = defaultdict(Counter)
    for row in rows:
        dilution_counts[int(row["log_dilution"])][row["score"]] += 1
    return dict(dilution_counts)


def has_repeated_within_dilution_lysis(dilution_counts: Mapping[int, Counter[str]]) -> bool:
    return any(counts["1"] >= 2 for counts in dilution_counts.values())


def compute_autoresearch_training_weight(
    *,
    label_any_lysis: str,
    dilution_counts: Mapping[int, Counter[str]],
) -> float:
    if label_any_lysis != "1":
        return 1.0
    if has_repeated_within_dilution_lysis(dilution_counts):
        return 1.0
    return float(BORDERLINE_NOISE_WEIGHT)


def assign_bacteria_splits(
    bacteria_ids: Sequence[str],
    *,
    holdout_fraction: float,
    inner_val_fraction: float,
    split_salt: str,
) -> Dict[str, str]:
    if not 0 < holdout_fraction < 1:
        raise ValueError("holdout_fraction must be between 0 and 1.")
    if not 0 < inner_val_fraction < 1:
        raise ValueError("inner_val_fraction must be between 0 and 1.")
    if holdout_fraction + inner_val_fraction >= 1:
        raise ValueError("holdout_fraction + inner_val_fraction must be < 1.")
    if len(bacteria_ids) < 3:
        raise ValueError("Need at least 3 bacteria to predeclare train/inner_val/holdout splits.")

    scored = sorted((normalized_hash_01(f"{split_salt}|{bacteria}"), bacteria) for bacteria in sorted(bacteria_ids))
    total = len(scored)
    n_holdout = max(1, int(round(total * holdout_fraction)))
    n_inner_val = max(1, int(round(total * inner_val_fraction)))
    if n_holdout + n_inner_val >= total:
        raise ValueError("Split fractions leave no bacteria for train.")

    assignments: Dict[str, str] = {}
    for _, bacteria in scored[:n_holdout]:
        assignments[bacteria] = HOLDOUT_SPLIT
    for _, bacteria in scored[n_holdout : n_holdout + n_inner_val]:
        assignments[bacteria] = INNER_VAL_SPLIT
    for _, bacteria in scored[n_holdout + n_inner_val :]:
        assignments[bacteria] = TRAIN_SPLIT

    for split_name in SPLIT_ORDER:
        if split_name not in assignments.values():
            raise ValueError(f"Split '{split_name}' ended up empty.")
    return assignments


def to_str_column(row: Dict[str, object], key: str) -> str:
    value = row[key]
    return "" if value is None else str(value)


def build_locked_comparator_benchmark(benchmark: Mapping[str, object]) -> Dict[str, object]:
    resolved_benchmark = dict(benchmark)
    artifact_checksums: Dict[str, str] = {}
    for path_key in COMPARATOR_BENCHMARK_PATH_KEYS:
        artifact_path = Path(str(benchmark[path_key]))
        if not artifact_path.exists():
            raise FileNotFoundError(f"Locked comparator artifact not found: {artifact_path}")
        if not artifact_path.is_file():
            raise FileNotFoundError(f"Locked comparator artifact is not a file: {artifact_path}")
        artifact_checksums[path_key] = sha256_file(artifact_path)
    resolved_benchmark["artifact_checksums"] = artifact_checksums
    return resolved_benchmark


def aggregate_raw_pairs(
    raw_interactions_path: Path,
) -> Tuple[List[Dict[str, str]], Dict[Tuple[str, str], List[Dict[str, str]]]]:
    # SX05: drop log_dilution rows that the paper excludes from MLC (currently the unreplicated
    # 5 x 10^4 pfu/ml observation). See EXCLUDED_LOG_DILUTIONS rationale in build_track_a_foundation.py.
    raw_rows = [
        row
        for row in iter_raw_interactions(raw_interactions_path)
        if int(row["log_dilution"]) not in EXCLUDED_LOG_DILUTIONS
    ]
    if not raw_rows:
        raise ValueError(f"No rows found in {raw_interactions_path}")
    by_pair: Dict[Tuple[str, str], List[Dict[str, str]]] = defaultdict(list)
    for row in raw_rows:
        by_pair[(row["bacteria"], row["phage"])].append(row)
    return raw_rows, by_pair


def resolve_host_fasta_map(
    *,
    raw_interactions_path: Path,
    host_assembly_dir: Path,
    skip_host_assembly_resolution: bool,
) -> Dict[str, Path]:
    if skip_host_assembly_resolution:
        LOGGER.info("Skipping Picard assembly resolution; using existing FASTAs in %s", host_assembly_dir)
    else:
        LOGGER.info("Resolving Picard host assemblies with download_picard_assemblies()")
        download_picard_assemblies(assembly_dir=host_assembly_dir, raw_interactions_path=raw_interactions_path)
    fasta_paths = list_fasta_files(host_assembly_dir)
    if not fasta_paths:
        raise ValueError(f"No host FASTA files found in {host_assembly_dir}")
    return discover_fasta_map(host_assembly_dir)


def build_pair_rows(
    *,
    pair_rows: Mapping[Tuple[str, str], Sequence[Dict[str, str]]],
    host_fasta_map: Mapping[str, Path],
    phage_fasta_map: Mapping[str, Path],
    bacteria_split_map: Mapping[str, str],
) -> List[Dict[str, object]]:
    policy = LabelPolicyV1()
    out_rows: List[Dict[str, object]] = []

    for bacteria, phage in sorted(pair_rows):
        rows = pair_rows[(bacteria, phage)]
        total_obs = len(rows)
        score_counts = Counter(row["score"] for row in rows)
        dilution_counts = build_dilution_counts(rows)
        label_v1 = compute_label_v1(
            score_1_count=score_counts["1"],
            score_0_count=score_counts["0"],
            score_n_count=score_counts["n"],
            total_obs=total_obs,
            dilution_counts=dilution_counts,
            policy=policy,
        )

        host_fasta = host_fasta_map.get(bacteria)
        phage_fasta = phage_fasta_map.get(phage)
        exclusion_reasons: List[str] = []
        if label_v1["any_lysis"] == "":
            exclusion_reasons.append(EXCLUSION_REASON_UNRESOLVED_LABEL)
        if host_fasta is None:
            exclusion_reasons.append(EXCLUSION_REASON_MISSING_HOST_FASTA)
        if phage_fasta is None:
            exclusion_reasons.append(EXCLUSION_REASON_MISSING_PHAGE_FASTA)

        training_weight = compute_autoresearch_training_weight(
            label_any_lysis=str(label_v1["any_lysis"]),
            dilution_counts=dilution_counts,
        )
        retained_for_search = int(not exclusion_reasons)
        split_name = bacteria_split_map[bacteria]

        row = {
            "pair_id": f"{bacteria}__{phage}",
            "bacteria": bacteria,
            "phage": phage,
            "split": split_name,
            "label_any_lysis": label_v1["any_lysis"],
            "label_reason": label_v1["label_reason"],
            "training_weight_v3": stable_float_string(training_weight),
            "label_read_only": 1,
            "retained_for_autoresearch": retained_for_search,
            "exclusion_reasons": "|".join(exclusion_reasons),
            "host_fasta_path": "" if host_fasta is None else str(host_fasta),
            "phage_fasta_path": "" if phage_fasta is None else str(phage_fasta),
            "obs_total": total_obs,
            "obs_score_1_count": score_counts["1"],
            "obs_score_0_count": score_counts["0"],
            "obs_score_n_count": score_counts["n"],
            "obs_interpretable_count": score_counts["1"] + score_counts["0"],
            "obs_uncertainty_flags": "|".join(label_v1["uncertainty_flags"]),
            "obs_uncertainty_flag_count": len(label_v1["uncertainty_flags"]),
            "obs_best_lysis_dilution": label_v1["best_lysis_dilution"],
            "obs_dilution_potency": label_v1["dilution_potency"],
            "obs_potency_support": label_v1["potency_support"],
            "obs_has_repeated_within_dilution_lysis": int(has_repeated_within_dilution_lysis(dilution_counts)),
        }
        out_rows.append(row)

    for row in out_rows:
        retained = row["retained_for_autoresearch"] == 1
        if retained and (row["host_fasta_path"] == "" or row["phage_fasta_path"] == ""):
            raise ValueError(
                "Retained AUTORESEARCH pair is missing a FASTA: "
                f"{row['bacteria']} / {row['phage']} (split={row['split']})"
            )
    return out_rows


def build_label_policy_manifest() -> Dict[str, object]:
    policy = LabelPolicyV1()
    return {
        "policy_id": LABEL_POLICY_ID,
        "policy_version": "v1",
        "source_policy": "Track A label_set_v1 semantics plus AUTORESEARCH raw-contract training_weight_v3",
        "labels_read_only": True,
        "pair_key": ["bacteria", "phage"],
        "pair_table_id": PAIR_TABLE_ID,
        "track_a_any_lysis_rules": {
            "positive": "label_any_lysis=1 if any raw observation has score='1'",
            "negative": "label_any_lysis=0 if no score='1' and score_0_count >= 5",
            "unresolved": "label_any_lysis is blank otherwise",
        },
        "score_n_handling": {
            "training_label_effect": (
                "score='n' never creates a positive and never counts toward the hard-negative threshold"
            ),
            "uncertainty_effect": (
                "score='n' contributes to has_uninterpretable and high_uninterpretable_fraction flags"
            ),
            "resolution_effect": (
                "pairs with no score='1' and fewer than 5 score='0' observations remain unresolved and are excluded"
            ),
        },
        "training_weight_v3": {
            "default_weight": 1.0,
            "downweighted_weight": float(BORDERLINE_NOISE_WEIGHT),
            "downweight_rule": (
                "If label_any_lysis=1 and no dilution has repeated lysis support (>=2 score='1' observations at the "
                "same dilution), set training_weight_v3=0.1; otherwise 1.0."
            ),
            "rationale": (
                "This is the raw-only AUTORESEARCH equivalent of TA11's borderline noise-positive downweighting. "
                "It preserves the intent of training_weight_v3 while removing the hidden interaction_matrix.csv "
                "dependency from the AUTORESEARCH input contract."
            ),
        },
        "thresholds": {
            "expected_observations_per_pair": policy.expected_observations_per_pair,
            "min_interpretable_obs_for_hard_negative": policy.min_interpretable_obs_for_hard_negative,
            "high_uninterpretable_fraction_threshold": policy.high_uninterpretable_fraction_threshold,
        },
    }


def build_split_summary(pair_table_rows: Sequence[Dict[str, object]]) -> Dict[str, object]:
    split_bacteria: Dict[str, set[str]] = {split: set() for split in SPLIT_ORDER}
    split_phages: Dict[str, set[str]] = {split: set() for split in SPLIT_ORDER}
    split_rows: Dict[str, List[str]] = {split: [] for split in SPLIT_ORDER}
    split_retained_rows: Dict[str, List[str]] = {split: [] for split in SPLIT_ORDER}

    for row in pair_table_rows:
        split = to_str_column(row, "split")
        bacteria = to_str_column(row, "bacteria")
        phage = to_str_column(row, "phage")
        pair_id = to_str_column(row, "pair_id")

        split_bacteria[split].add(bacteria)
        split_phages[split].add(phage)
        split_rows[split].append(pair_id)
        if to_str_column(row, "retained_for_autoresearch") == "1":
            split_retained_rows[split].append(pair_id)

    bacteria_overlap_counts = {}
    for left in SPLIT_ORDER:
        for right in SPLIT_ORDER:
            if left >= right:
                continue
            bacteria_overlap_counts[f"{left}__{right}"] = len(split_bacteria[left] & split_bacteria[right])

    return {
        "split_contract_id": SPLIT_CONTRACT_ID,
        "assignment_unit": "bacteria",
        "sealed_split": HOLDOUT_SPLIT,
        "row_counts": {split: len(split_rows[split]) for split in SPLIT_ORDER},
        "retained_row_counts": {split: len(split_retained_rows[split]) for split in SPLIT_ORDER},
        "bacteria_counts": {split: len(split_bacteria[split]) for split in SPLIT_ORDER},
        "phage_counts": {split: len(split_phages[split]) for split in SPLIT_ORDER},
        "bacteria_ids": {split: sorted(split_bacteria[split]) for split in SPLIT_ORDER},
        "split_hashes": {
            split: {
                "bacteria_sha256": sha256_strings(sorted(split_bacteria[split])),
                "pair_ids_sha256": sha256_strings(sorted(split_rows[split])),
                "retained_pair_ids_sha256": sha256_strings(sorted(split_retained_rows[split])),
            }
            for split in SPLIT_ORDER
        },
        "leakage_checks": {
            "bacteria_overlap_counts": bacteria_overlap_counts,
            "max_overlap_count": max(bacteria_overlap_counts.values(), default=0),
        },
    }


def validate_split_contract(split_summary: Mapping[str, object]) -> None:
    leakage_checks = split_summary["leakage_checks"]
    if not isinstance(leakage_checks, Mapping) or leakage_checks["max_overlap_count"] != 0:
        raise ValueError("Bacterium overlap detected across AUTORESEARCH splits.")

    retained_row_counts = split_summary["retained_row_counts"]
    if not isinstance(retained_row_counts, Mapping):
        raise ValueError("AUTORESEARCH split summary is missing retained_row_counts.")

    empty_retained_splits = [split for split in SPLIT_ORDER if retained_row_counts[split] == 0]
    if empty_retained_splits:
        raise ValueError("AUTORESEARCH retained split(s) empty after exclusions: " + ", ".join(empty_retained_splits))


def build_contract_manifest(
    *,
    raw_interactions_path: Path,
    host_assembly_dir: Path,
    phage_fasta_dir: Path,
    output_dir: Path,
    pair_table_rows: Sequence[Dict[str, object]],
    pair_table_path: Path,
    input_checksums_path: Path,
    label_policy_path: Path,
    split_salt: str,
    holdout_fraction: float,
    inner_val_fraction: float,
    skip_comparator_lock: bool = False,
) -> Dict[str, object]:
    split_summary = build_split_summary(pair_table_rows)
    validate_split_contract(split_summary)

    return {
        "task_id": "AR01",
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "input_contract_id": INPUT_CONTRACT_ID,
        "pair_table_id": PAIR_TABLE_ID,
        "input_contract": {
            "raw_interactions_path": str(raw_interactions_path),
            "host_assemblies_resolver": "lyzortx/pipeline/deployment_paired_features/download_picard_assemblies.py",
            "host_assembly_dir": str(host_assembly_dir),
            "phage_fasta_dir": str(phage_fasta_dir),
            "notes": (
                "AUTORESEARCH training inputs are reproducible from raw interactions plus resolved host/phage FASTAs. "
                "No DEPLOY feature CSV is an AUTORESEARCH source-of-truth input."
            ),
        },
        "pair_table": {
            "path": str(pair_table_path),
            "sha256": sha256_file(pair_table_path),
            "row_count": len(pair_table_rows),
        },
        "label_policy_manifest_path": str(label_policy_path),
        "input_checksums_manifest_path": str(input_checksums_path),
        "split_contract": {
            "split_salt": split_salt,
            "holdout_fraction": holdout_fraction,
            "inner_val_fraction": inner_val_fraction,
            "train_fraction": 1.0 - holdout_fraction - inner_val_fraction,
            **split_summary,
        },
        "current_locked_comparator_benchmark": (
            dict(COMPARATOR_BENCHMARK, artifact_checksums=None, artifact_checksums_deferred=True)
            if skip_comparator_lock
            else build_locked_comparator_benchmark(COMPARATOR_BENCHMARK)
        ),
        "output_dir": str(output_dir),
    }


def main(argv: Optional[Sequence[str]] = None) -> int:
    setup_logging()
    args = parse_args(argv)
    LOGGER.info("AR01 starting: build AUTORESEARCH contract")
    ensure_directory(args.output_dir)

    raw_rows, pair_rows = aggregate_raw_pairs(args.raw_interactions_path)
    bacteria_ids = sorted({row["bacteria"] for row in raw_rows})

    host_fasta_map = resolve_host_fasta_map(
        raw_interactions_path=args.raw_interactions_path,
        host_assembly_dir=args.host_assembly_dir,
        skip_host_assembly_resolution=args.skip_host_assembly_resolution,
    )
    phage_fasta_map = discover_fasta_map(args.phage_fasta_dir)
    bacteria_split_map = assign_bacteria_splits(
        bacteria_ids,
        holdout_fraction=args.holdout_fraction,
        inner_val_fraction=args.inner_val_fraction,
        split_salt=args.split_salt,
    )

    pair_table_rows = build_pair_rows(
        pair_rows=pair_rows,
        host_fasta_map=host_fasta_map,
        phage_fasta_map=phage_fasta_map,
        bacteria_split_map=bacteria_split_map,
    )
    pair_table_path = args.output_dir / PAIR_TABLE_FILENAME
    write_csv(pair_table_path, fieldnames=list(pair_table_rows[0].keys()), rows=pair_table_rows)

    input_checksums = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "raw_interactions": {
            "path": str(args.raw_interactions_path),
            "sha256": sha256_file(args.raw_interactions_path),
        },
        "host_fastas": build_fasta_checksums(host_fasta_map, root=args.host_assembly_dir),
        "phage_fastas": build_fasta_checksums(phage_fasta_map, root=args.phage_fasta_dir),
    }
    input_checksums_path = args.output_dir / INPUT_CHECKSUMS_FILENAME
    write_json(input_checksums_path, input_checksums)

    label_policy_manifest = build_label_policy_manifest()
    label_policy_manifest_path = args.output_dir / LABEL_POLICY_FILENAME
    write_json(label_policy_manifest_path, label_policy_manifest)

    contract_manifest = build_contract_manifest(
        raw_interactions_path=args.raw_interactions_path,
        host_assembly_dir=args.host_assembly_dir,
        phage_fasta_dir=args.phage_fasta_dir,
        output_dir=args.output_dir,
        pair_table_rows=pair_table_rows,
        pair_table_path=pair_table_path,
        input_checksums_path=input_checksums_path,
        label_policy_path=label_policy_manifest_path,
        split_salt=args.split_salt,
        holdout_fraction=args.holdout_fraction,
        inner_val_fraction=args.inner_val_fraction,
        skip_comparator_lock=args.skip_comparator_lock,
    )
    write_json(args.output_dir / CONTRACT_MANIFEST_FILENAME, contract_manifest)

    LOGGER.info("AR01 completed: wrote %s", pair_table_path)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
