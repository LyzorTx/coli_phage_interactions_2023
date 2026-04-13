#!/usr/bin/env python3
"""Build Track A artifacts: IDs, integrity checks, cohorts, labels, and plaque-QC queue."""

from __future__ import annotations

import argparse
import csv
import hashlib
import re
from collections import Counter, defaultdict
from dataclasses import asdict, dataclass
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
from typing import Dict, Iterable, List, Mapping, Optional, Sequence, Set, Tuple

from lyzortx.pipeline.steel_thread_v0.io.write_outputs import ensure_directory, write_csv, write_json

RAW_REQUIRED_COLUMNS: Tuple[str, ...] = (
    "bacteria",
    "bacteria_index",
    "phage",
    "image",
    "replicate",
    "plate",
    "log_dilution",
    "X",
    "Y",
    "score",
)

HOST_REQUIRED_COLUMNS: Tuple[str, ...] = (
    "bacteria",
    "Pathotype",
    "Clermont_Phylo",
    "Origin",
    "LPS_type",
    "O-type",
    "H-type",
    "Collection",
    "ABC_serotype",
)

PHAGE_REQUIRED_COLUMNS: Tuple[str, ...] = (
    "phage",
    "Morphotype",
    "Family",
    "Genus",
    "Species",
    "Genome_size",
    "Phage_host",
    "Phage_host_phylo",
)

CV_REQUIRED_COLUMNS: Tuple[str, ...] = (
    "bacteria",
    "group",
)

ALLOWED_RAW_SCORES: Set[str] = {"0", "1", "n"}

# Per Gaborieau 2024 Methods ("Evaluating phage-bacteria interaction outcomes by plaque assay experiments"):
# "The outcome of interaction at 5 x 10^4 pfu/ml was not taken into account in the calculation of the MLC score
# because it was not verified by a replicate." SX05 aligns our pipeline with that protocol by dropping
# log_dilution=-4 (5 x 10^4 pfu/ml) from the MLC weight map and filtering those rows out of every downstream
# scoring path. The paper's MLC=4 is a morphological distinction ("entire lysis of the bacterial lawn at
# 5 x 10^6") that our binary 0/1/n raw data physically cannot capture, so the fix also collapses our MLC range
# from {0, 1, 2, 3, 4} to {0, 1, 2, 3}.
EXCLUDED_LOG_DILUTIONS: frozenset[int] = frozenset({-4})

DILUTION_WEIGHT_MAP: Dict[int, int] = {
    0: 1,
    -1: 2,
    -2: 3,
}

DILUTION_POTENCY_LABEL_MAP: Dict[int, str] = {
    0: "low",
    -1: "moderate",
    -2: "high",
}

IMAGE_EXTENSIONS: Tuple[str, ...] = (
    ".jpg",
    ".jpeg",
    ".png",
    ".tif",
    ".tiff",
)


@dataclass(frozen=True)
class LabelPolicyV1:
    # 8 = R1/R2/R3 at log_dilution=0 + R2/R3 at -1 + R1/R2/R3 at -2
    # (SX05 drops log_dilution=-4, the unreplicated 5 x 10^4 pfu/ml observation.)
    expected_observations_per_pair: int = 8
    min_interpretable_obs_for_hard_negative: int = 5
    high_uninterpretable_fraction_threshold: float = 0.25


@dataclass(frozen=True)
class LabelPolicyV2:
    expected_observations_per_pair: int = 8
    min_positive_obs: int = 2
    min_positive_fraction_interpretable: float = 0.4
    min_negative_obs: int = 7
    max_uninterpretable_obs: int = 1
    high_uninterpretable_fraction_threshold: float = 0.25


def parse_args(argv: Optional[List[str]] = None) -> argparse.Namespace:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--raw-interactions-path",
        type=Path,
        default=Path("data/interactions/raw/raw_interactions.csv"),
        help="Input raw interactions CSV (semicolon-delimited).",
    )
    parser.add_argument(
        "--interaction-matrix-path",
        type=Path,
        default=Path("data/interactions/interaction_matrix.csv"),
        help="Input interaction matrix CSV (semicolon-delimited).",
    )
    parser.add_argument(
        "--host-metadata-path",
        type=Path,
        default=Path("data/genomics/bacteria/picard_collection.csv"),
        help="Input host metadata CSV (semicolon-delimited).",
    )
    parser.add_argument(
        "--phage-metadata-path",
        type=Path,
        default=Path("data/genomics/phages/guelin_collection.csv"),
        help="Input phage metadata CSV (semicolon-delimited).",
    )
    parser.add_argument(
        "--cv-groups-path",
        type=Path,
        default=Path("data/metadata/370+host_cross_validation_groups_1e-4.csv"),
        help="Input cross-validation group CSV (semicolon-delimited).",
    )
    parser.add_argument(
        "--alias-overrides-path",
        type=Path,
        default=Path("lyzortx/pipeline/track_a/config/alias_overrides.csv"),
        help="CSV containing manual alias overrides.",
    )
    parser.add_argument(
        "--image-search-roots",
        type=str,
        default="paper/raw_plaque_images,data/interactions/raw,paper/_extract",
        help="Comma-separated roots for plaque image lookup.",
    )
    parser.add_argument(
        "--output-dir",
        type=Path,
        default=Path("lyzortx/generated_outputs/track_a"),
        help="Output directory for Track A artifacts.",
    )
    parser.add_argument(
        "--expected-observations-per-pair",
        type=int,
        default=8,
        help=(
            "Expected observation count per bacteria-phage pair after excluding "
            "EXCLUDED_LOG_DILUTIONS (SX05 drops the unreplicated 5 x 10^4 row, leaving 8 per pair)."
        ),
    )
    parser.add_argument(
        "--v1-min-interpretable-neg",
        type=int,
        default=5,
        help="V1 minimum interpretable no-lysis count for hard negative.",
    )
    parser.add_argument(
        "--v1-high-uninterpretable-frac",
        type=float,
        default=0.25,
        help="V1 high-uninterpretable fraction threshold.",
    )
    parser.add_argument(
        "--v2-min-positive-obs",
        type=int,
        default=2,
        help="V2 minimum score=1 count for strict positive.",
    )
    parser.add_argument(
        "--v2-min-positive-fraction",
        type=float,
        default=0.4,
        help="V2 minimum positive fraction among interpretable observations for strict positive.",
    )
    parser.add_argument(
        "--v2-min-negative-obs",
        type=int,
        default=7,
        help="V2 minimum score=0 count for strict negative.",
    )
    parser.add_argument(
        "--v2-max-uninterpretable-obs",
        type=int,
        default=1,
        help="V2 maximum score=n count for strict labels.",
    )
    parser.add_argument(
        "--v2-high-uninterpretable-frac",
        type=float,
        default=0.25,
        help="V2 high-uninterpretable fraction threshold.",
    )
    return parser.parse_args(argv)


def filter_excluded_dilutions(raw_rows: Sequence[Dict[str, str]]) -> List[Dict[str, str]]:
    """Drop raw rows whose log_dilution is in EXCLUDED_LOG_DILUTIONS (SX05)."""
    kept: List[Dict[str, str]] = []
    for row in raw_rows:
        dilution = parse_int_or_none(row["log_dilution"])
        if dilution is None:
            raise ValueError(f"Invalid log_dilution value: {row['log_dilution']!r}")
        if dilution in EXCLUDED_LOG_DILUTIONS:
            continue
        kept.append(row)
    return kept


def read_delimited_rows(path: Path, delimiter: str = ";") -> Tuple[List[Dict[str, str]], List[str]]:
    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter=delimiter)
        if reader.fieldnames is None:
            raise ValueError(f"No header found in {path}.")
        rows = []
        for row in reader:
            rows.append({k: (v.strip() if isinstance(v, str) else "") for k, v in row.items()})
        return rows, list(reader.fieldnames)


def read_alias_overrides(path: Path) -> Dict[str, Dict[str, Dict[str, str]]]:
    out: Dict[str, Dict[str, Dict[str, str]]] = {
        "bacteria": {},
        "phage": {},
    }
    if not path.exists():
        return out

    with path.open("r", newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames is None:
            return out
        required = {"entity_type", "alias_name", "canonical_name", "reason", "apply"}
        missing = sorted(required - set(reader.fieldnames))
        if missing:
            raise ValueError(f"Missing alias override columns in {path}: {', '.join(missing)}")
        for row in reader:
            entity_type = row["entity_type"].strip().lower()
            if entity_type not in out:
                continue
            apply_raw = row["apply"].strip().lower()
            apply_flag = apply_raw in {"1", "true", "yes", "y"}
            if not apply_flag:
                continue
            alias_name = row["alias_name"].strip()
            canonical_name = row["canonical_name"].strip()
            if alias_name == "" or canonical_name == "":
                continue
            out[entity_type][alias_name] = {
                "canonical_name": canonical_name,
                "reason": row["reason"].strip(),
            }
    return out


def normalize_token(value: str) -> str:
    return re.sub(r"[^a-z0-9]+", "", value.lower().strip())


def resolve_alias(
    name: str,
    entity_type: str,
    override_map: Dict[str, Dict[str, Dict[str, str]]],
) -> Tuple[str, bool, str]:
    alias_map = override_map.get(entity_type, {})
    if name in alias_map:
        return alias_map[name]["canonical_name"], True, alias_map[name]["reason"]
    return name, False, ""


def sha256_of_file(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def parse_int_or_none(value: str) -> Optional[int]:
    if value == "":
        return None
    try:
        return int(value)
    except ValueError:
        return None


def parse_float_or_none(value: str) -> Optional[float]:
    if value == "":
        return None
    try:
        return float(value)
    except ValueError:
        return None


def build_alias_candidates(
    names_by_source: Dict[str, Set[str]],
    entity_type: str,
    override_map: Dict[str, Dict[str, Dict[str, str]]],
) -> List[Dict[str, object]]:
    all_names = sorted(set().union(*names_by_source.values()))
    source_by_name: Dict[str, Set[str]] = defaultdict(set)
    for source, names in names_by_source.items():
        for name in names:
            source_by_name[name].add(source)

    candidates: List[Dict[str, object]] = []
    seen_pairs: Set[Tuple[str, str]] = set()
    alias_names = set(override_map.get(entity_type, {}).keys())

    for idx, left in enumerate(all_names):
        left_norm = normalize_token(left)
        if left_norm == "":
            continue
        for right in all_names[idx + 1 :]:
            if left == right:
                continue
            pair = (left, right)
            if pair in seen_pairs:
                continue
            seen_pairs.add(pair)

            right_norm = normalize_token(right)
            if right_norm == "":
                continue

            confidence = 0.0
            reason = ""
            ratio = 0.0
            if left_norm == right_norm:
                confidence = 1.0
                reason = "normalized_exact"
                ratio = 1.0
            else:
                ratio = SequenceMatcher(None, left_norm, right_norm).ratio()
                min_len = min(len(left_norm), len(right_norm))
                if min_len >= 5 and (left_norm in right_norm or right_norm in left_norm):
                    confidence = 0.93
                    reason = "normalized_substring"
                elif min_len >= 6 and ratio >= 0.88:
                    confidence = ratio
                    reason = "string_similarity"

            if reason == "":
                continue

            left_sources = sorted(source_by_name[left])
            right_sources = sorted(source_by_name[right])
            if set(left_sources) == set(right_sources):
                continue

            if left in alias_names or right in alias_names:
                continue

            candidates.append(
                {
                    "entity_type": entity_type,
                    "left_name": left,
                    "right_name": right,
                    "left_sources": "|".join(left_sources),
                    "right_sources": "|".join(right_sources),
                    "match_reason": reason,
                    "similarity_ratio": round(ratio, 6),
                    "confidence_score": round(confidence, 6),
                    "review_status": "pending",
                }
            )

    candidates.sort(
        key=lambda row: (
            -float(row["confidence_score"]),
            -float(row["similarity_ratio"]),
            str(row["left_name"]),
            str(row["right_name"]),
        )
    )
    return candidates


def potency_label_from_dilution(dilution: Optional[int]) -> str:
    if dilution is None or dilution in EXCLUDED_LOG_DILUTIONS:
        return "none"
    return DILUTION_POTENCY_LABEL_MAP.get(dilution, f"dilution_{dilution}")


def potency_rank_from_dilution(dilution: Optional[int]) -> int:
    if dilution is None or dilution in EXCLUDED_LOG_DILUTIONS:
        return 0
    return DILUTION_WEIGHT_MAP.get(dilution, 0)


def find_best_dilution_any_lysis(dilution_counts: Mapping[int, Counter[str]]) -> Optional[int]:
    positive_dilutions = [
        d for d, counts in dilution_counts.items() if counts["1"] > 0 and d not in EXCLUDED_LOG_DILUTIONS
    ]
    if not positive_dilutions:
        return None
    return min(positive_dilutions)


def build_uncertainty_flags(
    *,
    total_obs: int,
    score_1_count: int,
    score_0_count: int,
    score_n_count: int,
    dilution_counts: Mapping[int, Counter[str]],
    expected_observations_per_pair: int,
    min_interpretable_obs_for_negative: int,
    high_uninterpretable_fraction_threshold: float,
) -> List[str]:
    flags: List[str] = []
    interpretable = score_1_count + score_0_count

    if score_n_count > 0:
        flags.append("has_uninterpretable")
    if total_obs > 0 and score_n_count / total_obs >= high_uninterpretable_fraction_threshold:
        flags.append("high_uninterpretable_fraction")
    if score_1_count > 0 and score_0_count > 0:
        flags.append("conflicting_interpretable_observations")
    if total_obs != expected_observations_per_pair:
        flags.append("incomplete_observation_grid")
    if interpretable < min_interpretable_obs_for_negative:
        flags.append("low_interpretable_support")

    within_dilution_conflict = False
    for counts in dilution_counts.values():
        if counts["1"] > 0 and counts["0"] > 0:
            within_dilution_conflict = True
            break
    if within_dilution_conflict:
        flags.append("within_dilution_conflict")

    return flags


def compute_label_v1(
    *,
    score_1_count: int,
    score_0_count: int,
    score_n_count: int,
    total_obs: int,
    dilution_counts: Mapping[int, Counter[str]],
    policy: LabelPolicyV1,
) -> Dict[str, object]:
    interpretable = score_1_count + score_0_count
    positive_fraction_interpretable = score_1_count / interpretable if interpretable > 0 else 0.0

    if score_1_count > 0:
        any_lysis = "1"
        label_reason = "at_least_one_lysis_observed"
    elif score_0_count >= policy.min_interpretable_obs_for_hard_negative:
        any_lysis = "0"
        label_reason = "sufficient_interpretable_no_lysis_observed"
    else:
        any_lysis = ""
        label_reason = "insufficient_interpretable_support_for_hard_negative"

    if any_lysis == "1":
        if score_1_count >= 6 or positive_fraction_interpretable >= 0.66:
            lysis_strength = "high"
        elif score_1_count >= 3 or positive_fraction_interpretable >= 0.33:
            lysis_strength = "medium"
        else:
            lysis_strength = "low"
    elif any_lysis == "0":
        lysis_strength = "none"
    else:
        lysis_strength = "unresolved"

    best_dilution = find_best_dilution_any_lysis(dilution_counts)
    if any_lysis == "1" and best_dilution is not None:
        potency_label = potency_label_from_dilution(best_dilution)
        potency_rank = potency_rank_from_dilution(best_dilution)
    elif any_lysis == "0":
        potency_label = "none"
        potency_rank = 0
        best_dilution = None
    else:
        potency_label = "unresolved"
        potency_rank = ""
        best_dilution = None

    uncertainty_flags = build_uncertainty_flags(
        total_obs=total_obs,
        score_1_count=score_1_count,
        score_0_count=score_0_count,
        score_n_count=score_n_count,
        dilution_counts=dilution_counts,
        expected_observations_per_pair=policy.expected_observations_per_pair,
        min_interpretable_obs_for_negative=policy.min_interpretable_obs_for_hard_negative,
        high_uninterpretable_fraction_threshold=policy.high_uninterpretable_fraction_threshold,
    )
    if any_lysis == "":
        uncertainty_flags.append("unresolved_label")

    return {
        "any_lysis": any_lysis,
        "label_reason": label_reason,
        "lysis_strength": lysis_strength,
        "lysis_strength_score": round(positive_fraction_interpretable, 6),
        "dilution_potency": potency_label,
        "dilution_potency_rank": potency_rank,
        "best_lysis_dilution": "" if best_dilution is None else best_dilution,
        "potency_support": "any_positive_at_dilution" if any_lysis == "1" else "",
        "uncertainty_flags": uncertainty_flags,
    }


def compute_label_v2(
    *,
    score_1_count: int,
    score_0_count: int,
    score_n_count: int,
    total_obs: int,
    dilution_counts: Mapping[int, Counter[str]],
    policy: LabelPolicyV2,
) -> Dict[str, object]:
    interpretable = score_1_count + score_0_count
    positive_fraction_interpretable = score_1_count / interpretable if interpretable > 0 else 0.0

    if (
        score_1_count >= policy.min_positive_obs
        and positive_fraction_interpretable >= policy.min_positive_fraction_interpretable
        and score_n_count <= policy.max_uninterpretable_obs
    ):
        any_lysis = "1"
        label_reason = "strict_positive_rule"
    elif score_0_count >= policy.min_negative_obs and score_n_count <= policy.max_uninterpretable_obs:
        any_lysis = "0"
        label_reason = "strict_negative_rule"
    else:
        any_lysis = ""
        label_reason = "strict_rules_not_met"

    weighted_pos = 0
    weighted_interpretable = 0
    strong_support_dilutions: List[int] = []
    any_positive_dilutions: List[int] = []

    for dilution, counts in dilution_counts.items():
        interpretable_at_d = counts["1"] + counts["0"]
        if counts["1"] > 0:
            any_positive_dilutions.append(dilution)
        if counts["1"] >= 2 and interpretable_at_d >= 2:
            strong_support_dilutions.append(dilution)

        weight = DILUTION_WEIGHT_MAP.get(dilution, 1)
        weighted_pos += counts["1"] * weight
        weighted_interpretable += interpretable_at_d * weight

    weighted_strength = weighted_pos / weighted_interpretable if weighted_interpretable > 0 else 0.0

    if any_lysis == "1":
        if weighted_strength >= 0.65:
            lysis_strength = "high"
        elif weighted_strength >= 0.40:
            lysis_strength = "medium"
        else:
            lysis_strength = "low"
    elif any_lysis == "0":
        lysis_strength = "none"
    else:
        lysis_strength = "unresolved"

    best_dilution: Optional[int]
    potency_support: str
    if any_lysis == "1":
        if strong_support_dilutions:
            best_dilution = min(strong_support_dilutions)
            potency_support = "strong_support"
        elif any_positive_dilutions:
            best_dilution = min(any_positive_dilutions)
            potency_support = "weak_support"
        else:
            best_dilution = None
            potency_support = "unresolved"
    elif any_lysis == "0":
        best_dilution = None
        potency_support = "none"
    else:
        best_dilution = None
        potency_support = "unresolved"

    if any_lysis == "1" and best_dilution is not None:
        potency_label = potency_label_from_dilution(best_dilution)
        potency_rank = potency_rank_from_dilution(best_dilution)
    elif any_lysis == "0":
        potency_label = "none"
        potency_rank = 0
    else:
        potency_label = "unresolved"
        potency_rank = ""

    uncertainty_flags = build_uncertainty_flags(
        total_obs=total_obs,
        score_1_count=score_1_count,
        score_0_count=score_0_count,
        score_n_count=score_n_count,
        dilution_counts=dilution_counts,
        expected_observations_per_pair=policy.expected_observations_per_pair,
        min_interpretable_obs_for_negative=policy.min_negative_obs,
        high_uninterpretable_fraction_threshold=policy.high_uninterpretable_fraction_threshold,
    )
    if score_n_count > policy.max_uninterpretable_obs:
        uncertainty_flags.append("strict_high_uninterpretable")
    if any_lysis == "":
        uncertainty_flags.append("strict_policy_ambiguous")

    return {
        "any_lysis": any_lysis,
        "label_reason": label_reason,
        "lysis_strength": lysis_strength,
        "lysis_strength_score": round(weighted_strength, 6),
        "dilution_potency": potency_label,
        "dilution_potency_rank": potency_rank,
        "best_lysis_dilution": "" if best_dilution is None else best_dilution,
        "potency_support": potency_support,
        "uncertainty_flags": uncertainty_flags,
    }


def discover_images(search_roots: Sequence[Path]) -> Dict[str, str]:
    image_map: Dict[str, str] = {}
    for root in search_roots:
        if not root.exists() or not root.is_dir():
            continue
        for path in root.rglob("*"):
            if not path.is_file():
                continue
            if path.suffix.lower() not in IMAGE_EXTENSIONS:
                continue
            name = path.name
            if name not in image_map:
                image_map[name] = str(path)
    return image_map


def build_checks() -> List[Dict[str, object]]:
    return []


def add_check(
    checks: List[Dict[str, object]],
    *,
    check_name: str,
    passed: bool,
    severity_on_fail: str,
    pass_detail: str,
    fail_detail: str,
) -> None:
    checks.append(
        {
            "check_name": check_name,
            "status": "pass" if passed else "fail",
            "severity_on_fail": severity_on_fail,
            "detail": pass_detail if passed else fail_detail,
        }
    )


def to_pipe_sorted(values: Iterable[str]) -> str:
    return "|".join(sorted({v for v in values if v != ""}))


def main(argv: Optional[List[str]] = None) -> None:
    args = parse_args(argv)

    policy_v1 = LabelPolicyV1(
        expected_observations_per_pair=args.expected_observations_per_pair,
        min_interpretable_obs_for_hard_negative=args.v1_min_interpretable_neg,
        high_uninterpretable_fraction_threshold=args.v1_high_uninterpretable_frac,
    )
    policy_v2 = LabelPolicyV2(
        expected_observations_per_pair=args.expected_observations_per_pair,
        min_positive_obs=args.v2_min_positive_obs,
        min_positive_fraction_interpretable=args.v2_min_positive_fraction,
        min_negative_obs=args.v2_min_negative_obs,
        max_uninterpretable_obs=args.v2_max_uninterpretable_obs,
        high_uninterpretable_fraction_threshold=args.v2_high_uninterpretable_frac,
    )

    output_dir = args.output_dir
    id_map_dir = output_dir / "id_map"
    integrity_dir = output_dir / "integrity"
    cohort_dir = output_dir / "cohort"
    labels_dir = output_dir / "labels"
    qc_dir = output_dir / "qc"

    for path in (id_map_dir, integrity_dir, cohort_dir, labels_dir, qc_dir):
        ensure_directory(path)

    raw_rows_all, raw_header = read_delimited_rows(args.raw_interactions_path, delimiter=";")
    raw_rows = filter_excluded_dilutions(raw_rows_all)
    matrix_rows, matrix_header = read_delimited_rows(args.interaction_matrix_path, delimiter=";")
    host_rows, host_header = read_delimited_rows(args.host_metadata_path, delimiter=";")
    phage_rows, phage_header = read_delimited_rows(args.phage_metadata_path, delimiter=";")
    cv_rows, cv_header = read_delimited_rows(args.cv_groups_path, delimiter=";")

    missing_raw = sorted(set(RAW_REQUIRED_COLUMNS) - set(raw_header))
    missing_host = sorted(set(HOST_REQUIRED_COLUMNS) - set(host_header))
    missing_phage = sorted(set(PHAGE_REQUIRED_COLUMNS) - set(phage_header))
    missing_cv = sorted(set(CV_REQUIRED_COLUMNS) - set(cv_header))
    if missing_raw:
        raise ValueError(f"Missing raw columns: {', '.join(missing_raw)}")
    if missing_host:
        raise ValueError(f"Missing host metadata columns: {', '.join(missing_host)}")
    if missing_phage:
        raise ValueError(f"Missing phage metadata columns: {', '.join(missing_phage)}")
    if missing_cv:
        raise ValueError(f"Missing CV group columns: {', '.join(missing_cv)}")
    if "bacteria" not in matrix_header:
        raise ValueError("interaction_matrix.csv is missing bacteria column")

    matrix_phage_columns = [col for col in matrix_header if col != "bacteria"]

    override_map = read_alias_overrides(args.alias_overrides_path)

    raw_bacteria = {row["bacteria"] for row in raw_rows}
    raw_phages = {row["phage"] for row in raw_rows}
    matrix_bacteria = {row["bacteria"] for row in matrix_rows}
    matrix_phages = set(matrix_phage_columns)
    host_bacteria = {row["bacteria"] for row in host_rows}
    phage_metadata_phages = {row["phage"] for row in phage_rows}
    cv_bacteria = {row["bacteria"] for row in cv_rows}

    bacteria_by_source: Dict[str, Set[str]] = {
        "raw_interactions": raw_bacteria,
        "interaction_matrix": matrix_bacteria,
        "host_metadata": host_bacteria,
        "cv_groups": cv_bacteria,
    }
    phage_by_source: Dict[str, Set[str]] = {
        "raw_interactions": raw_phages,
        "interaction_matrix": matrix_phages,
        "phage_metadata": phage_metadata_phages,
    }

    bacteria_resolution_rows: List[Dict[str, object]] = []
    phage_resolution_rows: List[Dict[str, object]] = []

    canonical_bacteria_by_source_name: Dict[str, Dict[str, str]] = defaultdict(dict)
    canonical_phage_by_source_name: Dict[str, Dict[str, str]] = defaultdict(dict)

    for source, names in bacteria_by_source.items():
        for name in sorted(names):
            canonical, applied, reason = resolve_alias(name, "bacteria", override_map)
            canonical_bacteria_by_source_name[source][name] = canonical
            bacteria_resolution_rows.append(
                {
                    "entity_type": "bacteria",
                    "source": source,
                    "original_name": name,
                    "canonical_name": canonical,
                    "alias_applied": int(applied),
                    "alias_reason": reason,
                }
            )

    for source, names in phage_by_source.items():
        for name in sorted(names):
            canonical, applied, reason = resolve_alias(name, "phage", override_map)
            canonical_phage_by_source_name[source][name] = canonical
            phage_resolution_rows.append(
                {
                    "entity_type": "phage",
                    "source": source,
                    "original_name": name,
                    "canonical_name": canonical,
                    "alias_applied": int(applied),
                    "alias_reason": reason,
                }
            )

    canonical_bacteria_names = sorted(
        {canonical for per_source in canonical_bacteria_by_source_name.values() for canonical in per_source.values()}
    )
    canonical_phage_names = sorted(
        {canonical for per_source in canonical_phage_by_source_name.values() for canonical in per_source.values()}
    )

    bacteria_id_map = {name: f"BAC{idx:04d}" for idx, name in enumerate(canonical_bacteria_names, start=1)}
    phage_id_map = {name: f"PHG{idx:04d}" for idx, name in enumerate(canonical_phage_names, start=1)}

    raw_obs_count_by_bacteria: Counter[str] = Counter()
    raw_pair_count_by_bacteria: Dict[str, Set[str]] = defaultdict(set)
    raw_obs_count_by_phage: Counter[str] = Counter()
    raw_pair_count_by_phage: Dict[str, Set[str]] = defaultdict(set)

    for row in raw_rows:
        bacteria_raw = row["bacteria"]
        phage_raw = row["phage"]
        bacteria_can = canonical_bacteria_by_source_name["raw_interactions"][bacteria_raw]
        phage_can = canonical_phage_by_source_name["raw_interactions"][phage_raw]
        raw_obs_count_by_bacteria[bacteria_can] += 1
        raw_pair_count_by_bacteria[bacteria_can].add(phage_can)
        raw_obs_count_by_phage[phage_can] += 1
        raw_pair_count_by_phage[phage_can].add(bacteria_can)

    bacteria_sources_by_canonical: Dict[str, Set[str]] = defaultdict(set)
    bacteria_raw_names_by_canonical: Dict[str, Set[str]] = defaultdict(set)
    for source, mapping in canonical_bacteria_by_source_name.items():
        for raw_name, canonical in mapping.items():
            bacteria_sources_by_canonical[canonical].add(source)
            bacteria_raw_names_by_canonical[canonical].add(raw_name)

    phage_sources_by_canonical: Dict[str, Set[str]] = defaultdict(set)
    phage_raw_names_by_canonical: Dict[str, Set[str]] = defaultdict(set)
    for source, mapping in canonical_phage_by_source_name.items():
        for raw_name, canonical in mapping.items():
            phage_sources_by_canonical[canonical].add(source)
            phage_raw_names_by_canonical[canonical].add(raw_name)

    bacteria_id_rows: List[Dict[str, object]] = []
    for canonical in canonical_bacteria_names:
        sources = bacteria_sources_by_canonical[canonical]
        bacteria_id_rows.append(
            {
                "canonical_bacteria": canonical,
                "canonical_bacteria_id": bacteria_id_map[canonical],
                "source_presence": to_pipe_sorted(sources),
                "raw_names": to_pipe_sorted(bacteria_raw_names_by_canonical[canonical]),
                "in_raw_interactions": int("raw_interactions" in sources),
                "in_interaction_matrix": int("interaction_matrix" in sources),
                "in_host_metadata": int("host_metadata" in sources),
                "in_cv_groups": int("cv_groups" in sources),
                "raw_observation_count": raw_obs_count_by_bacteria.get(canonical, 0),
                "raw_unique_phage_count": len(raw_pair_count_by_bacteria.get(canonical, set())),
            }
        )

    phage_id_rows: List[Dict[str, object]] = []
    for canonical in canonical_phage_names:
        sources = phage_sources_by_canonical[canonical]
        phage_id_rows.append(
            {
                "canonical_phage": canonical,
                "canonical_phage_id": phage_id_map[canonical],
                "source_presence": to_pipe_sorted(sources),
                "raw_names": to_pipe_sorted(phage_raw_names_by_canonical[canonical]),
                "in_raw_interactions": int("raw_interactions" in sources),
                "in_interaction_matrix": int("interaction_matrix" in sources),
                "in_phage_metadata": int("phage_metadata" in sources),
                "raw_observation_count": raw_obs_count_by_phage.get(canonical, 0),
                "raw_unique_bacteria_count": len(raw_pair_count_by_phage.get(canonical, set())),
            }
        )

    bacteria_alias_candidates = build_alias_candidates(
        bacteria_by_source,
        "bacteria",
        override_map,
    )
    phage_alias_candidates = build_alias_candidates(
        phage_by_source,
        "phage",
        override_map,
    )

    write_csv(
        id_map_dir / "bacteria_id_map.csv",
        fieldnames=list(bacteria_id_rows[0].keys()) if bacteria_id_rows else [],
        rows=bacteria_id_rows,
    )
    write_csv(
        id_map_dir / "phage_id_map.csv",
        fieldnames=list(phage_id_rows[0].keys()) if phage_id_rows else [],
        rows=phage_id_rows,
    )
    write_csv(
        id_map_dir / "bacteria_alias_resolution.csv",
        fieldnames=list(bacteria_resolution_rows[0].keys()) if bacteria_resolution_rows else [],
        rows=bacteria_resolution_rows,
    )
    write_csv(
        id_map_dir / "phage_alias_resolution.csv",
        fieldnames=list(phage_resolution_rows[0].keys()) if phage_resolution_rows else [],
        rows=phage_resolution_rows,
    )

    if bacteria_alias_candidates:
        write_csv(
            id_map_dir / "bacteria_alias_candidates.csv",
            fieldnames=list(bacteria_alias_candidates[0].keys()),
            rows=bacteria_alias_candidates,
        )
    else:
        write_csv(
            id_map_dir / "bacteria_alias_candidates.csv",
            fieldnames=[
                "entity_type",
                "left_name",
                "right_name",
                "left_sources",
                "right_sources",
                "match_reason",
                "similarity_ratio",
                "confidence_score",
                "review_status",
            ],
            rows=[],
        )

    if phage_alias_candidates:
        write_csv(
            id_map_dir / "phage_alias_candidates.csv",
            fieldnames=list(phage_alias_candidates[0].keys()),
            rows=phage_alias_candidates,
        )
    else:
        write_csv(
            id_map_dir / "phage_alias_candidates.csv",
            fieldnames=[
                "entity_type",
                "left_name",
                "right_name",
                "left_sources",
                "right_sources",
                "match_reason",
                "similarity_ratio",
                "confidence_score",
                "review_status",
            ],
            rows=[],
        )

    checks: List[Dict[str, object]] = build_checks()

    raw_score_counter = Counter(row["score"] for row in raw_rows)
    invalid_scores = sorted(set(raw_score_counter.keys()) - ALLOWED_RAW_SCORES)
    add_check(
        checks,
        check_name="raw_allowed_scores",
        passed=len(invalid_scores) == 0,
        severity_on_fail="error",
        pass_detail="raw score values are limited to {0,1,n}",
        fail_detail=f"invalid score values found: {', '.join(invalid_scores)}",
    )

    raw_duplicate_counter: Counter[Tuple[str, ...]] = Counter()
    for row in raw_rows:
        key = (
            row["bacteria"],
            row["phage"],
            row["image"],
            row["replicate"],
            row["plate"],
            row["log_dilution"],
            row["X"],
            row["Y"],
            row["score"],
        )
        raw_duplicate_counter[key] += 1
    duplicate_obs_count = sum(count - 1 for count in raw_duplicate_counter.values() if count > 1)
    add_check(
        checks,
        check_name="raw_duplicate_observation_rows",
        passed=duplicate_obs_count == 0,
        severity_on_fail="warning",
        pass_detail="no exact duplicate raw observation rows detected",
        fail_detail=f"duplicate raw observation row count={duplicate_obs_count}",
    )

    raw_pair_counter: Counter[Tuple[str, str]] = Counter()
    for row in raw_rows:
        bacteria_can = canonical_bacteria_by_source_name["raw_interactions"][row["bacteria"]]
        phage_can = canonical_phage_by_source_name["raw_interactions"][row["phage"]]
        raw_pair_counter[(bacteria_can, phage_can)] += 1

    raw_pairs_with_non9 = sorted(
        [pair for pair, n in raw_pair_counter.items() if n != args.expected_observations_per_pair]
    )
    add_check(
        checks,
        check_name="raw_pair_observation_grid_count",
        passed=len(raw_pairs_with_non9) == 0,
        severity_on_fail="warning",
        pass_detail=f"all raw pairs have {args.expected_observations_per_pair} observations",
        fail_detail=(
            f"{len(raw_pairs_with_non9)} pairs deviate from expected observation count "
            f"{args.expected_observations_per_pair}"
        ),
    )

    expected_raw_rows = len(raw_bacteria) * len(raw_phages) * args.expected_observations_per_pair
    add_check(
        checks,
        check_name="raw_total_rows_vs_full_grid",
        passed=len(raw_rows) == expected_raw_rows,
        severity_on_fail="warning",
        pass_detail=(
            f"raw rows={len(raw_rows)} matches full grid {len(raw_bacteria)}*{len(raw_phages)}*"
            f"{args.expected_observations_per_pair}"
        ),
        fail_detail=(f"raw rows={len(raw_rows)} does not match expected full grid rows={expected_raw_rows}"),
    )

    add_check(
        checks,
        check_name="raw_phages_present_in_matrix_columns",
        passed=len(raw_phages - matrix_phages) == 0,
        severity_on_fail="error",
        pass_detail="all raw phages are present in interaction_matrix columns",
        fail_detail=f"raw phages missing from matrix columns={len(raw_phages - matrix_phages)}",
    )

    add_check(
        checks,
        check_name="raw_phages_present_in_phage_metadata",
        passed=len(raw_phages - phage_metadata_phages) == 0,
        severity_on_fail="error",
        pass_detail="all raw phages are present in phage metadata",
        fail_detail=f"raw phages missing from phage metadata={len(raw_phages - phage_metadata_phages)}",
    )

    add_check(
        checks,
        check_name="raw_bacteria_present_in_host_metadata",
        passed=len(raw_bacteria - host_bacteria) == 0,
        severity_on_fail="error",
        pass_detail="all raw bacteria are present in host metadata",
        fail_detail=f"raw bacteria missing from host metadata={len(raw_bacteria - host_bacteria)}",
    )

    add_check(
        checks,
        check_name="raw_bacteria_present_in_cv_groups",
        passed=len(raw_bacteria - cv_bacteria) == 0,
        severity_on_fail="error",
        pass_detail="all raw bacteria are present in CV groups",
        fail_detail=f"raw bacteria missing from cv groups={len(raw_bacteria - cv_bacteria)}",
    )

    matrix_duplicate_bacteria = len(matrix_rows) - len(matrix_bacteria)
    host_duplicate_bacteria = len(host_rows) - len(host_bacteria)
    cv_duplicate_bacteria = len(cv_rows) - len(cv_bacteria)
    phage_duplicate_ids = len(phage_rows) - len(phage_metadata_phages)

    add_check(
        checks,
        check_name="matrix_unique_bacteria_rows",
        passed=matrix_duplicate_bacteria == 0,
        severity_on_fail="warning",
        pass_detail="no duplicate bacteria rows in interaction_matrix",
        fail_detail=f"duplicate bacteria rows in interaction_matrix={matrix_duplicate_bacteria}",
    )
    add_check(
        checks,
        check_name="host_unique_bacteria_rows",
        passed=host_duplicate_bacteria == 0,
        severity_on_fail="warning",
        pass_detail="no duplicate bacteria rows in host metadata",
        fail_detail=f"duplicate bacteria rows in host metadata={host_duplicate_bacteria}",
    )
    add_check(
        checks,
        check_name="cv_unique_bacteria_rows",
        passed=cv_duplicate_bacteria == 0,
        severity_on_fail="warning",
        pass_detail="no duplicate bacteria rows in cv groups",
        fail_detail=f"duplicate bacteria rows in cv groups={cv_duplicate_bacteria}",
    )
    add_check(
        checks,
        check_name="phage_unique_rows",
        passed=phage_duplicate_ids == 0,
        severity_on_fail="warning",
        pass_detail="no duplicate phage rows in phage metadata",
        fail_detail=f"duplicate phage rows in phage metadata={phage_duplicate_ids}",
    )

    matrix_parse_errors = 0
    for row in matrix_rows:
        for phage in matrix_phage_columns:
            value = row.get(phage, "")
            if value == "":
                continue
            if parse_float_or_none(value) is None:
                matrix_parse_errors += 1
    add_check(
        checks,
        check_name="matrix_numeric_cells_parseable",
        passed=matrix_parse_errors == 0,
        severity_on_fail="warning",
        pass_detail="all non-empty matrix cells parse as float",
        fail_detail=f"matrix parse failures={matrix_parse_errors}",
    )

    issues = [check for check in checks if check["status"] == "fail"]
    error_issue_count = sum(1 for issue in issues if issue["severity_on_fail"] == "error")
    warning_issue_count = sum(1 for issue in issues if issue["severity_on_fail"] == "warning")

    integrity_summary = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "status": "pass" if error_issue_count == 0 else "fail",
        "input_counts": {
            "raw_rows": len(raw_rows),
            "raw_bacteria": len(raw_bacteria),
            "raw_phages": len(raw_phages),
            "matrix_rows": len(matrix_rows),
            "matrix_bacteria": len(matrix_bacteria),
            "matrix_phages": len(matrix_phages),
            "host_rows": len(host_rows),
            "host_bacteria": len(host_bacteria),
            "phage_metadata_rows": len(phage_rows),
            "cv_rows": len(cv_rows),
            "cv_bacteria": len(cv_bacteria),
        },
        "failed_check_count": len(issues),
        "failed_error_check_count": error_issue_count,
        "failed_warning_check_count": warning_issue_count,
        "checks": checks,
    }

    write_csv(
        integrity_dir / "integrity_checks.csv",
        fieldnames=list(checks[0].keys()) if checks else ["check_name", "status", "severity_on_fail", "detail"],
        rows=checks,
    )
    write_json(integrity_dir / "integrity_report.json", integrity_summary)

    # Cohort contracts and denominator definitions.
    cohorts = {
        "raw369": {
            "description": "Bacteria present in raw_interactions.csv",
            "bacteria": sorted({canonical_bacteria_by_source_name["raw_interactions"][name] for name in raw_bacteria}),
        },
        "matrix402": {
            "description": "Bacteria present in interaction_matrix.csv",
            "bacteria": sorted(
                {canonical_bacteria_by_source_name["interaction_matrix"][name] for name in matrix_bacteria}
            ),
        },
        "features404": {
            "description": "Bacteria present in CV group table (feature-ready denominator)",
            "bacteria": sorted({canonical_bacteria_by_source_name["cv_groups"][name] for name in cv_bacteria}),
        },
    }

    host_bacteria_canonical = {canonical_bacteria_by_source_name["host_metadata"][name] for name in host_bacteria}
    matrix_bacteria_canonical = {
        canonical_bacteria_by_source_name["interaction_matrix"][name] for name in matrix_bacteria
    }
    raw_bacteria_canonical = {canonical_bacteria_by_source_name["raw_interactions"][name] for name in raw_bacteria}
    cv_bacteria_canonical = {canonical_bacteria_by_source_name["cv_groups"][name] for name in cv_bacteria}

    cohort_rows: List[Dict[str, object]] = []
    for cohort_name, cohort in cohorts.items():
        bacteria_set = set(cohort["bacteria"])
        cohort_rows.append(
            {
                "cohort_name": cohort_name,
                "description": cohort["description"],
                "bacteria_count": len(bacteria_set),
                "phage_panel_count": len(canonical_phage_names),
                "pair_grid_size": len(bacteria_set) * len(canonical_phage_names),
                "in_raw_bacteria_count": len(bacteria_set & raw_bacteria_canonical),
                "in_matrix_bacteria_count": len(bacteria_set & matrix_bacteria_canonical),
                "in_host_metadata_bacteria_count": len(bacteria_set & host_bacteria_canonical),
                "in_cv_groups_bacteria_count": len(bacteria_set & cv_bacteria_canonical),
            }
        )

    cohort_contracts = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "cohort_contracts": cohorts,
        "denominator_rules": {
            "primary_training_reporting_cohort": "raw369",
            "matrix_auxiliary_reporting_cohort": "matrix402",
            "feature_universe_reporting_cohort": "features404",
            "phage_denominator": "all canonical phages from union of raw/matrix/phage_metadata",
        },
    }

    write_csv(
        cohort_dir / "cohort_contracts.csv",
        fieldnames=list(cohort_rows[0].keys()) if cohort_rows else [],
        rows=cohort_rows,
    )
    write_json(cohort_dir / "cohort_contracts.json", cohort_contracts)

    # Preserve observation-level replicate/dilution structure with canonical IDs.
    observation_rows: List[Dict[str, object]] = []
    pair_dilution_counts: Dict[Tuple[str, str], Dict[int, Counter[str]]] = defaultdict(lambda: defaultdict(Counter))
    pair_grid_scores: Dict[Tuple[str, str], Dict[Tuple[str, int], Counter[str]]] = defaultdict(
        lambda: defaultdict(Counter)
    )

    for idx, row in enumerate(raw_rows, start=1):
        bacteria_raw = row["bacteria"]
        phage_raw = row["phage"]
        bacteria_can = canonical_bacteria_by_source_name["raw_interactions"][bacteria_raw]
        phage_can = canonical_phage_by_source_name["raw_interactions"][phage_raw]
        pair_id = f"{bacteria_can}__{phage_can}"
        dilution = parse_int_or_none(row["log_dilution"])
        if dilution is None:
            raise ValueError(f"Invalid log_dilution value: {row['log_dilution']!r}")

        score = row["score"]
        pair_dilution_counts[(bacteria_can, phage_can)][dilution][score] += 1
        pair_grid_scores[(bacteria_can, phage_can)][(row["replicate"], dilution)][score] += 1

        observation_rows.append(
            {
                "observation_id": idx,
                "pair_id": pair_id,
                "bacteria_raw": bacteria_raw,
                "bacteria": bacteria_can,
                "bacteria_id": bacteria_id_map[bacteria_can],
                "phage_raw": phage_raw,
                "phage": phage_can,
                "phage_id": phage_id_map[phage_can],
                "bacteria_index": row["bacteria_index"],
                "image": row["image"],
                "replicate": row["replicate"],
                "plate": row["plate"],
                "log_dilution": dilution,
                "x_coord": row["X"],
                "y_coord": row["Y"],
                "score": score,
            }
        )

    write_csv(
        labels_dir / "track_a_observations_with_ids.csv",
        fieldnames=list(observation_rows[0].keys()) if observation_rows else [],
        rows=observation_rows,
    )

    dilution_rows: List[Dict[str, object]] = []
    grid_rows: List[Dict[str, object]] = []

    for pair_key in sorted(pair_dilution_counts.keys()):
        bacteria_can, phage_can = pair_key
        pair_id = f"{bacteria_can}__{phage_can}"
        dilution_map = pair_dilution_counts[pair_key]
        for dilution in sorted(dilution_map.keys()):
            counts = dilution_map[dilution]
            interpretable = counts["0"] + counts["1"]
            pos_frac = counts["1"] / interpretable if interpretable > 0 else 0.0
            dilution_rows.append(
                {
                    "pair_id": pair_id,
                    "bacteria": bacteria_can,
                    "phage": phage_can,
                    "log_dilution": dilution,
                    "n_obs": sum(counts.values()),
                    "score_1_count": counts["1"],
                    "score_0_count": counts["0"],
                    "score_n_count": counts["n"],
                    "interpretable_count": interpretable,
                    "positive_fraction_interpretable": round(pos_frac, 6),
                }
            )

        grid_map = pair_grid_scores[pair_key]
        for rep, dilution in sorted(grid_map.keys(), key=lambda key: (key[0], key[1])):
            counts = grid_map[(rep, dilution)]
            dominant_score = ""
            if counts:
                dominant_score = sorted(counts.items(), key=lambda item: (-item[1], item[0]))[0][0]
            grid_rows.append(
                {
                    "pair_id": pair_id,
                    "bacteria": bacteria_can,
                    "phage": phage_can,
                    "replicate": rep,
                    "log_dilution": dilution,
                    "n_obs": sum(counts.values()),
                    "score_1_count": counts["1"],
                    "score_0_count": counts["0"],
                    "score_n_count": counts["n"],
                    "dominant_score": dominant_score,
                }
            )

    write_csv(
        labels_dir / "track_a_pair_dilution_summary.csv",
        fieldnames=list(dilution_rows[0].keys()) if dilution_rows else [],
        rows=dilution_rows,
    )
    write_csv(
        labels_dir / "track_a_pair_observation_grid.csv",
        fieldnames=list(grid_rows[0].keys()) if grid_rows else [],
        rows=grid_rows,
    )

    pair_rows_v1: List[Dict[str, object]] = []
    pair_rows_v2: List[Dict[str, object]] = []

    for pair_key in sorted(pair_dilution_counts.keys()):
        bacteria_can, phage_can = pair_key
        pair_id = f"{bacteria_can}__{phage_can}"
        dilution_map = pair_dilution_counts[pair_key]

        score_1_count = sum(counts["1"] for counts in dilution_map.values())
        score_0_count = sum(counts["0"] for counts in dilution_map.values())
        score_n_count = sum(counts["n"] for counts in dilution_map.values())
        total_obs = score_1_count + score_0_count + score_n_count
        interpretable = score_1_count + score_0_count
        positive_fraction = score_1_count / interpretable if interpretable > 0 else 0.0

        v1 = compute_label_v1(
            score_1_count=score_1_count,
            score_0_count=score_0_count,
            score_n_count=score_n_count,
            total_obs=total_obs,
            dilution_counts=dilution_map,
            policy=policy_v1,
        )
        v2 = compute_label_v2(
            score_1_count=score_1_count,
            score_0_count=score_0_count,
            score_n_count=score_n_count,
            total_obs=total_obs,
            dilution_counts=dilution_map,
            policy=policy_v2,
        )

        base = {
            "pair_id": pair_id,
            "bacteria": bacteria_can,
            "bacteria_id": bacteria_id_map[bacteria_can],
            "phage": phage_can,
            "phage_id": phage_id_map[phage_can],
            "total_obs": total_obs,
            "score_1_count": score_1_count,
            "score_0_count": score_0_count,
            "score_n_count": score_n_count,
            "interpretable_count": interpretable,
            "positive_fraction_interpretable": round(positive_fraction, 6),
            "n_distinct_replicates": len({rep for rep, _ in pair_grid_scores[pair_key].keys()}),
            "n_distinct_dilutions": len(dilution_map),
        }

        row_v1 = dict(base)
        row_v1.update(
            {
                "any_lysis": v1["any_lysis"],
                "label_reason": v1["label_reason"],
                "lysis_strength": v1["lysis_strength"],
                "lysis_strength_score": v1["lysis_strength_score"],
                "dilution_potency": v1["dilution_potency"],
                "dilution_potency_rank": v1["dilution_potency_rank"],
                "best_lysis_dilution": v1["best_lysis_dilution"],
                "potency_support": v1["potency_support"],
                "uncertainty_flags": "|".join(v1["uncertainty_flags"]),
                "uncertainty_flag_count": len(v1["uncertainty_flags"]),
                "include_in_training": int(v1["any_lysis"] in {"0", "1"}),
            }
        )

        row_v2 = dict(base)
        row_v2.update(
            {
                "any_lysis": v2["any_lysis"],
                "label_reason": v2["label_reason"],
                "lysis_strength": v2["lysis_strength"],
                "lysis_strength_score": v2["lysis_strength_score"],
                "dilution_potency": v2["dilution_potency"],
                "dilution_potency_rank": v2["dilution_potency_rank"],
                "best_lysis_dilution": v2["best_lysis_dilution"],
                "potency_support": v2["potency_support"],
                "uncertainty_flags": "|".join(v2["uncertainty_flags"]),
                "uncertainty_flag_count": len(v2["uncertainty_flags"]),
                "include_in_training": int(v2["any_lysis"] in {"0", "1"}),
            }
        )

        pair_rows_v1.append(row_v1)
        pair_rows_v2.append(row_v2)

    write_csv(
        labels_dir / "label_set_v1_pairs.csv",
        fieldnames=list(pair_rows_v1[0].keys()) if pair_rows_v1 else [],
        rows=pair_rows_v1,
    )
    write_csv(
        labels_dir / "label_set_v2_pairs.csv",
        fieldnames=list(pair_rows_v2[0].keys()) if pair_rows_v2 else [],
        rows=pair_rows_v2,
    )

    label_set_v1_policy = {
        "policy_name": "track_a_label_set_v1",
        "policy_version": "v1",
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "any_lysis_rule": {
            "positive": "any score=1 across 9 observations",
            "negative": ("no score=1 and score_0_count >= min_interpretable_obs_for_hard_negative"),
            "unresolved": "otherwise",
        },
        "lysis_strength_rule": {
            "high": "positive and (score_1_count >= 6 or positive_fraction_interpretable >= 0.66)",
            "medium": "positive and (score_1_count >= 3 or positive_fraction_interpretable >= 0.33)",
            "low": "positive otherwise",
            "none": "hard negative",
            "unresolved": "unresolved any_lysis",
        },
        "dilution_potency_rule": ("best dilution where any lysis observed (most diluted/highest potency retained)"),
        "uncertainty_flags": [
            "has_uninterpretable",
            "high_uninterpretable_fraction",
            "conflicting_interpretable_observations",
            "incomplete_observation_grid",
            "low_interpretable_support",
            "within_dilution_conflict",
            "unresolved_label",
        ],
        "thresholds": asdict(policy_v1),
    }

    label_set_v2_policy = {
        "policy_name": "track_a_label_set_v2",
        "policy_version": "v1",
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "any_lysis_rule": {
            "positive": (
                "score_1_count >= min_positive_obs and positive_fraction_interpretable >= "
                "min_positive_fraction_interpretable and score_n_count <= max_uninterpretable_obs"
            ),
            "negative": ("score_0_count >= min_negative_obs and score_n_count <= max_uninterpretable_obs"),
            "unresolved": "otherwise",
        },
        "lysis_strength_rule": ("weighted positive fraction across dilutions using potency weights"),
        "dilution_potency_rule": (
            "best dilution with strong support (>=2 positives), else weak support if positive label"
        ),
        "uncertainty_flags": [
            "has_uninterpretable",
            "high_uninterpretable_fraction",
            "conflicting_interpretable_observations",
            "incomplete_observation_grid",
            "low_interpretable_support",
            "within_dilution_conflict",
            "strict_high_uninterpretable",
            "strict_policy_ambiguous",
        ],
        "thresholds": asdict(policy_v2),
    }

    write_json(labels_dir / "label_set_v1_policy.json", label_set_v1_policy)
    write_json(labels_dir / "label_set_v2_policy.json", label_set_v2_policy)

    v1_label_counts = Counter(row["any_lysis"] if row["any_lysis"] != "" else "unresolved" for row in pair_rows_v1)
    v2_label_counts = Counter(row["any_lysis"] if row["any_lysis"] != "" else "unresolved" for row in pair_rows_v2)

    v1_summary = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "pair_count": len(pair_rows_v1),
        "label_counts": {
            "positive_1": v1_label_counts["1"],
            "negative_0": v1_label_counts["0"],
            "unresolved": v1_label_counts["unresolved"],
        },
        "include_in_training_count": sum(row["include_in_training"] for row in pair_rows_v1),
        "strict_uncertainty_counts": Counter(
            flag for row in pair_rows_v1 for flag in row["uncertainty_flags"].split("|") if flag
        ),
    }

    v2_summary = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "pair_count": len(pair_rows_v2),
        "label_counts": {
            "positive_1": v2_label_counts["1"],
            "negative_0": v2_label_counts["0"],
            "unresolved": v2_label_counts["unresolved"],
        },
        "include_in_training_count": sum(row["include_in_training"] for row in pair_rows_v2),
        "strict_uncertainty_counts": Counter(
            flag for row in pair_rows_v2 for flag in row["uncertainty_flags"].split("|") if flag
        ),
    }

    write_json(labels_dir / "label_set_v1_summary.json", v1_summary)
    write_json(labels_dir / "label_set_v2_summary.json", v2_summary)

    v2_by_pair = {row["pair_id"]: row for row in pair_rows_v2}
    comparison_rows: List[Dict[str, object]] = []
    for row_v1 in pair_rows_v1:
        row_v2 = v2_by_pair[row_v1["pair_id"]]
        comparison_rows.append(
            {
                "pair_id": row_v1["pair_id"],
                "bacteria": row_v1["bacteria"],
                "phage": row_v1["phage"],
                "v1_any_lysis": row_v1["any_lysis"],
                "v2_any_lysis": row_v2["any_lysis"],
                "v1_lysis_strength": row_v1["lysis_strength"],
                "v2_lysis_strength": row_v2["lysis_strength"],
                "v1_dilution_potency": row_v1["dilution_potency"],
                "v2_dilution_potency": row_v2["dilution_potency"],
                "changed_any_lysis": int(row_v1["any_lysis"] != row_v2["any_lysis"]),
                "changed_lysis_strength": int(row_v1["lysis_strength"] != row_v2["lysis_strength"]),
                "changed_dilution_potency": int(row_v1["dilution_potency"] != row_v2["dilution_potency"]),
                "v1_uncertainty_flags": row_v1["uncertainty_flags"],
                "v2_uncertainty_flags": row_v2["uncertainty_flags"],
            }
        )

    write_csv(
        labels_dir / "label_set_v1_v2_comparison.csv",
        fieldnames=list(comparison_rows[0].keys()) if comparison_rows else [],
        rows=comparison_rows,
    )

    # Plaque-image-assisted QC queue for ambiguous/conflicting pairs.
    v1_by_pair = {row["pair_id"]: row for row in pair_rows_v1}
    v2_by_pair = {row["pair_id"]: row for row in pair_rows_v2}

    qc_reason_by_pair: Dict[str, Set[str]] = defaultdict(set)
    for pair_id, row in v1_by_pair.items():
        flags = set(flag for flag in str(row["uncertainty_flags"]).split("|") if flag)
        if "conflicting_interpretable_observations" in flags:
            qc_reason_by_pair[pair_id].add("conflicting_interpretable_observations")
        if "high_uninterpretable_fraction" in flags:
            qc_reason_by_pair[pair_id].add("high_uninterpretable_fraction")
        if "unresolved_label" in flags:
            qc_reason_by_pair[pair_id].add("v1_unresolved")

    for pair_id, row in v2_by_pair.items():
        flags = set(flag for flag in str(row["uncertainty_flags"]).split("|") if flag)
        if "strict_policy_ambiguous" in flags:
            qc_reason_by_pair[pair_id].add("v2_ambiguous")
        if v1_by_pair[pair_id]["any_lysis"] != row["any_lysis"]:
            qc_reason_by_pair[pair_id].add("v1_v2_label_disagreement")

    raw_obs_by_pair: Dict[str, List[Dict[str, object]]] = defaultdict(list)
    for row in observation_rows:
        raw_obs_by_pair[str(row["pair_id"])].append(row)

    image_search_roots = [Path(part.strip()) for part in args.image_search_roots.split(",") if part.strip() != ""]
    image_name_to_path = discover_images(image_search_roots)

    qc_rows: List[Dict[str, object]] = []
    for pair_id in sorted(qc_reason_by_pair.keys()):
        reasons = sorted(qc_reason_by_pair[pair_id])
        for row in sorted(
            raw_obs_by_pair.get(pair_id, []),
            key=lambda item: (int(item["observation_id"]),),
        ):
            image_name = str(row["image"])
            image_path = image_name_to_path.get(image_name, "")
            qc_rows.append(
                {
                    "pair_id": pair_id,
                    "bacteria": row["bacteria"],
                    "phage": row["phage"],
                    "qc_reasons": "|".join(reasons),
                    "observation_id": row["observation_id"],
                    "image": image_name,
                    "image_found": int(image_path != ""),
                    "image_path": image_path,
                    "replicate": row["replicate"],
                    "plate": row["plate"],
                    "log_dilution": row["log_dilution"],
                    "score": row["score"],
                }
            )

    write_csv(
        qc_dir / "plaque_image_qc_queue.csv",
        fieldnames=list(qc_rows[0].keys())
        if qc_rows
        else [
            "pair_id",
            "bacteria",
            "phage",
            "qc_reasons",
            "observation_id",
            "image",
            "image_found",
            "image_path",
            "replicate",
            "plate",
            "log_dilution",
            "score",
        ],
        rows=qc_rows,
    )

    qc_summary = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "image_search_roots": [str(path) for path in image_search_roots],
        "indexed_image_file_count": len(image_name_to_path),
        "qc_pair_count": len(qc_reason_by_pair),
        "qc_observation_count": len(qc_rows),
        "qc_rows_with_found_image": sum(int(row["image_found"]) for row in qc_rows),
        "qc_rows_missing_image": sum(1 for row in qc_rows if int(row["image_found"]) == 0),
    }
    write_json(qc_dir / "plaque_image_qc_summary.json", qc_summary)

    # Write ID map summary and manifest.
    id_map_summary = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "canonical_bacteria_count": len(canonical_bacteria_names),
        "canonical_phage_count": len(canonical_phage_names),
        "manual_alias_counts": {
            "bacteria": len(override_map.get("bacteria", {})),
            "phage": len(override_map.get("phage", {})),
        },
        "auto_alias_candidate_counts": {
            "bacteria": len(bacteria_alias_candidates),
            "phage": len(phage_alias_candidates),
        },
    }
    write_json(id_map_dir / "id_map_summary.json", id_map_summary)

    manifest = {
        "generated_at_utc": datetime.now(tz=timezone.utc).isoformat(),
        "step_name": "build_track_a_foundation",
        "inputs": {
            "raw_interactions_path": str(args.raw_interactions_path),
            "interaction_matrix_path": str(args.interaction_matrix_path),
            "host_metadata_path": str(args.host_metadata_path),
            "phage_metadata_path": str(args.phage_metadata_path),
            "cv_groups_path": str(args.cv_groups_path),
            "alias_overrides_path": str(args.alias_overrides_path),
        },
        "input_hashes_sha256": {
            "raw_interactions": sha256_of_file(args.raw_interactions_path),
            "interaction_matrix": sha256_of_file(args.interaction_matrix_path),
            "host_metadata": sha256_of_file(args.host_metadata_path),
            "phage_metadata": sha256_of_file(args.phage_metadata_path),
            "cv_groups": sha256_of_file(args.cv_groups_path),
            "alias_overrides": sha256_of_file(args.alias_overrides_path) if args.alias_overrides_path.exists() else "",
        },
        "outputs": {
            "id_map_summary_json": str(id_map_dir / "id_map_summary.json"),
            "integrity_report_json": str(integrity_dir / "integrity_report.json"),
            "cohort_contracts_json": str(cohort_dir / "cohort_contracts.json"),
            "label_set_v1_summary_json": str(labels_dir / "label_set_v1_summary.json"),
            "label_set_v2_summary_json": str(labels_dir / "label_set_v2_summary.json"),
            "plaque_image_qc_summary_json": str(qc_dir / "plaque_image_qc_summary.json"),
        },
        "policies": {
            "label_set_v1": asdict(policy_v1),
            "label_set_v2": asdict(policy_v2),
        },
    }
    write_json(output_dir / "track_a_manifest.json", manifest)

    print("Track A foundation build completed.")
    print(f"- Output dir: {output_dir}")
    print(f"- Canonical bacteria: {len(canonical_bacteria_names)}")
    print(f"- Canonical phages: {len(canonical_phage_names)}")
    print(f"- Integrity status: {integrity_summary['status']}")
    print(f"- V1 pair labels: {len(pair_rows_v1)}")
    print(f"- V2 pair labels: {len(pair_rows_v2)}")
    print(f"- QC pair count: {len(qc_reason_by_pair)}")


if __name__ == "__main__":
    main()
