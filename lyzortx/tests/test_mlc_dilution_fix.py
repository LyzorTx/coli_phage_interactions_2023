"""SX05: Regression tests for MLC dilution correction (drop log_dilution=-4).

These guard the invariants the paper protocol requires:
  - DILUTION_WEIGHT_MAP contains exactly {0, -1, -2}, matching the three replicated concentrations.
  - A pair whose only lysis observation sits at log_dilution=-4 yields MLC=0 (not MLC=4).
  - find_best_dilution_any_lysis and potency helpers defensively skip EXCLUDED_LOG_DILUTIONS.
  - The regeneration helpers collapse old MLC=4 pairs (lysis at -4, none at -2) to MLC=0 and
    pairs with lysis at -2 AND -4 to MLC=3.
"""

from __future__ import annotations

from collections import Counter

from lyzortx.pipeline.autoresearch.regenerate_interaction_matrix import (
    build_dilution_counts_by_pair,
    compute_mlc_for_pair,
)
from lyzortx.pipeline.track_a.steps.build_track_a_foundation import (
    DILUTION_POTENCY_LABEL_MAP,
    DILUTION_WEIGHT_MAP,
    EXCLUDED_LOG_DILUTIONS,
    filter_excluded_dilutions,
    find_best_dilution_any_lysis,
    potency_label_from_dilution,
    potency_rank_from_dilution,
)


def test_dilution_weight_map_has_exactly_three_keys() -> None:
    assert set(DILUTION_WEIGHT_MAP) == {0, -1, -2}
    assert DILUTION_WEIGHT_MAP == {0: 1, -1: 2, -2: 3}


def test_potency_label_map_aligns_with_weight_map() -> None:
    assert set(DILUTION_POTENCY_LABEL_MAP) == set(DILUTION_WEIGHT_MAP)


def test_excluded_dilutions_contains_minus_four() -> None:
    assert -4 in EXCLUDED_LOG_DILUTIONS


def test_filter_excluded_dilutions_drops_minus_four_rows() -> None:
    rows = [
        {"bacteria": "B1", "phage": "P1", "log_dilution": "0", "score": "1"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-1", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-2", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-4", "score": "1"},
    ]
    kept = filter_excluded_dilutions(rows)
    assert len(kept) == 3
    assert all(int(r["log_dilution"]) != -4 for r in kept)


def test_find_best_dilution_ignores_excluded_even_if_present() -> None:
    # Defensive: even if a -4 count slipped through, the helper must not pick it.
    counts = {
        0: Counter({"0": 3}),
        -1: Counter({"0": 2}),
        -2: Counter({"0": 3}),
        -4: Counter({"1": 1}),  # would have been MIN positive under old code
    }
    assert find_best_dilution_any_lysis(counts) is None


def test_find_best_dilution_uses_min_non_excluded_positive() -> None:
    counts = {
        0: Counter({"0": 3}),
        -1: Counter({"1": 2}),
        -2: Counter({"1": 3}),
    }
    assert find_best_dilution_any_lysis(counts) == -2


def test_potency_helpers_zero_for_excluded_dilution() -> None:
    assert potency_rank_from_dilution(-4) == 0
    assert potency_label_from_dilution(-4) == "none"


def test_regeneration_collapses_pair_with_only_minus_four_lysis_to_mlc_zero() -> None:
    raw_rows = [
        {"bacteria": "B1", "phage": "P1", "log_dilution": "0", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-1", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-2", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-4", "score": "1"},
    ]
    pair_counts = build_dilution_counts_by_pair(raw_rows)
    assert compute_mlc_for_pair(pair_counts[("B1", "P1")]) == 0


def test_regeneration_collapses_old_mlc_four_to_three_when_lysis_also_at_minus_two() -> None:
    # Old pipeline returned MLC=4 (min positive = -4); new pipeline must return MLC=3 (min = -2).
    raw_rows = [
        {"bacteria": "B1", "phage": "P1", "log_dilution": "0", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-1", "score": "0"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-2", "score": "1"},
        {"bacteria": "B1", "phage": "P1", "log_dilution": "-4", "score": "1"},
    ]
    pair_counts = build_dilution_counts_by_pair(raw_rows)
    assert compute_mlc_for_pair(pair_counts[("B1", "P1")]) == 3


def test_regeneration_yields_full_mlc_range_across_inputs() -> None:
    raw_rows = [
        # MLC=0 pair — no lysis anywhere
        {"bacteria": "B1", "phage": "P0", "log_dilution": "0", "score": "0"},
        {"bacteria": "B1", "phage": "P0", "log_dilution": "-2", "score": "0"},
        # MLC=1 pair — lysis at highest titre only
        {"bacteria": "B1", "phage": "P1", "log_dilution": "0", "score": "1"},
        # MLC=2 pair — lysis at middle titre (and not at 0)
        {"bacteria": "B1", "phage": "P2", "log_dilution": "0", "score": "0"},
        {"bacteria": "B1", "phage": "P2", "log_dilution": "-1", "score": "1"},
        # MLC=3 pair — lysis at lowest replicated titre only
        {"bacteria": "B1", "phage": "P3", "log_dilution": "0", "score": "0"},
        {"bacteria": "B1", "phage": "P3", "log_dilution": "-1", "score": "0"},
        {"bacteria": "B1", "phage": "P3", "log_dilution": "-2", "score": "1"},
    ]
    pair_counts = build_dilution_counts_by_pair(raw_rows)
    assert compute_mlc_for_pair(pair_counts[("B1", "P0")]) == 0
    assert compute_mlc_for_pair(pair_counts[("B1", "P1")]) == 1
    assert compute_mlc_for_pair(pair_counts[("B1", "P2")]) == 2
    assert compute_mlc_for_pair(pair_counts[("B1", "P3")]) == 3
