"""SX06: Unit tests for the BASEL TL17 slot-merge helper.

The full projection path requires pyrodigal + mmseqs and 52 BASEL FASTAs, so the test
exercises only the deterministic CSV-merge logic on a synthetic fixture.
"""

from __future__ import annotations

from pathlib import Path


from lyzortx.pipeline.autoresearch.project_basel_tl17_features import (
    SLOT_HIT_COUNT_COLUMN,
    merge_basel_features_into_slot,
)


def test_merge_overwrites_basel_rows_and_preserves_guelin(tmp_path: Path) -> None:
    slot_path = tmp_path / "features.csv"
    slot_path.write_text(
        "\n".join(
            [
                "phage,phage_projection__tl17_phage_rbp_family_1_present,"
                f"phage_projection__tl17_phage_rbp_family_2_present,{SLOT_HIT_COUNT_COLUMN}",
                "GuelinA,42.0,0.0,3",
                "GuelinB,0.0,17.0,1",
                "Bas01,0.0,0.0,0",
                "Bas02,0.0,0.0,0",
                "",
            ]
        )
    )

    basel_feature_rows = [
        {
            "phage": "Bas01",
            "tl17_phage_rbp_family_1_present": 88.5,
            "tl17_phage_rbp_family_2_present": 0.0,
            "tl17_rbp_reference_hit_count": 2,
        },
        {
            "phage": "Bas02",
            "tl17_phage_rbp_family_1_present": 0.0,
            "tl17_phage_rbp_family_2_present": 71.2,
            "tl17_rbp_reference_hit_count": 4,
        },
    ]

    merged = merge_basel_features_into_slot(slot_path, basel_feature_rows)

    # Guelin rows unchanged
    guelin_a = merged[merged["phage"] == "GuelinA"].iloc[0]
    assert guelin_a["phage_projection__tl17_phage_rbp_family_1_present"] == 42.0
    assert guelin_a["phage_projection__tl17_phage_rbp_family_2_present"] == 0.0
    assert guelin_a[SLOT_HIT_COUNT_COLUMN] == 3

    # BASEL rows now carry projected values
    bas01 = merged[merged["phage"] == "Bas01"].iloc[0]
    assert bas01["phage_projection__tl17_phage_rbp_family_1_present"] == 88.5
    assert bas01[SLOT_HIT_COUNT_COLUMN] == 2
    bas02 = merged[merged["phage"] == "Bas02"].iloc[0]
    assert bas02["phage_projection__tl17_phage_rbp_family_2_present"] == 71.2
    assert bas02[SLOT_HIT_COUNT_COLUMN] == 4


def test_merge_errors_on_unknown_basel_phage(tmp_path: Path) -> None:
    slot_path = tmp_path / "features.csv"
    slot_path.write_text(
        f"phage,phage_projection__tl17_phage_rbp_family_1_present,{SLOT_HIT_COUNT_COLUMN}\nBas01,0.0,0\n"
    )
    basel_feature_rows = [
        {"phage": "BasUNKNOWN", "tl17_phage_rbp_family_1_present": 5.0, "tl17_rbp_reference_hit_count": 1}
    ]
    try:
        merge_basel_features_into_slot(slot_path, basel_feature_rows)
    except KeyError as exc:
        assert "BasUNKNOWN" in str(exc)
        return
    raise AssertionError("Expected KeyError for unknown BASEL phage not in slot CSV")


def test_merge_errors_on_unknown_runtime_column(tmp_path: Path) -> None:
    slot_path = tmp_path / "features.csv"
    slot_path.write_text("phage,phage_projection__tl17_phage_rbp_family_1_present\nBas01,0.0\n")
    basel_feature_rows = [{"phage": "Bas01", "tl17_phage_rbp_family_UNKNOWN_present": 7.0}]
    try:
        merge_basel_features_into_slot(slot_path, basel_feature_rows)
    except KeyError as exc:
        assert "UNKNOWN" in str(exc)
        return
    raise AssertionError("Expected KeyError when runtime column has no matching slot column")
