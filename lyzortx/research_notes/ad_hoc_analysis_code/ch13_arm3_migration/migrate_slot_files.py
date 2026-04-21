"""CH13: migrate Arm 3 per-receptor-fraction slot to canonical phage_projection.

Touches two cache dirs:
  (a) .scratch/basel/feature_slots/                    -- unified 148-phage panel (CH05/CH07/CH08)
  (b) lyzortx/generated_outputs/autoresearch/search_cache_v1/feature_slots/
      -- Guelin-only 96-phage cache (CH04 reads this directly)

For each: move current phage_projection/ to phage_projection_tl17/ (sensitivity fallback),
then materialize a new phage_projection/ from Arm 3's source CSV. For (b), filter Arm 3's
148 phages to the 96 that appear in the cache's entity_index. Rewrite schema_manifest.json.
"""

from __future__ import annotations

import json
import shutil
from pathlib import Path

import pandas as pd

ARM3_SRC = Path(".scratch/basel/feature_slots_arm3/phage_projection")
UNIFIED_SLOT = Path(".scratch/basel/feature_slots/phage_projection")
UNIFIED_SLOT_TL17 = Path(".scratch/basel/feature_slots/phage_projection_tl17")
CACHE_SLOT = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1/feature_slots/phage_projection")
CACHE_SLOT_TL17 = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1/feature_slots/phage_projection_tl17")

ARM3_FEATURE_COLS = [
    "phage_projection__recep_frac_GluI",
    "phage_projection__recep_frac_HepI",
    "phage_projection__recep_frac_HepII",
    "phage_projection__recep_frac_Kdo",
    "phage_projection__recep_frac_NGR",
    "phage_projection__recep_frac_btuB",
    "phage_projection__recep_frac_fhuA",
    "phage_projection__recep_frac_lamB",
    "phage_projection__recep_frac_lptD",
    "phage_projection__recep_frac_ompA",
    "phage_projection__recep_frac_ompC",
    "phage_projection__recep_frac_ompF",
    "phage_projection__recep_frac_tsx",
]


def move_tl17_aside(src: Path, dst: Path) -> None:
    if dst.exists():
        print(f"  {dst} already exists; skipping move")
        return
    if src.exists():
        print(f"  moving {src} -> {dst}")
        shutil.move(str(src), str(dst))


def materialize_arm3_slot(target_dir: Path, arm3_features_csv: Path, entity_index_restrict: list[str] | None) -> None:
    """Write Arm 3 features.csv + schema_manifest.json + entity_index.csv into target_dir."""
    target_dir.mkdir(parents=True, exist_ok=True)
    df = pd.read_csv(arm3_features_csv)
    expected_cols = ["phage"] + ARM3_FEATURE_COLS
    if list(df.columns) != expected_cols:
        raise ValueError(f"Arm 3 source columns mismatch.\n  expected: {expected_cols}\n  found: {list(df.columns)}")
    if entity_index_restrict is not None:
        missing = set(entity_index_restrict) - set(df["phage"])
        if missing:
            raise ValueError(f"Arm 3 source missing phages required by cache entity_index: {sorted(missing)[:5]}...")
        df = df[df["phage"].isin(entity_index_restrict)].copy()
        df = df.set_index("phage").loc[entity_index_restrict].reset_index()
    features_out = target_dir / "features.csv"
    df.to_csv(features_out, index=False)
    print(f"  wrote {features_out}: {len(df)} phages, {len(df.columns)} cols")

    manifest = {
        "block_role": "phage",
        "cache_contract_id": "autoresearch_search_cache_v1",
        "column_family_prefix": "phage_projection__",
        "composability_contract": {
            "column_ownership": (
                "Future columns for phage_projection must start with phage_projection__ "
                "and may only be added inside this slot."
            ),
            "join_type": "left",
            "row_granularity": "one_row_per_phage",
        },
        "description": (
            "Phage projection features: Moriniere 2026 per-receptor-class k-mer fractions. "
            "Each column is |kmers(R) ∩ kmers(P)| / |kmers(R)| for receptor class R. "
            "Panel-independent (trained on 260 non-Guelin reference phages). Canonical "
            "since CH13 (2026-04-21); supersedes TL17 BLAST phage_rbp_family presence "
            "vectors, which are retained at phage_projection_tl17/ as an opt-in "
            "sensitivity fallback."
        ),
        "entity_index_row_count": len(df),
        "entity_key": "phage",
        "join_keys": ["phage"],
        "reserved_feature_column_count": len(ARM3_FEATURE_COLS),
        "reserved_feature_columns": ARM3_FEATURE_COLS,
        "schema_manifest_id": "autoresearch_feature_schema_v1",
        "slot_name": "phage_projection",
        "task_id": "CH13",
    }
    manifest_out = target_dir / "schema_manifest.json"
    with open(manifest_out, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"  wrote {manifest_out}")

    entity_index_out = target_dir / "entity_index.csv"
    df[["phage"]].to_csv(entity_index_out, index=False)
    print(f"  wrote {entity_index_out}: {len(df)} rows")


def main() -> None:
    print("=== CH13 migration: Arm 3 -> canonical phage_projection ===\n")

    print("(1) Unified panel slot (.scratch/basel/feature_slots/)")
    move_tl17_aside(UNIFIED_SLOT, UNIFIED_SLOT_TL17)
    materialize_arm3_slot(
        target_dir=UNIFIED_SLOT,
        arm3_features_csv=ARM3_SRC / "features.csv",
        entity_index_restrict=None,  # keep all 148 phages
    )

    print("\n(2) Autoresearch cache slot (search_cache_v1/feature_slots/)")
    # Read cache's existing entity_index (96 Guelin phages) before moving.
    cache_entity_csv = CACHE_SLOT / "entity_index.csv" if CACHE_SLOT.exists() else CACHE_SLOT_TL17 / "entity_index.csv"
    cache_phages = pd.read_csv(cache_entity_csv)["phage"].tolist()
    print(f"  cache entity_index has {len(cache_phages)} phages")
    move_tl17_aside(CACHE_SLOT, CACHE_SLOT_TL17)
    materialize_arm3_slot(
        target_dir=CACHE_SLOT,
        arm3_features_csv=ARM3_SRC / "features.csv",
        entity_index_restrict=cache_phages,
    )

    print("\n=== Migration complete ===")


if __name__ == "__main__":
    main()
