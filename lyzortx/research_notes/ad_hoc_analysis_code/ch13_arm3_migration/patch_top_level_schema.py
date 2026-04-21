"""Patch the top-level autoresearch cache schema manifest to point phage_projection
at the Arm 3 per-receptor-fraction columns instead of TL17 BLAST families.
"""

from __future__ import annotations

import json
from pathlib import Path

MANIFEST_PATH = Path("lyzortx/generated_outputs/autoresearch/search_cache_v1/ar02_schema_manifest_v1.json")

ARM3_COLS = [
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


def main() -> None:
    with open(MANIFEST_PATH, encoding="utf-8") as f:
        manifest = json.load(f)

    slot = manifest["feature_slots"]["phage_projection"]
    prev_cols = slot["reserved_feature_columns"]
    print(
        f"Top-level phage_projection had {len(prev_cols)} TL17 columns; replacing with {len(ARM3_COLS)} Arm 3 columns."
    )
    slot["reserved_feature_column_count"] = len(ARM3_COLS)
    slot["reserved_feature_columns"] = ARM3_COLS
    # Preserve block_role / column_family_prefix / entity_key / join_keys.

    with open(MANIFEST_PATH, "w", encoding="utf-8") as f:
        json.dump(manifest, f, indent=2)
    print(f"Wrote {MANIFEST_PATH}")


if __name__ == "__main__":
    main()
