# CH13 Arm 3 canonical `phage_projection` migration

Single-shot scripts to materialize the Arm 3 Moriniere per-receptor-class k-mer
fractions as the canonical `phage_projection` slot. Run these once, in order,
after regenerating the autoresearch cache on a fresh clone (or whenever the
CH06 Arm 3 source CSV at `.scratch/basel/feature_slots_arm3/phage_projection/features.csv`
is regenerated).

## Usage

```
python lyzortx/research_notes/ad_hoc_analysis_code/ch13_arm3_migration/migrate_slot_files.py
python lyzortx/research_notes/ad_hoc_analysis_code/ch13_arm3_migration/patch_top_level_schema.py
```

## What they do

`migrate_slot_files.py`:

1. Moves current `phage_projection/` slot artifacts aside to
   `phage_projection_tl17/` under both `.scratch/basel/feature_slots/` (unified
   148-phage panel) and `lyzortx/generated_outputs/autoresearch/search_cache_v1/feature_slots/`
   (Guelin-only 96-phage cache).
2. Materializes the Arm 3 per-receptor-class fraction slot (13 `phage_projection__recep_frac_*`
   columns) at the canonical path in both locations, filtering the 148-phage
   source to the 96 Guelin phages for the autoresearch cache.
3. Rewrites per-slot `schema_manifest.json` on both paths.

`patch_top_level_schema.py`:

1. Updates the top-level `ar02_schema_manifest_v1.json`'s
   `feature_slots.phage_projection.reserved_feature_columns` list from the 33
   TL17 family-presence column names to the 13 Arm 3 receptor-class
   column names. Without this, `lyzortx/autoresearch/train.py::load_slot_artifact`
   raises `ValueError: slot phage_projection bypassed the frozen top-level
   cache schema` at CH04 startup.

## Why this lives in ad_hoc_analysis_code/

The autoresearch cache and `.scratch/` slots are both gitignored (see root
`AGENTS.md` "Generated outputs" rule), so the migration cannot be committed
as an artifact. It has to be reproducible via scripts that agents run after
regenerating the cache. These two scripts are idempotent: re-running them on
a clone that already has the migration applied is a no-op (the TL17 side-dirs
just already exist and the per-slot schema is already Arm 3).

## Follow-up

`prepare.py` (the only path from raw inputs to the search cache) should
eventually be rewired to produce the Arm 3 phage_projection natively, at which
point these scripts can be deleted. Until then, a fresh-clone bootstrap is
`prepare.py` + these two scripts, in that order.
