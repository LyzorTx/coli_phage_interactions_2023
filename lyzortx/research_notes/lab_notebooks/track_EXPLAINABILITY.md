# Track EXPLAINABILITY

## 2026-04-22 07:47 CEST: EX01 — per-axis feature importance emission

### Executive summary

Patched `ch05_eval.py` so both `run_bacteria_axis` and `run_phage_axis` persist
`ch05_{bacteria,phage}_axis_feature_importance.csv` with the same schema as
`ch04_feature_importance.csv` (`feature, is_concentration_feature, mean_importance,
n_folds_selected`, sorted descending). Extracted the aggregation as
`_aggregate_feature_importance` and added unit tests covering the schema contract,
descending order, and the empty-input failure path. Full CH05 rerun (bacteria-axis +
phage-axis) re-emitted all existing artifacts plus the two new CSVs. No model or
evaluation logic changed; this is a pure artifact-emission patch.

#### Why

The explainability UI (EX02) needs per-axis feature importance to render the top-30 bar
chart and per-slot breakdown tabs. CH05 was already training 30 boosters per axis
(10 folds × 3 seeds) and calling the same `fit_seeds` machinery CH04 uses, but it was
discarding each seed's `feature_importance_` frame (`_fi` underscore-prefix on the
seed_results unpacking). Persisting it was a 5-line change inside each axis function
plus a shared aggregation helper.

#### Design choices

- **Shared `_aggregate_feature_importance` helper, not axis-specific code**: both axes
  produce the same per-seed frames from `fit_seeds`; one helper keeps the CSV shape
  identical between axes and identical to CH04. The helper raises `ValueError` on empty
  input (CH04 does not; CH05 should — an empty aggregate means zero fits ran, which is
  a pipeline bug, not a recoverable state).
- **No change to `fit_seeds` or `_fit_one_seed`**: the `feature_importance` frame is
  already populated with `is_concentration_feature` in `ch04_parallel.py:230`, so CH05
  just needs to stop discarding it.
- **Drop the CH04 pattern of attaching `fold_id`/`seed` to the per-fit frame**: CH04
  does `fi["fold_id"] = fold_id; fi["seed"] = seed` before appending, but its aggregate
  groups only on `(feature, is_concentration_feature)` so the extra columns are
  never used. Mirroring CH04 exactly preserves the option to emit long-form frames in a
  future ticket, so I kept the fold/seed tagging even though it's unused at aggregate
  time.

#### Verification

- Unit tests: `pytest lyzortx/tests/test_ch05_unified_kfold.py` → 10 passed
  (2 new covering the aggregator).
- Ruff: `ruff check lyzortx/pipeline/autoresearch/ch05_eval.py
  lyzortx/tests/test_ch05_unified_kfold.py` → clean.
- CH05 rerun: `python -m lyzortx.pipeline.autoresearch.ch05_eval --device-type cpu`.
  Elapsed + headline numbers + CSV row counts recorded in the PR description once the
  run completes.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch05_eval.py` — aggregation helper + axis functions
  now emit the per-axis feature-importance CSV.
- `lyzortx/tests/test_ch05_unified_kfold.py` — two new tests on the aggregator.
- `lyzortx/generated_outputs/ch05_unified_kfold/ch05_bacteria_axis_feature_importance.csv`
  — new canonical artifact (gitignored, regenerated per run).
- `lyzortx/generated_outputs/ch05_unified_kfold/ch05_phage_axis_feature_importance.csv`
  — new canonical artifact (gitignored, regenerated per run).

#### Next

EX02: wire the two CSVs into the UI snapshot extractor and render them in a top-30
feature-importance bar chart + per-slot breakdown small-multiples.
