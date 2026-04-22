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

### 2026-04-22 08:53 CEST: EX02 — static UI scaffold + snapshot extractor

#### Executive summary

Built the committed UI scaffold under `lyzortx/explainability_ui/` (`index.html` +
`main.js` + `style.css`) and a `build_snapshot.py` extractor that reads CH05/CH09
artifacts and emits seven normalized JSONs. Plotly.js + Alpine.js loaded via CDN — no
npm, no bundler, no build step. Six tabs render: Overview (AUC/Brier lollipops with
CIs + stat cards), Cross-source (Guelin vs BASEL bars + table), Calibration
(reliability diagrams per axis with raw/calibrated overlays + per-variant ECE),
Feature importance (top-30 horizontal bar colored by slot), Per-slot breakdown
(small-multiples sorted by cumulative importance), Predictions table (filterable +
sortable, 36,643 rows per axis). Eleven unit tests cover parsers, assemblers, and an
end-to-end integration pass against the live CH05/CH09 artifacts.

#### Why

EX01 gave the UI its missing input (per-axis feature importance). EX02 delivers the
extractor + renderer so a developer can inspect the CHISEL model end-to-end from one
HTML page. Overview-only MVP; per-pair SHAP drill-down is deferred to EX04.

#### Design choices

- **Two concerns separated**: UI code (HTML/JS/CSS — committed, safe to edit) vs data
  snapshots (JSON under `.scratch/` — generated, never committed). The committed side
  follows `lyzortx/AGENTS.md` generated-outputs policy; the generated side ships as
  GitHub Release assets in EX03.
- **No build step**: Plotly + Alpine via CDN, vanilla JS everywhere else. The
  intended dev workflow is "open the file in an editor and reload the page"; keeping
  that property forever is more valuable than code-organization gains from a bundler.
- **Bact-axis cross-source CIs left null**: CH05 bootstraps phage-axis cross-source
  but not bact-axis (the paired bootstrap is phage-axis-scoped by design). Rather
  than fabricate CIs, the snapshot carries CH09's point estimates with explicit
  `null` CI bounds; the UI renders them with no error bars and a clear label.
- **Slot coloring**: features prefixed `{slot}__` map to a fixed palette derived from
  `runtime_contract.SLOT_SPECS` plus two pseudo-slots (`pair_depo_capsule`,
  `pair_receptor_omp`) and `pair_concentration`. Features without a recognized prefix
  fall through to "other".

#### Verification

- `pytest lyzortx/tests/test_explainability_ui_build_snapshot.py` → 11 passed
  (10 parser/assembler unit tests + 1 end-to-end integration test that skips when
  CH05/CH09 artifacts are absent).
- `python -m lyzortx.explainability_ui.build_snapshot --out .scratch/explainability_ui/data`
  → 7 JSONs materialize in 0.4 s (total 11.5 MB; gzip ~2 MB over the wire). Headline
  numbers verified: bact-axis AUC 0.807921 [0.793432, 0.822287], phage-axis AUC
  0.887042 [0.865808, 0.905509], 36,643 bact pairs.
- `python -m http.server -d .scratch/explainability_ui 8765` serves HTML + data JSONs;
  `curl` against `/` and `/data/ch05_summary.json` both return 200 with correct
  content. `node --check main.js` reports zero syntax errors.
- `ruff check` clean on all new code.

#### Artifacts

- `lyzortx/explainability_ui/{index.html, main.js, style.css, build_snapshot.py,
  AGENTS.md, CLAUDE.md}` — committed UI scaffold + extractor.
- `lyzortx/tests/test_explainability_ui_build_snapshot.py` — 11 tests.
- `.scratch/explainability_ui/data/*.json` — generated snapshot (regenerated per run,
  gitignored).

#### Next

EX03: wire the static site to GitHub Pages (Actions source mode) + publish data
JSONs as GitHub Release assets keyed to the `latest` tag; `main.js` fetches from the
release URL, owner/repo resolved from `window.location` at runtime.
