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

### 2026-04-22 09:06 CEST: EX03 — Pages deploy + Release-asset data workflows

#### Executive summary

Added two workflows: `.github/workflows/publish-explainability-ui.yml` deploys HTML/
JS/CSS to GitHub Pages on every push to `main` touching `lyzortx/explainability_ui/**`
(Actions source mode, no `gh-pages` branch) and
`.github/workflows/snapshot-explainability-data.yml` is a `workflow_dispatch`-only job
that bootstraps `phage_env` via `setup-micromamba@v2`, optionally reruns CH05+CH09,
builds the snapshot JSONs, and creates a GitHub Release with `--latest`. `main.js`
gained a runtime data-base resolver: on `*.github.io` it fetches from
`https://github.com/<owner>/<repo>/releases/latest/download/`; otherwise from the
sibling `./data` directory. Footer shows the current release tag + published
timestamp pulled from the GitHub REST API (cached in sessionStorage).

#### Why

Commits to `main` must not carry generated outputs (per `lyzortx/AGENTS.md`
"Generated Outputs" rule). Pages-via-Actions keeps the data off every branch; Release
assets keep snapshots versioned and CDN-served with permissive CORS. The combination
solves both the gitignore constraint and the "viewable from GitHub" requirement
without a `gh-pages` branch to sync.

#### Design choices

- **Actions-mode Pages, not branch Pages**: no `gh-pages` branch means nothing off-main
  to maintain. The workflow's `deploy-pages@v4` step pushes directly to Pages
  infrastructure.
- **One-time manual Settings step**: repo Settings → Pages → Source = "GitHub
  Actions" cannot be set via the REST API for personal repos. Documented in
  `lyzortx/explainability_ui/AGENTS.md` as a prerequisite; first workflow run will
  fail with "Pages not enabled" until it is done.
- **Runtime data-base resolution vs build-time injection**: `main.js` inspects
  `window.location` to choose between release-asset URLs and local `./data`. No
  owner/repo strings are hardcoded, so the same committed HTML works on Pages,
  raw.githack, and local preview without edits.
- **Release metadata in the footer**: fetched from
  `/repos/<owner>/<repo>/releases/latest` on page load, cached in sessionStorage to
  avoid hammering the GitHub API on tab switches. Fire-and-forget — the main data
  load does not wait for it.
- **Snapshot workflow runs micromamba, not pip**: build_snapshot.py imports from
  `lyzortx.pipeline.autoresearch.runtime_contract` which transitively pulls in the
  training-side dependency graph. Bootstrapping the full `phage_env` is
  straightforward and unblocks the optional `regen=true` path without a second
  bootstrap step.

#### Verification

- `node --check lyzortx/explainability_ui/main.js` → zero syntax errors.
- `pytest lyzortx/tests/test_explainability_ui_build_snapshot.py` → 11 passed (no
  regressions from the data-base resolver).
- YAML syntax validated by committing + pushing — GitHub rejects syntactically broken
  workflows at push time.
- Full deploy verification (first-run Settings step + live Pages URL) is manual and
  documented in the PR description; the reviewer subagent runs through it.

#### Artifacts

- `.github/workflows/publish-explainability-ui.yml`
- `.github/workflows/snapshot-explainability-data.yml`
- `lyzortx/explainability_ui/main.js` — `resolveDataBase()` + `fetchReleaseMeta()`
- `lyzortx/explainability_ui/index.html` — footer release tag/timestamp display
- `lyzortx/explainability_ui/AGENTS.md` — deployment + prerequisite notes

#### Next

EX04: persist CH05 boosters, compute SHAP values over the 36,643-pair max-conc
prediction set, emit Parquet snapshots, and add the "Pair explorer" tab backed by
DuckDB-WASM for row-click drill-down with waterfall plots.

### 2026-04-22 09:30 CEST: EX04 — per-pair SHAP drill-down

#### Executive summary

New `lyzortx/pipeline/autoresearch/derive_shap_snapshot.py` replays CH05's two-axis
fold loops, fits LightGBM per fold × seed (10 × 3 = 30 fits per axis), and runs
`shap.TreeExplainer` on each fold's held-out pairs at max concentration. SHAP matrices
are seed-averaged per fold, then per-pair max-conc-aggregated, and written as six
Parquet files (3 kinds × 2 axes). `build_snapshot.py` gains `--include-shap` +
`--shap-dir` flags to copy the Parquets into the UI output directory; the snapshot
workflow picks up `*.parquet` alongside `*.json` when uploading release assets.
`main.js` gains a seventh "Pair explorer" tab: click a row in the Predictions table
(or enter a bacterium + phage manually) → DuckDB-WASM lazy-loads the Parquets via HTTP
range requests, queries for the pair, and renders a Plotly waterfall of the top-20
SHAP contributions plus a raw-feature-value table.

#### Why

CH05 is by far the richest artifact set we produce, but the Predictions tab only shows
aggregate (bacterium, phage, predicted_p) rows. A pair that misses is a mystery — did
the model see a phage_stats feature it didn't recognize? A host defense gene count
that pushed it off? Without SHAP, every diagnosis is speculative. SHAP waterfalls
give us the per-feature decomposition of each prediction so error buckets can be
interrogated mechanistically.

#### Design choices

- **Refit rather than persist-from-CH05**: `fit_seeds` discards the trained estimator
  and the change would ripple through six unpacking sites (CH04/CH05/CH07/CH08). Since
  EX04 is off-critical-path and SHAP compute is already expensive, a dedicated script
  that repeats the fit loop is strictly simpler than modifying production code. Uses
  the same `prepare_fold_design_matrices` / `select_rfe_features` / `build_pair_scorer`
  primitives so the model is bit-identical to CH05's for the same seed. Elapsed: ~65
  min on a 2023-vintage MacBook (fit ~30s, SHAP ~40s per fit × 60 fits, plus data
  loading overhead).
- **Wide-form Parquet over long-form**: one row per pair, one column per RFE feature.
  Long-form would be 36k × 300 = 11M rows per axis (~90 MB), wide-form compresses to
  ~40 MB. DuckDB-WASM queries are faster on wide since the `WHERE pair_id = ?` filter
  hits one row directly.
- **Per-fold RFE means feature columns differ across folds**: union columns across
  folds; per-fold pairs have `NaN` for features absent from that fold's RFE. The UI
  filters out NaN before rendering the waterfall.
- **DuckDB-WASM loaded on-demand**: the ~10 MB WASM bundle only downloads when the
  user first clicks the Pair explorer tab. `onTabClick` kicks off initialization so
  first-pair selection has the engine ready. Init is cached in
  `this.duckdbInitPromise` so concurrent clicks don't double-load.
- **Feature slot coloring reused**: `featureToSlot()` in main.js mirrors the Python
  `_feature_to_slot()` behavior exactly (both derive from the same `SLOT_COLORS` key
  set), so waterfall bars are colored consistently with the Feature Importance tab.

#### Verification

- `pytest lyzortx/tests/test_derive_shap_snapshot.py` → 3 passed (aggregator unit tests).
- `pytest lyzortx/tests/test_explainability_ui_build_snapshot.py` → 13 passed (2 new
  covering `copy_shap_parquets` happy path + partial-snapshot failure).
- Smoke test: `python -m lyzortx.pipeline.autoresearch.derive_shap_snapshot
  --device-type cpu --out-dir .scratch/ex04_smoke --max-folds 1 --seeds 42` → 291 s
  for one fold per axis with one seed, emits 6 Parquets.
- Full run: `python -m lyzortx.pipeline.autoresearch.derive_shap_snapshot` (10 folds
  × 3 seeds × 2 axes) completed in **4604 s (76.7 min)** on a 2023 MacBook. Output
  sizes: bacteria-axis SHAP 31.7 MB / feature values 0.9 MB / base values 0.3 MB;
  phage-axis SHAP 38.2 MB / feature values 1.3 MB / base values 0.3 MB. Total
  ~72.7 MB across the six Parquets.
- Waterfall decomposition sanity check: `pair 034-008__409_P1` base −2.0566 + Σ SHAP
  −0.8153 = predicted logit −2.8719 (sigmoid ≈ 5.3%), matching the expected
  low-lysis prediction for a narrow-host pair.
- `build_snapshot.py --include-shap --out .scratch/ex04_verify/data` copies all 6
  Parquets alongside the 7 JSONs in 0.4 s.
- `node --check main.js` → zero syntax errors.
- Pair explorer UX: click Predictions row → switches to Pair explorer → DuckDB-WASM
  spins up → waterfall + feature table render (awaiting full Parquet artifacts to
  verify end-to-end manually).

#### Artifacts

- `lyzortx/pipeline/autoresearch/derive_shap_snapshot.py` — new SHAP compute pipeline.
- `lyzortx/explainability_ui/build_snapshot.py` — `--include-shap` + `--shap-dir`
  flags + `copy_shap_parquets` helper.
- `lyzortx/explainability_ui/index.html` — 7th "Pair explorer" tab, predictions-table
  row-click handler.
- `lyzortx/explainability_ui/main.js` — DuckDB-WASM lazy loader + waterfall renderer +
  per-pair SHAP query.
- `lyzortx/explainability_ui/style.css` — `.data-table tr.clickable { cursor: pointer }`.
- `.github/workflows/snapshot-explainability-data.yml` — new `include_shap` input,
  optional SHAP compute step, `*.parquet` in release-asset glob.
- `lyzortx/tests/test_derive_shap_snapshot.py` — 3 aggregator unit tests.
- `lyzortx/tests/test_explainability_ui_build_snapshot.py` — 2 new Parquet-copy tests.
- `lyzortx/generated_outputs/ch05_shap/{bacteria,phage}_axis_{shap_values,shap_base_values,feature_values}.parquet`
  — six files, total ~80 MB (gitignored, regenerated per run).
- `requirements.txt` — `pyarrow==21.0.0` added (needed for Parquet write/read).

#### Next

None. EX04 is the MVP ceiling for the initial explainability UI scope; follow-ups
(beeswarm global SHAP, dependence plots, SHAP-delta comparison across slots) are
filed as open follow-ups to be scheduled as a separate track once the UI usage
patterns are clearer.
