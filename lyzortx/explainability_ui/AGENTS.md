# Explainability UI Directory

## Purpose

Static HTML + Plotly.js dashboard over the CHISEL CH05 unified baseline. Overview-only
MVP: headline AUC/Brier per axis, Guelin-vs-BASEL cross-source, CH09 reliability diagrams,
feature importance, per-slot breakdown, and a filterable predictions table.

## Separation of concerns

- **UI code** (`index.html`, `main.js`, `style.css`) — committed, never generated. Safe
  to edit in place.
- **Snapshot extractor** (`build_snapshot.py`) — committed. Reads gitignored CH05/CH09
  artifacts under `lyzortx/generated_outputs/` and writes normalized JSONs.
- **Snapshot data** (JSON files) — generated, never committed. Default output is
  `.scratch/explainability_ui/data/` for local dev; in production (per EX03) they ship
  as GitHub Release assets and the HTML fetches from the `latest` release URL.

## Running locally

```bash
python -m lyzortx.explainability_ui.build_snapshot --out .scratch/explainability_ui/data
cp lyzortx/explainability_ui/{index.html,main.js,style.css} .scratch/explainability_ui/
python -m http.server -d .scratch/explainability_ui 8765
```

Then open `http://localhost:8765/`. `fetch('./data/*.json')` requires same-origin, so
`file://` opens do not work — always serve via `python -m http.server`.

## Dependencies

Only CDN-loaded JS libraries — no npm, no bundler, no build step. Version pinned in
`index.html`:

- Plotly.js 2.35.2 (via `cdn.plot.ly`)
- Alpine.js 3.14.1 (via `unpkg.com`)
- DuckDB-WASM 1.29.0 (via `cdn.jsdelivr.net`, dynamically imported from `main.js`
  only when the Pair explorer tab is first opened — keeps the ~10 MB WASM bundle out
  of the initial page load for users who don't drill down)

If a new library is needed, pin its version in the `<script>` tag (or `import()`
specifier) and document it here. Do NOT introduce a build step without explicit user
approval — "open the file and edit it" is the intended dev workflow.

## Pair explorer (SHAP drill-down)

The 7th tab ("Pair explorer") runs on six SHAP Parquet assets that ship alongside the
JSON snapshots as GitHub Release assets:

- `{bacteria,phage}_axis_shap_values.parquet` — per-pair SHAP contributions
- `{bacteria,phage}_axis_feature_values.parquet` — matching raw feature values
- `{bacteria,phage}_axis_shap_base_values.parquet` — per-pair model log-odds base

Produced by `lyzortx/pipeline/autoresearch/derive_shap_snapshot.py` (see
`lyzortx/pipeline/CLAUDE.md` for pipeline conventions). Run this once per CH05
refresh; takes ~65 min wallclock. `build_snapshot.py --include-shap` copies the
Parquets into the deploy directory.

DuckDB-WASM queries these over HTTP range requests — only the row matching the
selected `(bacteria, phage)` pair is transferred, so latency is bounded by one HTTP
round-trip per pair regardless of Parquet size.

## Deployment

- HTML/JS/CSS **and the current release's data assets** deploy together to GitHub
  Pages via `.github/workflows/publish-explainability-ui.yml` (fires on pushes to
  `main` that touch `lyzortx/explainability_ui/**`, plus `workflow_dispatch`). Uses
  `actions/deploy-pages@v4` in "GitHub Actions" source mode — no `gh-pages` branch is
  ever created.
- The deploy workflow runs `gh release view --json tagName -q .tagName` to discover
  the `latest` tag, then `gh release download $tag --pattern '*.json' --pattern
  '*.parquet' --dir _site/explainability/data/`. The UI fetches from `./data/…`
  same-origin, bypassing GitHub's missing CORS headers on release-asset URLs.
- **Publishing a data snapshot**: run locally.

  ```bash
  python -m lyzortx.pipeline.autoresearch.ch05_eval --device-type cpu           # ~55 min
  python -m lyzortx.pipeline.autoresearch.ch09_calibration_layer                # ~5 min
  python -m lyzortx.pipeline.autoresearch.derive_shap_snapshot --device-type cpu  # ~75 min (optional, needed for the Pair explorer tab)
  python -m lyzortx.explainability_ui.build_snapshot --include-shap --out .scratch/release_snapshot
  tag="explainability-data-$(date -u +%Y%m%d-%H%M)"
  gh release create "$tag" --latest --title "Explainability data snapshot $tag" \
    .scratch/release_snapshot/*.json .scratch/release_snapshot/*.parquet
  ```

  Any push to `lyzortx/explainability_ui/**` after that re-bakes the fresh assets
  into the next Pages deploy. To force a redeploy without a code change, trigger
  `publish-explainability-ui.yml` via `workflow_dispatch`.
- **Prerequisite one-time manual step** (not automatable): repo Settings → Pages →
  Source = "GitHub Actions". GitHub does not expose this via API on personal repos.
  Without it, the first `publish-explainability-ui` run fails with a "Pages not
  enabled" error.

## Policies

- Do not commit any generated JSON under this directory or under
  `lyzortx/explainability_ui/data/` — follow the `lyzortx/AGENTS.md` "Generated Outputs"
  rule.
- Do not hardcode owner/repo strings in `main.js` — resolve from `window.location` at
  runtime so the same HTML works on Pages, raw.githack, and local preview.
