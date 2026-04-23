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

### Hosting and data flow

- Do not commit any generated JSON under this directory or under
  `lyzortx/explainability_ui/data/` — follow the `lyzortx/AGENTS.md` "Generated Outputs"
  rule.
- Do not hardcode owner/repo strings in `main.js` — resolve from `window.location` at
  runtime so the same HTML works on Pages, raw.githack, and local preview.
- **GitHub Release asset URLs do not serve `Access-Control-Allow-Origin`.** Do not
  `fetch()` from `github.com/.../releases/latest/download/...` from the UI; the browser
  will block the request. The workflow bundles release assets into the Pages deploy
  artifact under `_site/explainability/data/` so the UI fetches same-origin via
  `./data/`. This was the EX05 bug (PR #473) — do not re-introduce runtime release URLs.
- **api.github.com DOES send `Access-Control-Allow-Origin: *`.** It's safe to fetch
  release metadata (the footer tag) from `api.github.com/repos/:owner/:repo/releases/
  latest` at runtime.

### URL-first state

- **The URL is the source of truth for all view state.** Every tab, axis toggle,
  filter, and pair selection must sync to `window.location.search` via the
  `URL_BINDINGS` table in `main.js`. Any reload or shared link must reconstruct the
  exact view. Do not introduce ad-hoc state that does not have a URL binding.
- Add new URL-backed state by appending one row to `URL_BINDINGS` (param name, state
  field, default, allowlist, tab gate). Do NOT open-code `new URLSearchParams(...)`
  reads or writes elsewhere in `main.js`.
- Defaults are omitted from URLs. A view at the default tab / default axes has a bare
  URL with no query string; shared links only carry the non-default bits.
- Use `history.replaceState`, not `pushState`, for filter/toggle changes. We don't
  spam back-button history for every axis click.

### DuckDB-WASM gotchas

- **Use absolute URLs, not relative, when the DuckDB worker reads Parquets.** Inside
  the worker, `self.location` is the worker's blob URL, so `./data/foo.parquet`
  resolves to `blob:...` and 404s. Compute absolute URLs via `new URL(rel,
  window.location.href).href` before passing to `read_parquet(...)`. This was the
  EX06 bug (PR #474).
- Prefer `read_parquet('<abs-url>')` inline in queries over the old
  `registerFileURL(...)` pattern. The registration API silently fails on 404s and
  produces confusing downstream "No files found that match the pattern" errors.

### Third-party library posture

- **marked.js does not sanitize HTML output and dropped `headerIds` in v13.** When
  rendering markdown into `x-html`, always:
  1. Register the `marked-gfm-heading-id` extension (`marked.use(gfmHeadingId())`)
     so heading anchors work for deep-links.
  2. Wrap the parse result in `DOMPurify.sanitize(...)` to prevent XSS from any
     untrusted content that could leak in via `/sleeponit` rewriting `KNOWLEDGE.md`
     from lab-notebook prose.
  Both deps are lazy-loaded from jsDelivr on first Glossary/Knowledge visit so
  initial page load stays light.

## Development workflow

### Local feedback loop beats deploy-and-check

- **Always preview locally before pushing.** The Pages deploy pipeline is ~2 min per
  iteration; a local `python3 -m http.server` iteration is seconds. Do not commit
  UI changes without serving them locally first.
- Canonical preview layout matches the Pages deploy structure:

  ```bash
  mkdir -p .scratch/ex_preview/explainability/{data,docs}
  cp lyzortx/explainability_ui/{index.html,main.js,style.css} .scratch/ex_preview/explainability/
  cp .scratch/release_snapshot/*.{json,parquet} .scratch/ex_preview/explainability/data/
  cp lyzortx/research_notes/GLOSSARY.md lyzortx/KNOWLEDGE.md .scratch/ex_preview/explainability/docs/
  python3 -m http.server 8767 --directory .scratch/ex_preview
  # → http://localhost:8767/explainability/
  ```

- Verify all URL routing permutations locally before pushing — `?tab=X`, `?tab=X&bact=
  Y&phage=Z`, anchor deep-links like `?tab=glossary#slot-feature-slot`.

### Snapshot provenance check

Before publishing a new release, verify the snapshot's AUC+Brier headline matches the
current canonical `chisel-unified-kfold-baseline` numbers in `lyzortx/KNOWLEDGE.md`
(the post-CH11 reference under the Arm 3 phage_projection slot):

```bash
python -c "
import json, pathlib
s = json.load(open('.scratch/release_snapshot/ch05_summary.json'))
for axis in ('bacteria', 'phage'):
    a = s['axes'][axis]
    print(f'{axis}: AUC {a[\"auc\"][\"point\"]:.4f} '
          f'[{a[\"auc\"][\"lo\"]:.4f}, {a[\"auc\"][\"hi\"]:.4f}], '
          f'Brier {a[\"brier\"][\"point\"]:.4f}')
"
```

Expected match (CH11 headline, rerun 2026-04-21):

- bacteria-axis AUC 0.8079 [0.7934, 0.8223], Brier 0.1763 [0.1688, 0.1840]
- phage-axis AUC 0.8870 [0.8658, 0.9055], Brier 0.1352 [0.1227, 0.1489]

If the snapshot point estimate falls outside the canonical CI, the snapshot was
computed on a pre-CH11 CH05 run — re-run the publish flow from scratch (CH05 → CH09
→ derive_shap_snapshot → build_snapshot). Drift within the CI is acceptable (stems
from concentration_mean_importance-level fluctuation in LightGBM RFE across reruns,
not canonical change).

When `chisel-unified-kfold-baseline` shifts (e.g. a future CH track retrains the
model), update the expected numbers above in the same PR that bumps the knowledge
unit — otherwise this check silently stops matching.

### Release publishing

- Releases publish from the local laptop via `gh release create`, not from CI. The
  upstream pipeline (ST01 → CH05 → derive_shap_snapshot) depends on gitignored
  `lyzortx/generated_outputs/` artifacts, which CI runners don't have. The
  `snapshot-explainability-data.yml` workflow was removed for exactly this reason
  (EX05).
- The Pages deploy fires automatically on any push to `lyzortx/explainability_ui/**`
  and re-downloads the latest release assets. Publishing a new release does NOT
  auto-redeploy — trigger `publish-explainability-ui.yml` via `workflow_dispatch` to
  force a fresh bundle of the new assets.

### Commit cadence

- UI iteration commits can go directly to `main` when the user explicitly invites it
  — UI changes don't have the orchestrator-ticket reviewer requirement from the root
  `AGENTS.md`. Split large changes into logically named commits (glossary entry,
  workflow update, UI bundle, cleanup refactor) so `git log --oneline` tells the
  change story.
- Spawn the self-review subagent only for logic-heavy commits (URL state machines,
  SHAP aggregation, DuckDB query refactors). Trivial edits (docs, workflow cp
  additions, single-line fixes) don't warrant it. Use judgement: "could a future
  reader reasonably misunderstand this?" → yes means review.
