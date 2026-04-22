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

If a new library is needed, pin its version in the `<script>` tag and document it
here. Do NOT introduce a build step without explicit user approval — "open the file and
edit it" is the intended dev workflow.

## Deployment (EX03)

- HTML/JS/CSS deploy to GitHub Pages via `actions/deploy-pages@v4` (GitHub Actions source
  mode — no `gh-pages` branch is ever created).
- Data snapshots publish as GitHub Release assets. Release tag
  `explainability-data-YYYYMMDD-HHMM`; the HTML fetches from
  `https://github.com/<owner>/<repo>/releases/latest/download/<file>.json` so the page
  doesn't need to know the tag.
- Prerequisite one-time step: repo Settings → Pages → Source = "GitHub Actions".
  GitHub does not expose this via API on personal repos.

## Policies

- Do not commit any generated JSON under this directory or under
  `lyzortx/explainability_ui/data/` — follow the `lyzortx/AGENTS.md` "Generated Outputs"
  rule.
- Do not hardcode owner/repo strings in `main.js` — resolve from `window.location` at
  runtime so the same HTML works on Pages, raw.githack, and local preview.
