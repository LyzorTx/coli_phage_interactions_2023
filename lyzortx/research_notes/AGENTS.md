# Research Notes Directory

## How lab notebooks are written

Lab notebooks live under `lab_notebooks/`. Each file is either track-specific (`track_ST.md`, `track_B.md`, etc.) or
cross-track (`project.md`).

Entries are written by agents dispatched through the orchestrator pipeline:

1. `lyzortx/orchestration/orchestrator.py` creates GitHub issues for pending plan tasks.
2. `.github/workflows/orchestrator.yml` triggers the orchestrator on issue close and workflow dispatch.
3. When an agent implements a task, the orchestrator issue body instructs it to write findings to
   `lyzortx/research_notes/lab_notebooks/track_{TRACK}.md` following the existing entry format.
4. `project.md` is used for cross-track strategic decisions and plan-level notes that are not specific to one track.
5. `devops.md` records coding infrastructure changes: CI/CD workflows, tooling, automation, token/cost analysis,
   developer experience improvements, and other non-science engineering decisions.

## When to write replanning entries

When the `/replan` skill is used and results in plan changes, write entries to both the affected track notebook and
`project.md`. The track entry should document the technical findings (e.g., specific leaked features, nondeterminism
evidence, metric comparisons). The project entry should document the strategic decision (e.g., which tracks were killed
or restructured, and why). Include source citations (URLs + quotes) for any claims about external library behavior.

## Track closeout recaps

When a track finishes (all its tickets are merged and no further work is planned under that track), add a
`lab_notebooks/track_{TRACK}_recap.md` file that summarizes what the track accomplished. The recap is separate from the
per-ticket entries in `track_{TRACK}.md`: the track file is the episodic log, the recap is the digest.

A recap should cover, in roughly this order:

- **Goal** the track was chartered to answer, in one or two sentences.
- **Headline outcomes** — the load-bearing numbers (e.g., new canonical baseline, deltas vs prior canonical, closure of
  any named gap) with CI where meaningful.
- **What changed in the canonical pipeline** — training policy, feature slots, evaluation protocol, scorecard metrics.
- **Dead ends and null arms** — what was tested and found not to lift, with one-line reasons (for later tracks to not
  re-litigate).
- **Open follow-ups** — concrete items that the track surfaced but deferred, so a future replan can pick them up.
- **Artifact pointers** — the canonical generated-outputs directories and any scripts that reproduce headline numbers.

Keep the recap short — a reader should be able to understand the track's net contribution in a few minutes without
reading every ticket entry. Link back to notebook entries for detail rather than duplicating them.

## In-flight entry editing

When working on a feature branch or open PR, treat notebook entries touched in that branch as mutable working
documents. If subsequent work in the same branch invalidates statements in an in-flight entry (e.g., a metric was
wrong, a method was replaced, the split contract changed), delete the stale statements and replace them with the final
correct ones. Do not preserve churn with "initially we thought X" hedging — the git history already records the
evolution. The entry as merged should read as a clean, accurate record of the final state.

This applies only to entries being developed in the current branch/PR. Entries that were already merged to `main` are
historical records and must not be modified (see Done Task Immutability in the orchestration AGENTS.md).

## Entry format

- Each entry starts with `### YYYY-MM-DD HH:MM TZ: Title` (date heading level 3, local timezone with zone label e.g.
  CEST; timestamp helps resolve merge conflicts when two entries are added on the same day).
- Entries within a track file are ordered by task code, earliest first.
- Every entry must begin with an `#### Executive summary` section: 2-4 sentences covering what changed, why, and the
  key outcome or metric. A reader should be able to skip the rest of the entry and still understand the decision.
- Subsequent sections typically include: problem statement, design decisions, interpretation, and next steps.
- Entries should reference generated output paths and script paths so findings are traceable.
- Do not list files changed — that is what git is for.
- **"Future:" notes** — When an agent discovers something worth revisiting later (a deferred cleanup, a tool adoption
  trigger, a feature idea), add a section to `project.md` with a heading starting with `#### Future:`. Include the
  trigger condition ("revisit when...") so the `/replan` skill can evaluate whether the condition is now met. These notes
  are seeds for future plan tasks — they should be concrete enough to act on, not vague aspirations.

## Relationship to the knowledge model

Lab notebooks are the **episodic** record — specific experiments, dates, metrics, and intermediate findings. The
knowledge model (`lyzortx/orchestration/knowledge.yml` → `lyzortx/KNOWLEDGE.md`) is the **semantic** distillation —
validated facts, dead-end lessons, and active assumptions extracted from notebooks via the `/sleeponit` skill.

- Notebooks are the source of truth for *what happened*. The knowledge model is the source of truth for *what we know*.
- When writing a notebook entry that invalidates existing knowledge (e.g., a feature that was "active" is now proven
  leaky), note the invalidation explicitly. This helps the next `/sleeponit` run update the knowledge model.
- Source references in the knowledge model (e.g., `[TG02]`) point back to notebook entries. Keep notebook entry IDs
  (task codes) stable so these references remain valid.

## GLOSSARY.md

`GLOSSARY.md` defines recurring terms — metrics (AUC, Brier, nDCG, Spearman), losses (LambdaRank,
hurdle, ordinal all-threshold), evaluation splits (bacteria-axis, phage-axis), framework pieces
(strata, MLC), etc. Entries are alphabetical and keep definitions short, concrete, and grounded in
how the term is used in *this* project — not textbook generality.

Each entry should link back to the knowledge unit or notebook entry where the term is load-bearing
(e.g. `See ordinal-regression-not-better`) so readers can find where we learned what we know.

Add a new entry when a term starts showing up in notebook entries, PR discussions, or PR review
comments and can't be understood from context alone. Update an existing entry when a later ticket
sharpens what the term means for us (e.g. a new result changes the "rough interpretation scale").

Do not duplicate knowledge units — the glossary explains terminology; the knowledge model records
findings.

## LITERATURE.md separation of concerns

`LITERATURE.md` documents **what external papers and resources offer**: findings, methods, data, tools, caveats about
the source itself. It must not contain:

- References to our internal tracks, experiments, or task codes (TI09, AX03, etc.)
- Our project's plans for integrating a resource ("directed cross-terms", "recommended role: ...")
- Our past mistakes or detection bugs
- Gap-to-solution maps, retired gaps, or other project planning
- Speculation about how features would work in our pipeline

Those belong in the plan (`plan.yml`), knowledge model (`knowledge.yml`), lab notebooks, or gists. LITERATURE.md is a
reference library, not a project roadmap. A reader should be able to understand each entry without knowing anything
about our internal codebase.

## Other contents

- `PLAN.md` is auto-generated from `lyzortx/orchestration/plan.yml` by `lyzortx/orchestration/render_plan.py`. Do not
  edit it by hand.
- `KNOWLEDGE.md` (at `lyzortx/KNOWLEDGE.md`) is auto-generated from `lyzortx/orchestration/knowledge.yml` by
  `lyzortx/orchestration/render_knowledge.py`. Do not edit it by hand.
- `LITERATURE.md` is a curated reading list maintained manually.
- `GLOSSARY.md` is a curated glossary of project-specific terminology (see section above).
- `external_data/` contains the source registry and external dataset metadata.
- `ad_hoc_analysis_code/` contains one-off analysis scripts referenced by lab notebook entries.
- `TIER_BENCHMARK_DENOMINATOR_POLICY.md` documents denominator rules for benchmark reporting.
