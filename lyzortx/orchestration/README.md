# PLAN Orchestrator

This directory contains an event-driven orchestrator that executes tasks from `plan.yml` using GitHub Issues for
sequencing and agent dispatch.

## Automation Lifecycle

The issue-to-merge lifecycle is automated across two GitHub Actions workflows and local Claude sessions:

```mermaid
stateDiagram-v2
    [*] --> task_ready : orchestrator.yml selects <br> ready task from plan.yml

    task_ready --> issue_open : orchestrator creates <br> GitHub issue with <br> orchestrator-task label

    issue_open --> local_implement : local Claude <br> picks up issue

    local_implement --> pr_created : local Claude implements task, <br> writes findings to lab_notebooks/track_&lt;track&gt;.md, <br> and opens PR with Closes #issue

    pr_created --> claude_review : claude-pr-review.yml <br> auto-reviews PR

    claude_review --> has_feedback : review has <br> inline comments
    claude_review --> no_issues : review approved

    has_feedback --> local_fix : local Claude reads <br> feedback and pushes fixes

    local_fix --> claude_review : push triggers <br> re-review

    no_issues --> auto_merge : auto-merge enabled <br> wait for CI

    auto_merge --> auto_merged : CI passes — <br> squash merge
    auto_merge --> local_merge : auto-merge fails — <br> local Claude merges

    auto_merged --> issue_closed : PR merge <br> auto-closes issue

    local_merge --> issue_closed : PR merge <br> auto-closes issue

    issue_closed --> orchestrator_tick : orchestrator.yml <br> triggers on issue.closed

    orchestrator_tick --> plan_updated : mark task done <br> in plan.yml <br> regenerate PLAN.md

    plan_updated --> task_ready : commit updates <br> and dispatch next task
```

### Actors

<!-- pyml disable-num-lines 8 md013-->
| Actor | Type | Role |
|---|---|---|
| `orchestrator.yml` | GitHub Actions workflow | Selects ready tasks, creates issues, marks tasks done, commits plan updates |
| Local Claude | Laptop-based Claude session | Implements tasks, pushes PRs, reads review feedback, pushes fixes, merges if auto-merge fails |
| `claude-pr-review.yml` | GitHub Actions workflow | Auto-reviews every PR on open/push via `claude-code-action`. Claude (Opus 4.6) submits formal `APPROVE`/`COMMENT` reviews and manages thread resolution. Posts as `claude[bot]`. On approval, enables auto-merge. |
| `autoresearch-runpod.yml` | GitHub Actions workflow | Manual AUTORESEARCH search runner: stages a frozen cache bundle, provisions a locked RunPod pod behind an environment approval gate, executes one bounded `train.py` command, uploads candidate artifacts, and tears the pod down |

#### Disabled workflows

| Workflow | Reason |
|---|---|
| `codex-implement.yml` | Replaced by local Claude sessions for task implementation |
| `codex-pr-lifecycle.yml` | Replaced by local Claude sessions for review feedback fixes |

## Architecture

- **Source of truth:** `lyzortx/orchestration/plan.yml` — all tracks, tasks, dependencies, status, and acceptance
  criteria.
- **Rendered view:** `lyzortx/research_notes/PLAN.md` — auto-generated from `plan.yml` by `render_plan.py`. CI verifies
  it stays in sync.
- **Issue state:** GitHub issues labeled `orchestrator-task` are the authoritative progression signal. When an issue
  closes, the orchestrator marks the task `done` in `plan.yml` and regenerates `PLAN.md`.
- **Runtime state:** `lyzortx/generated_outputs/orchestration/runtime_state.json` — ephemeral per CI run, uploaded as
  artifact.

## Components

- `plan.yml` — task definitions (source of truth).
- `plan_parser.py` — pure functions: `load_plan`, `is_task_ready`, `resolve_task_dependencies`,
  `select_ready_tasks`, `mark_task_done`. Parses `model`, optional `depends_on_tasks`, and optional
  `ci_image_profile` fields from task entries.
- `ci_image_profiles.py` — shared enum/mapping helpers for `ci_image_profile`, `ci-image:*` labels, and GHCR image
  refs used by the orchestrator and workflows.
- `parse_model_directive.py` — extracts model ID from `<!-- model: ... -->` HTML comments in issue bodies. Used by CI
  workflows and available as a CLI: `echo "$BODY" | python -m lyzortx.orchestration.parse_model_directive`.
- `render_plan.py` — generates `PLAN.md` from `plan.yml` with Mermaid DAG and track checklists.
- `orchestrator.py` — CLI runner that dispatches tasks as GitHub issues.
- `review_threads.py` — fetches unresolved PR review threads via GitHub GraphQL, paginates across thread pages, filters
  to unresolved non-outdated threads, and formats them into a Codex feedback prompt. This remains the merge-gate helper
  used by `claude-pr-review.yml`.
- `pr_feedback.py` — fetches top-level PR conversation comments, inline review comments, and non-empty review bodies,
  then formats them into a Codex feedback prompt so lifecycle runs read every visible PR feedback surface.
- `verify_review_replies.py` — checks that PR review comments have been addressed with replies.
- `ci_token_usage.py` — CLI for token/cost analysis across all LLM-invoking workflows (Codex + Claude).
- `.github/workflows/orchestrator.yml` — CI trigger: task dispatch and plan updates.
- `.github/workflows/codex-implement.yml` — disabled; retained for reference.
- `.github/workflows/codex-pr-lifecycle.yml` — disabled; retained for reference.
- `.github/workflows/autoresearch-runpod.yml` — manual GPU experiment runner for AUTORESEARCH. It either builds a
  fresh host-side cache bundle or reuses one from a prior workflow run, then provisions a fixed RunPod pod, copies only
  the AUTORESEARCH runtime bundle plus frozen cache artifacts, runs one bounded `train.py` command, serializes
  workflow-dispatch metadata via `jq` instead of shell-built JSON, uploads candidate outputs, and deletes the pod.
- `.github/workflows/ci-duplicate-check.yml` — informational CI check: runs pylint `symilar` to detect duplicate code
  in `lyzortx/`. Does not block PRs (`continue-on-error: true`).

## Task Readiness

A task is ready when:

1. If the task does not declare `depends_on_tasks`, all prior tasks in the same track are `done` (sequential by
   default within track).
2. If the task does declare `depends_on_tasks`, only those explicit task IDs block it within the track. This is how
   the plan expresses intra-track parallelism such as "TL15/TL16/TL17 can start together, TL18 waits on all three."
3. All tasks in all prerequisite tracks (from `depends_on`) are `done`.

Task IDs are derived from track letter + ordinal (e.g., `TB03`, `TF01`). Gates use `GNG` prefix.

## Task Authoring Guidance

Plan authors should size tasks by boundary risk, not just by how small the diff sounds.

- For fragile tasks, write low-freedom acceptance criteria. State the exact contract that changed and the exact failure
  modes to avoid.
- When a task introduces a fallback, acceptance criteria should require both:
  - a positive test for the intended narrow use
  - a negative test proving strict failure still happens outside that use
- When a task consumes generated artifacts, acceptance criteria should say whether stale default artifacts must be
  regenerated or rejected.
- For AUTORESEARCH-style tasks, treat raw inputs plus frozen featurizer code as the source of truth. Checked-in feature
  CSVs may be optional warm caches only; acceptance criteria should require rebuildability from raw data and should
  exclude panel-only metadata or proxies that cannot run on unseen FASTAs.
- For AUTORESEARCH plan design, split cache-building tasks by runtime-risk boundary when the stages use different
  toolchains or cost profiles. In this repo that means separate tickets for host defense, host surface, host typing,
  and phage projection instead of one broad "implement prepare.py" task.
- For AUTORESEARCH critical-path design, make the adsorption-first minimum cache sufficient for the first baseline when
  that is the most credible early substrate. Slower optional blocks such as host defense can join later as additive
  ablations instead of gating the first runnable search.
## CLI Usage

```bash
# Show status with ready tasks
python -m lyzortx.orchestration.orchestrator --command status --plan-path lyzortx/orchestration/plan.yml

# Dispatch one ready task (creates GitHub issue when GITHUB_TOKEN is set)
python -m lyzortx.orchestration.orchestrator --command run_once --plan-path lyzortx/orchestration/plan.yml

# Pause/resume
python -m lyzortx.orchestration.orchestrator --command pause --note "maintenance"
python -m lyzortx.orchestration.orchestrator --command resume

# Regenerate PLAN.md from plan.yml
python -m lyzortx.orchestration.render_plan
```

## GitHub Actions Trigger Model

### orchestrator.yml

- `workflow_dispatch`: manual commands (`run_once`, `status`, `pause`, `resume`).
- `repository_dispatch`: API/CLI command trigger.
- `issues.closed`: when an `orchestrator-task` issue closes, marks the task done and dispatches the next ready task.

A concurrency group (`orchestrator`) queues runs instead of running in parallel, preventing duplicate issue creation
when multiple trigger events fire simultaneously.

On each tick the workflow commits `plan.yml` and `PLAN.md` changes back to the repo.

Default `max_active_tasks` is `1` (CLI) or `50` (CI workflow). The `orchestrator-task` label is created automatically on
first dispatch. Dispatched issues also receive a `model-{id}` label (e.g., `model-gpt-5.4-mini`) for at-a-glance model
visibility plus a mirrored `ci-image:{profile}` label so workflows can route the task to the matching prebaked
container image. Missing profile labels are treated as configuration errors, not as permission to fall back.

### codex-implement.yml (disabled)

Disabled — task implementation is now handled by local Claude sessions. The workflow file is retained for reference.

### claude-pr-review.yml

- `pull_request: [opened, synchronize]`: auto-reviews every PR on open or push.
- `issue_comment: [created]` / `pull_request_review_comment: [created]`: interactive `@claude` mentions.

Claude reads `AGENTS.md` review guidelines, submits formal `APPROVE` or `COMMENT` reviews via MCP GitHub tools, and is
the sole judge of thread resolution (can resolve/unresolve threads via GraphQL mutations). Requires the
`ANTHROPIC_API_KEY` repository secret. The workflow explicitly allows the repo's `czarphage` GitHub App bot to trigger
re-reviews after local Claude pushes, which would otherwise be blocked by `claude-code-action`'s default "no bots"
policy. After reviewing, it auto-merges only when Claude's latest review is `APPROVED` and the shared
`lyzortx.orchestration.review_threads` helper reports zero unresolved review threads. If Claude leaves a `COMMENTED`
review or any unresolved review threads remain, local Claude reads the feedback and pushes fixes.

### codex-pr-lifecycle.yml (disabled)

Disabled — review feedback is now addressed by local Claude sessions. The workflow file is retained for reference.

### autoresearch-runpod.yml

- `workflow_dispatch`: manual AUTORESEARCH GPU search only.

This workflow is deliberately outside the generic Codex issue/PR loop. It has two phases:

1. **Stage the frozen host-side bundle.** Either:
   - build a fresh AUTORESEARCH cache on the GitHub-hosted side with `prepare.py`, then package the minimal runtime
     bundle; or
   - download a previously staged `autoresearch-runpod-bundle` artifact from an earlier workflow run.
2. **Run one bounded pod-side experiment.** After environment approval, provision one locked single-GPU RunPod pod,
   copy only the staged bundle, create `phage_env` inside the pod, run one bounded `train.py` command, pull candidate
   outputs back as a workflow artifact, and delete the pod.

The workflow exists to keep expensive, infrequent host-side cache building separate from many short `train.py` search
runs. Repeated experiments should normally point at an existing staged bundle instead of rebuilding AR03-AR06 outputs.

## AUTORESEARCH RunPod Contract

This is a human-approved lock, not an auto-selected cloud default.

- **Required GitHub environment:** `runpod-autoresearch`
- **Required environment secret:** `RUNPOD_API_KEY`
- **Approval gate:** configure required reviewers on the `runpod-autoresearch` GitHub environment. The workflow's
  RunPod job will not start until that environment is approved.
- **Locked pod spec:** community-cloud single-GPU `NVIDIA A40`, `48 GB` VRAM, `1` GPU, `50 GB` container disk,
  `20 GB` volume, public IP enabled, image `runpod/pytorch:2.1.0-py3.10-cuda11.8.0-devel-ubuntu22.04`.
- **Locked hourly price point:** `$0.35/hr` on RunPod's official GPU pricing page at the time AR08 was implemented.
- **Why this pod:** `train.py` is now a thin cache consumer, so it does not need the cheapest possible auto-selected
  GPU. The A40 lock buys more memory headroom than the nearby 24 GB options for essentially the same hourly cost, while
  staying within a single-GPU community-cloud budget.
- **Secret boundary:** `RUNPOD_API_KEY` is referenced only on the provisioning, status-polling, and teardown steps in
  `autoresearch-runpod.yml`.
- **Host-to-RunPod handoff:** the host side packages exactly one bundle artifact containing:
  - `lyzortx/autoresearch/{prepare.py,train.py,README.md,program.md}`
  - minimal runtime support files (`environment.yml`, `requirements.txt`, `pyproject.toml`,
    `lyzortx/log_config.py`, `lyzortx/pipeline/autoresearch/runtime_contract.py`)
  - the frozen `lyzortx/generated_outputs/autoresearch/search_cache_v1/` tree
- **Pod-side responsibilities:** unpack the bundle, create `phage_env`, run the bounded `train.py` command, collect
  `ar07_baseline_summary.json`, `ar07_inner_val_predictions.csv`, the exact `train.py`, and RunPod metadata.
- **Out of bounds for the pod:** `prepare.py`, Picard assembly download, DefenseFinder execution, host typing calls,
  host-surface derivation, and phage projection rebuilding. Those stay on the host side or in prior staged bundles.

## CI Image Profiles

The `ci-image:*` label system is retained for `claude-pr-review.yml` routing and for reference if CI-based
implementation is re-enabled. Orchestrator issues still carry `ci-image:*` labels.

- `base` — `phage_env` only
- `host-typing` — `phage_env` plus `phylogroup_caller`, `serotype_caller`, and `sequence_type_caller`
- `full-bio` — `host-typing` plus `phage_annotation_tools`

## Agent Instructions in Dispatched Issues

Each dispatched issue includes:

- Task description and acceptance criteria (from `plan.yml`).
- Model directive as an HTML comment: `<!-- model: gpt-5.4-mini -->`. The model is set per-task in `plan.yml` and
  emitted by `orchestrator.py` when creating the issue. Both `model` and `acceptance_criteria` are required for all
  pending tasks — the orchestrator raises `ValueError` if either is missing. (The model directive is a historical
  artifact from when Codex CI consumed it; it is retained for traceability.)
- CI image profile directive as an HTML comment and mirrored GitHub label. The profile is set per-task in `plan.yml`
  via `ci_image_profile` and mirrored into `ci-image:{profile}` labels for issue/PR routing. Pending tasks must declare
  this explicitly; missing labels fail the workflow rather than silently falling back.
- Instruction to write findings to `lyzortx/research_notes/lab_notebooks/track_<track>.md`.
- PR creation instructions using `gh pr create` with `Closes #<issue>`.
