# Orchestration Directory

- When modifying files in this directory, update `lyzortx/orchestration/README.md` to reflect the changes.
- The README contains a Mermaid state diagram of the full automation lifecycle and descriptions of all components. Keep
  both in sync with the actual workflow logic.

# Implementation Workflow

- `codex-implement.yml`, `codex-pr-lifecycle.yml`, and `claude-pr-review.yml` are all disabled in GitHub Actions. Task
  implementation, review, and feedback fixes happen locally (laptop-based Claude sessions), not in CI.
- The orchestrator (`orchestrator.yml`) still runs in CI: it ticks on issue close events and `workflow_dispatch`, marks
  tasks done, and creates new issues.
- Full lifecycle:
  1. Orchestrator tick creates a GitHub issue for the next pending task.
  2. Local Claude implements the task and pushes a PR (`Closes #N`).
  3. Local Claude spawns a **reviewer subagent** (see below) to self-review the PR.
  4. If the self-review has findings, local Claude fixes them and re-runs the reviewer.
  5. When the reviewer approves, local Claude merges the PR.
  6. The merged PR closes the linked issue, which triggers an orchestrator tick to dispatch the next task.
  7. Local Claude continues to the next task in the track without waiting for human input.

# Self-Review via Subagent

- After pushing a PR, the implementing agent must spawn a **separate reviewer subagent** to review the diff.
- The reviewer subagent must NOT receive the implementing agent's conversation context — it starts fresh from the diff
  and acceptance criteria only. This prevents confirmation bias.
- The reviewer subagent assesses the PR in this priority order (reject on higher-priority issues before looking at
  lower ones):
  1. **Acceptance criteria met?** Does the PR implement what the ticket asks for? Are all criteria satisfied with real
     results (not scaffolding, not zero-row outputs)?
  2. **Scientific/biological/logical correctness.** Biological plausibility, statistical rigor, honest interpretation.
     Wrong biology or flawed statistics are worse than ugly code. **If the PR writes prediction or metric artifacts
     under `lyzortx/generated_outputs/` (CSV, JSON) or cites AUC / Brier / bootstrap CI numbers, the reviewer subagent
     MUST invoke the `review-ml-pr` skill as part of this step.** Reading the diff and acceptance criteria is not
     sufficient verification — the skill recomputes headline numbers from the artifacts, checks fold disjointness
     empirically, decomposes aggregates into strata, and builds reliability diagrams. A read-only review of an ML
     results PR is not an acceptable self-review; the implementing agent's spawn prompt must explicitly instruct the
     reviewer to run `review-ml-pr` in that case.
  3. **Code correctness.** Bugs, logic errors, off-by-one, wrong variable usage.
  4. **Code cleanliness.** Naming, structure, readability. Do NOT nitpick style — ruff handles formatting.

# Knowledge Model Rendering

- `knowledge.yml` is the source of truth for consolidated project knowledge. `render_knowledge.py` renders it to
  `lyzortx/KNOWLEDGE.md`.
- Every field in `knowledge.yml` must appear in the rendered output. If a field exists in the YAML, the renderer must
  surface it — statement, sources, status, confidence, context, and relates_to.
- When adding new fields to `knowledge.yml`, update `render_knowledge.py` to render them.
- Run `python -m lyzortx.orchestration.render_knowledge` after modifying `knowledge.yml` to regenerate `KNOWLEDGE.md`.
- The validator (`knowledge_parser.validate_knowledge()`) must pass before rendering. Invalid YAML is rejected.

# Plan Rendering Completeness

- Every field in `plan.yml` must appear in the rendered `PLAN.md`. If a field exists in the YAML, the renderer must
  surface it — acceptance criteria, model, implemented_in, baseline, description, and any future fields.
- When adding new fields to `plan.yml`, update `render_plan.py` to render them. Do not add YAML fields that are silently
  dropped during rendering.

# Task ID Numbering

- Task IDs use a track prefix followed by an integer with no letter suffix: `CH01`, `CH02`, `CH03`, never `CH02a` /
  `CH02b`. Letter suffixes hide the true ticket count, confuse grep, and imply a substructure the orchestrator doesn't
  model.
- If a task feels too big for one ID, split it into separate numbered tasks (`CH02`, `CH03`, `CH04`) and renumber
  subsequent tasks. Do not paper over the split with letters.
- Renumbering pending tasks is cheap. Renumbering done tasks is not allowed — done IDs are referenced across
  notebooks, knowledge units, and PR history (see Done Task Immutability below). If a task sequence needs to change
  after some entries are done, keep the done IDs stable and let the pending numbering flow around them.

# Pending Task Requirements

- Every pending task in `plan.yml` must have both a `model` field and non-empty `acceptance_criteria`.
- Done tasks may omit either — they are historical records, not dispatchable work.
- The orchestrator validates both fields at load time and raises `ValueError` if either is missing on a pending task.

# Acceptance Criteria for Artifact-Boundary Tasks

- If a task touches generated artifacts or evaluation semantics, the acceptance criteria must constrain boundary
  behavior explicitly. Include, when relevant:
  - a stale-artifact check: existing default outputs must be regenerated or rejected if their schema/provenance is old
  - a provenance-path check: validation must follow the actual CLI-provided artifact path, not a hardcoded default
  - a narrow-fallback check: permissive behavior must be limited to the intended context
- When adding a permissive fallback, require two tests:
  - one proving the fallback works in the intended narrow context
  - one proving strict failure still occurs outside that context
- Never accept criteria that say only "rerun downstream task on clean outputs" when the upstream change alters what rows
  or files exist. The task text must restate the changed contract.

# Task Status Management

- Never manually set a task's `status` to `done` in `plan.yml`. The orchestrator automatically marks tasks as done when
  their corresponding GitHub issue is closed as completed (via PR merge with `Closes #N`).
- Agents should only add new tasks or modify pending task fields (acceptance_criteria, model, title). Status transitions
  are the orchestrator's responsibility.

# Done Task Immutability

- Never modify the title, acceptance_criteria, or any other field on a done task. Done tasks are historical records —
  they document what was originally required, not the current state of the code.
- If a later task changes the codebase in ways that make a done task's criteria look stale (e.g., deleting features that
  a done task created), that is expected. The done task records what was true when it was completed; the later task
  records what changed. Git history tells the full story.

# Push Before Closing Issues

- Always push code changes to main before closing GitHub issues. Issue close events trigger orchestrator runs
  immediately. If the new code hasn't been pushed yet, the old code runs and can recreate the issue you just closed.
- This applies to both "completed" and "not planned" closures — any close event triggers a tick.

# Orchestration Robustness vs Fail-Fast

- The root AGENTS.md fail-fast rule applies to pipeline code. Orchestration code is different: the orchestrator should
  be robust against transient failures (GitHub API errors, network timeouts, race conditions) via retries.
- But robustness means retrying on transient errors, not tolerating missing data or silently skipping work. If the
  orchestrator cannot dispatch a task because its inputs are missing, that is an error to surface, not a condition to
  swallow.

# PR Lifecycle Feedback Contract

- `codex-pr-lifecycle.yml` and `claude-pr-review.yml` are both disabled. Review is handled by the local implementing
  agent's self-review subagent (see "Self-Review via Subagent" above), followed by human review.
- The self-review subagent replaces the CI-based Claude Auto Review. It runs locally, costs no CI minutes, and provides
  faster feedback loops.
