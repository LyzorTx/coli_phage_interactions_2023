# First Principles for Every Task

Before thinking, planning, or doing anything, follow these four steps in order:

1. **Question the requirement.** Should we be doing this at all? Delete the requirement if possible.
2. **Simplify.** Strip to bare bones. Remove every unnecessary constraint, layer, and edge case.
3. **Accelerate.** Make the implementation thinner, more elegant, and faster.
4. **Automate.** Only after steps 1-3. Automating a bad requirement locks in waste.

# Directory-Specific Policies

Detailed coding, testing, scientific review, CI, and orchestration policies live in subdirectory `AGENTS.md` files
(`lyzortx/`, `.github/workflows/`, `lyzortx/orchestration/`, etc.) and load automatically when working in those areas.

# Code Placement Policy

- All new content (code, data, and other resources) goes under `lyzortx/`. Only what must be outside goes outside
  (e.g., root `AGENTS.md`, `.github/workflows/`, environment manifests).
- This is a fork; keeping new work under `lyzortx/` separates it from upstream sources.
- Modify code outside `lyzortx/` only when fixing a bug in the original repository code.

# Environment Policy

- **Local development:** The `phage_env` micromamba env is auto-activated via direnv (`.envrc`). No `micromamba run`
  needed locally.
- **Codex shells:** If `python`/`pytest` are not from `phage_env`, source micromamba explicitly first, then
  `micromamba activate phage_env`.
- **GitHub Actions:** Bootstrap with `micromamba env create -f environment.yml -n phage_env`, then use
  `micromamba run -n phage_env ...` (always `-n`, never `-p <path>`).
- **Prebaked CI images:** If a tool is missing, add the dependency to the appropriate `environment*.yml` manifest.
  Never weaken validation to work around a missing tool.
- `micromamba run` does not source activation scripts — tools needing them (e.g., openjdk) require
  `micromamba activate` instead.
- Detect CI: check `GITHUB_ACTIONS=true`. Git identity is pre-configured in CI; do not set it.
- **Long-running local jobs:** Run `caffeinate -dims -w <pid> &` to prevent macOS idle sleep.
- **Generated outputs do not exist in CI.** Tasks must regenerate them or fail loudly — never silently produce empty
  results.

# Dependency Policies

- **Never run `pip install` or `micromamba install` directly.** Add dependencies to `requirements.txt` (pip) or
  `environment*.yml` (micromamba), then install from the manifest.
- Keep `requirements.txt` alphabetically sorted. Pin every dependency to an exact version (e.g., `ruff==0.11.6`).
- Prefer library-based approaches when they improve quality or maintainability — do not degrade implementations to
  avoid adding dependencies.

# Path Style in Commands

- Use **relative paths** in all shell commands. Never use absolute paths.
- Why: `.claude/settings.json` permissions use glob patterns that break with absolute paths.

# Git and PR Policies

- Committing directly to main is allowed. Use feature branches for orchestrator tasks or changes needing review.
- **On-ticket orchestrator work: push and open PRs without asking.** Implementing a dispatched orchestrator ticket is
  pre-authorised to push the feature branch, open the PR, run a fresh-context reviewer subagent, and merge if the
  subagent approves. This is the standard batch-ticket flow (one PR per ticket, sequential merges).
- **Off-ticket work on main: ask before pushing.** Opportunistic cleanup, refactors, or fixes that are not part of a
  dispatched orchestrator ticket require explicit confirmation before `git push`, because a premature push can trigger
  CI or auto-merge before the work is ready (as happened with PR #366).
- Always rebase on main before starting work and before every push:
  `git fetch origin main && git rebase origin/main` (or `gt sync` for Graphite stacks).
- Never use `git add -f`, `git add .`, or `git add -A`. Stage files by explicit path.
- One kind of change per commit. Do not mix unrelated changes.
- The user may have multiple sessions open. Never claim changes are uncommitted without running `git status` first.
- After every push to a PR branch, update the PR title and body to reflect all commits.
- Any PR addressing a tracked issue must include `Closes #<issue_number>` in the description.
- Create worktrees under `.claude/worktrees/<descriptive-name>` (the canonical location). Never place worktrees in
  sibling directories like `../foo` or `../<repo-name>-<branch>` — they end up next to the repo and are easy to lose
  track of. If you find an existing worktree outside this path, move it with `git worktree move` before continuing.
- Do not remove a worktree until its PR is merged or the user explicitly abandons it. After creating a worktree and
  pushing a PR, tell the user where it is.
- The Edit tool succeeds silently when `old_string` already equals `new_string`. Your edit appearing to succeed does
  **not** prove the change was uncommitted — another session may have committed it already.

# Plan-Driven Execution

- The main project driver is `lyzortx/orchestration/plan.yml`, rendered to `lyzortx/research_notes/PLAN.md`.
- Follow the plan for task sequencing; the orchestrator updates checklist states automatically.
- When scope is ambiguous, prefer alignment with the plan unless the user overrides.
- Use the `/replan` skill when completed work is questioned or assumptions are invalidated.

# Self-Review Before Requesting Human Review

- After pushing a PR, spawn a **reviewer subagent** to self-review the diff. Do NOT pass your conversation context to
  the reviewer — it must start fresh from the diff and acceptance criteria only (prevents confirmation bias).
- The reviewer assesses in priority order: (1) acceptance criteria met, (2) scientific/biological/logical correctness,
  (3) code correctness, (4) code cleanliness. See `lyzortx/orchestration/AGENTS.md` for full protocol.
- Fix any findings and re-run the reviewer until clean. When the reviewer approves, merge the PR and continue to the
  next task in the track.

# PR Creation for Orchestrator Tasks

- Create PRs with `gh pr create`. Title pattern: `[ORCH][TASK_ID] Brief description`.
- PR body MUST include `Closes #<issue_number>`. Add `--label orchestrator-task`.
- Always use a HEREDOC for the body. Use the `/gh create-pr` command for the canonical template.

# Issue Closure Policy

- Closing a GitHub issue signals completion to the orchestrator (`done`). Never close an orchestrator-task issue unless
  the task is genuinely finished and its acceptance criteria are met.
- To close **without** marking done: `gh issue close --reason "not planned"`.
- Never reopen a closed orchestrator issue — reopening triggers the implementation workflow. To invalidate, change close
  reason via API:
  `gh api repos/OWNER/REPO/issues/NUMBER -X PATCH -f state=closed -f state_reason=not_planned`.

# Graphite Stacked PRs

- Use plain `git` for single-branch workflows. Only use `gt` for stacked PRs.
- This repo is not synced with Graphite's remote. Use `git push` instead of `gt submit`.
- Use `gt create --no-interactive` for new branches in a stack.
- Prefer stacked PRs when the diff exceeds ~300 lines or has clear layered stages. Use the `/graphite` skill to create
  stacks.
- Each PR in a stack must be atomic and pass CI independently.

# Review Focus Areas

Reviews must prioritize in this order — reject on higher-priority issues before even looking at lower ones:

1. **Does the ticket actually work?** — Does the PR implement what the ticket asks for? Are acceptance criteria met?
   Does it produce real results (not scaffolding, not zero-row outputs)?
2. **Scientific/biological/logical correctness** — biological plausibility, statistical rigor, honest interpretation.
   Wrong biology or flawed statistics are worse than ugly code.
3. **Correctness** — bugs, logic errors, off-by-one, wrong variable usage.
4. **Test coverage** — are new/changed functions tested? Critical edge cases covered?
5. **Security** — no secrets committed, no injection risks.
6. **AGENTS.md compliance** — verify PR follows policies (code placement, dependency pinning, git staging, etc.).
7. **Clarity** — naming, structure, readability.
8. **Coding principles** — no magic numbers, constants defined and reused, progress logging with timestamps.

Do NOT nitpick style — ruff handles formatting. Focus on substantive issues only. Do not invent problems.

- Approve the PR based on code quality alone. Do not wait for CI checks — auto-merge already gates on CI.
- Every time an agent creates a PR or pushes an update, it must self-review against these guidelines before considering
  the PR ready.
- Before raising an issue, check existing review threads — do not re-raise concerns already addressed.
- Do not add "CI passes" as a test-plan item in PR descriptions — it is redundant noise.
- When addressing review feedback, push back on comments that are wrong, overcomplicated, or low-value.

# Requirement Challenge Policy

- For non-trivial requests, first question the requirement before implementing.
- Push back clearly if the request is unreasonable or overcomplicated. Suggest deletion or simplification.
- Applies to plan decisions and PR review feedback equally.

# Agent Transparency

- When posting on GitHub using human credentials, identify as an agent with model name and version.

# Knowledge Persistence Policy

- Write learnings and rules to `AGENTS.md` files, **not** user memory. User memory is only for truly personal info.

# AGENTS.md and CLAUDE.md Pairing

- Every `AGENTS.md` must have an accompanying `CLAUDE.md` that imports it with `@AGENTS.md`.

# Timestamps

- Use the computer's local timezone for all timestamps, including lab notebook entries.
- Get the time from `date '+%Y-%m-%d %H:%M %Z'` (local timezone with zone label).
- Lab notebook entry format: `### YYYY-MM-DD HH:MM TZ: Title` (e.g., `### 2026-04-12 00:30 CEST: Title`).

# Documentation Style Rules

- Line-length: 120 characters for prose. See `.pre-commit-config.yaml` and `.pymarkdown.yaml` for tool settings.
- On a feature branch or open PR, notebook-style research documents being developed in that branch are mutable until
  merge. Edit them in place so the merged version reflects the final state; once merged to `main`, they become
  historical records.

# Agent Scratch Space

- Write temp files to `.scratch/`, not `/tmp`. It is gitignored.

# Dead Code Policy

- Delete unused code immediately. No external consumers — all callers are in this repo.

# Custom Skills

- Skills live in `.agents/skills/` (canonical path). `.claude/skills` is a symlink.
- When searching for skill files, use `.agents/skills/` — Glob does not follow symlinks.
- Use the `/skill-creator` skill to create or modify skills.

# Commit Shortcut

- On `commit staged`: commit only currently staged files with a generated message. Do not stage, unstage, or modify
  anything. If no files are staged, report that and stop.
