### 2026-04-08: RunPod REST v1 API — schema keeps changing, GPU availability is unreliable

#### Executive summary

Continued RunPod API debugging after the initial workflow bash replacement. The REST v1 `/pods` endpoint changed its
schema validation between our GitHub Actions runs and local testing, requiring a second round of payload fixes. The
locked GPU type (NVIDIA A40) was unavailable on community cloud, forcing a pivot to local CPU execution.

#### Findings — API schema instability

The 2026-04-07 entry recorded that `gpuTypeId` (singular string) was the correct field and `gpuTypeIds` (array) caused
errors. When we ran the same payload locally on 2026-04-08, the API reversed its position:

- `gpuTypeId` (singular) now causes HTTP 400: `"key provided in request body which is not in input schema: 'gpuTypeId'
  (Did you mean 'gpuTypeIds'?)"`.
- `gpuTypeIds` (array of strings) is now the accepted form: `"gpuTypeIds": ["NVIDIA A40"]`.
- `dockerArgs` now causes HTTP 400: `"key provided in request body which is not in input schema: 'dockerArgs'"`.
  Previously it was required as an empty string.

This means the RunPod REST v1 schema is not stable across time. The unit tests in
`lyzortx/tests/test_runpod_orchestrator.py` now assert the latest known-good schema, but this may break again.

#### Findings — GPU availability

Queried available GPU types via the RunPod GraphQL API (`gpuTypes { id communityCloud secureCloud }`). Key finding:
NVIDIA A40 has `communityCloud=False` — it is only available on secure cloud, which may have different pricing and
provisioning. Available community alternatives include RTX 3090 (24GB, ~$0.20/hr), RTX 4090 (24GB, ~$0.35/hr), and
RTX A6000 (48GB). For our LightGBM workload (64 trees, 15 leaves), even the smallest GPU is massive overkill — the
entire training runs in 0.4 seconds on CPU.

#### Findings — local execution is viable

Since `train.py` supports `--device-type cpu` and the current LightGBM baseline is lightweight, the entire
AUTORESEARCH pipeline can run locally without RunPod:

- `prepare.py --skip-comparator-lock`: ~25 min first build (369 hosts), ~5 sec on resume with existing slots.
- `train.py --device-type cpu`: 0.4 seconds.
- `replicate.py replicate`: ~20 seconds (3 seeds × candidate + comparator + 1000 bootstrap samples).

The RunPod orchestrator and workflow remain useful for future heavier workloads but are not needed for the current
baseline.

#### Resolution

Updated `runpod_orchestrator.py` to use `gpuTypeIds` (array) and removed `dockerArgs`. The first honest baseline was
run entirely locally on CPU. RunPod provisioning is deferred until the workload justifies GPU rental.

### 2026-04-07: RunPod REST v1 API field format learnings (initial round)

#### Executive summary

First live run of the AUTORESEARCH RunPod workflow required several iterations to get the pod creation payload right.
The RunPod REST v1 `/pods` endpoint is picky about field names and types in ways that contradict their own docs. See
the 2026-04-08 entry above for subsequent corrections — the schema changed between this session and the next.

#### Findings

- `ports` must be an **array** (`["22/tcp"]`), not a string. The REST schema validation rejects `"22/tcp"`.
- `computeType` and `allowedCudaVersions` are not valid REST v1 fields — they cause 500 errors.
- The `env` field is a plain object (`{"SSH_PUBLIC_KEY": "..."}`) — not an array of key-value pairs.
- Using `curl --fail` swallows the response body on HTTP errors, making debugging impossible. Always capture the body
  separately with `--output` and `--write-out "%{http_code}"`.

#### Resolution

Replaced inline workflow bash with `lyzortx/pipeline/autoresearch/runpod_orchestrator.py` — a testable Python module
with unit tests covering the payload schema. Future field-format bugs will be caught by
`lyzortx/tests/test_runpod_orchestrator.py` before reaching CI.

### 2026-03-31: Codex CI image routing now follows task-level image profiles

#### Executive summary

The first cut of the prebaked Codex CI image put every environment into one default container image. That worked
mechanically but over-baked the default path, especially once `serotype_caller` pulled ECTyper's large species-ID
database during image build. The workflows now route tasks through a small set of named image profiles driven by
`ci_image_profile` in `lyzortx/orchestration/plan.yml`, mirrored as `ci-image:*` labels on orchestrator issues and
PRs.

#### Design decision

Keep `plan.yml` as the source of truth and GitHub labels as a derived routing surface. Tasks can now declare one of
three image profiles:

- `base`: `phage_env` only
- `host-typing`: `phage_env` plus `phylogroup_caller`, `serotype_caller`, and `sequence_type_caller`
- `full-bio`: `host-typing` plus `phage_annotation_tools`

The orchestrator mirrors that choice into a `ci-image:{profile}` label when it opens the task issue. The implement
workflow copies the same label onto the PR after Codex creates it. Both Codex workflows now fail loudly if the required
label is missing or duplicated instead of silently falling back to a default image.

#### Why this is better

This keeps `phage_env` lightweight and solver-friendly while still allowing specialized tasks to opt into heavier tool
stacks. It also avoids baking the heaviest bioinformatics and host-typing assets into every single Codex job by
default. The GHCR publish workflow now builds one image per profile instead of one universal image, and the Codex
workflows only refresh the envs that exist in the selected profile.

### 2026-03-31: Unit-test workflow now emits slow-test timings and runs pytest in parallel

#### Executive summary

The `Unit Tests` workflow was running `pytest -q`, which kept CI logs too quiet to diagnose where the 5-minute wall
time was going. Step timings showed the bottleneck was the pytest step itself, not checkout or dependency install. The
workflow now drops `-q`, adds `--durations=25` to surface the slowest tests directly in GitHub Actions logs, and runs
the suite with `pytest-xdist` using `-n auto --dist loadfile`.

#### Why this change was needed

For recent PRs, the `Run pytest with coverage` step was taking about 4 minutes on its own, while install and
pre-commit together were under a minute. With `-q`, the log only showed the final `[100%]` line, which made it
impossible to tell whether the time was coming from one pathological test, a small set of heavyweight integration
tests, or broad suite-wide slowdown on the runner.

#### Design decision

Prefer observability first, then parallelize with a conservative sharding mode. `--durations=25` gives actionable
timing data from CI, and removing `-q` restores visible progress output. For parallel execution, the workflow uses
`-n auto --dist loadfile` rather than per-test sharding so tests from the same file stay together, which reduces the
risk of fixture and file-system interference.

#### Local verification

After adding `pytest-xdist` to `requirements.txt`, ran the full suite locally with the exact workflow-style command:

`pytest -n auto --dist loadfile --durations=25 --cov=lyzortx --cov-report=xml:coverage.xml lyzortx/tests`

It passed cleanly (`382 passed`) and reduced local wall-clock time from the matching serial coverage run
(`pytest -q --cov=lyzortx --cov-report=xml:coverage.xml lyzortx/tests`, ~116s) to ~74s on this machine, so enabling
the same mode in CI is justified.

### 2026-03-31: Orchestrator now supports explicit task-level dependencies

#### Executive summary

The plan orchestrator used to assume every task within a track was strictly sequential. That was too rigid for the next
Track L replan, where `TL15`, `TL16`, and `TL17` should dispatch in parallel while `TL18` stays blocked on all three.
The orchestrator now accepts optional per-task `depends_on_tasks` entries in `plan.yml` and uses them for readiness
when present.

#### Design decision

Keep the old behavior as the default. If a task does not declare `depends_on_tasks`, earlier tasks in the same track
still block it. If a task does declare `depends_on_tasks`, those explicit task IDs replace the default same-track
serialization for that task only, while cross-track `depends_on` rules still apply. This keeps existing tracks stable
and makes parallelism opt-in rather than globally reinterpreting the plan.

#### Verification

- Added parser tests showing `TL15`, `TL16`, and `TL17` are ready together while `TL18` remains blocked until all three
  are done.
- Added an orchestrator `run_once` test proving the next tick would mark `TL15`/`TL16`/`TL17` in progress together and
  leave `TL18` pending.
- Re-rendered `PLAN.md` support so explicit task dependencies appear in the human-readable plan output.

### 2026-03-30: Claude review auto-merge was too permissive for commented approvals

#### Executive summary

PR #279 was closed too fast because `.github/workflows/claude-pr-review.yml` treated Claude's latest review state of
`APPROVED` as sufficient to enable auto-merge, even when review comments still needed addressing. The workflow now
requires both an `APPROVED` latest review and zero unresolved review threads before enabling auto-merge. If unresolved
review feedback remains, it dispatches the Codex lifecycle instead of merging.

#### Root cause

The previous gate looked only at the latest review state returned by `repos/{owner}/{repo}/pulls/{pr}/reviews`. That is
not a strong enough merge condition for this workflow because an approval and actionable inline comments can coexist on
the same PR. In practice, that let PR #279 merge before the commented feedback was addressed, which closed the PR while
review work was still outstanding.

#### Design decision

Interpret "zero comments" operationally as "zero unresolved review threads," not "no historical comments exist on the
PR." Historical comments remain visible on GitHub after they are addressed, so requiring literal zero comment objects
would make reviewed PRs permanently unmergeable. Unresolved-thread count matches the real requirement: do not merge
while any review feedback is still open.

#### Implementation

Updated `claude-pr-review.yml` to gate auto-merge on the unresolved-thread count returned by the shared
`lyzortx.orchestration.review_threads` helper. That keeps the Claude merge gate aligned with the Codex lifecycle, which
already uses the same helper to decide whether feedback remains. The workflow still accepts a clean Claude approval
path, but any unresolved review thread now forces the "request addressing feedback" path through the Codex PR
lifecycle workflow.

#### Follow-up after PR #283

PR #283 showed that author filtering was the wrong design entirely. The intended policy is that any unresolved review
thread blocks auto-merge, not only comments authored by Claude. The GraphQL author-login mismatch (`claude[bot]` in
REST, `claude` in GraphQL) was just the symptom that exposed this. The correct fix is to remove author filtering and
count every unresolved review thread on the PR. A second follow-up was needed after review pointed out that the
workflow's inline GraphQL query also disagreed with `review_threads.py` on outdated threads and silently capped itself
at 100 threads. Moving the count into the shared helper, and adding thread pagination there, keeps the automation
consistent.

### 2026-03-29: Fix conda environment solve failure — downgrade openjdk 25 to 24

#### Executive summary

The `codex-implement.yml` workflow (and `claude-pr-review.yml`, `codex-pr-lifecycle.yml`) failed on `conda env create`
because `openjdk=25.0.2` on linux-64 requires `harfbuzz >=12.3.2`, which pulls in `icu >=78`. Meanwhile,
`mmseqs2=18.8cc5c` depends on `aria2` → `libxml2 <2.14` → `icu <76`, creating an unsatisfiable constraint. Fixed by
pinning `openjdk=24.0.2`, whose linux-64 build requires only `harfbuzz >=11.4.5` (compatible with `icu 75`).

#### Root cause

Commit `1140205` ([ORCH][TL01]) added bioinformatics tools to `environment.yml` including `openjdk=25.0.2` and
`mmseqs2=18.8cc5c`. On macOS (osx-arm64), `openjdk` has minimal dependencies (`libzlib` only), so the environment
solved fine locally. On linux-64 (CI runners), `openjdk` has heavy GUI/font dependencies including `harfbuzz`, which
introduces transitive `icu` version constraints that conflict with `mmseqs2`'s `aria2` → `libxml2` chain.

#### Fix

Changed `openjdk=25.0.2` → `openjdk=24.0.2` in `environment.yml`. JDK 24 is the latest non-EA release before 25 and
its linux-64 conda-forge build requires `harfbuzz >=11.4.5` (not `>=12.3.2`), which is compatible with `icu 75.x`
shared by `libxml2`. JDK 25 (the bleeding-edge GA) added `harfbuzz >=12.3.2` which forced `icu >=78` — incompatible
with bioconda's `mmseqs2` package chain.

JDK 23 was also viable (`harfbuzz >=10.2`), and JDK 21 LTS (`harfbuzz >=8.2` for older builds) would have worked too,
but the latest 21.0.10 build also requires `harfbuzz >=12.3.2` — so the safest minimal downgrade is JDK 24.

Pharokka 1.9.1 uses Java only for minced (CRISPR detection) and tRNAscan, both of which work with JDK 24. No
pharokka-specific JDK 25 features are needed.

#### Bonus: strip conda from Claude PR review workflow

While investigating, found that `claude-pr-review.yml` (both auto-review and interactive jobs) was creating the full
`phage_env` conda environment — including pharokka, mmseqs2, openjdk — just to run Python orchestration helpers
(`review_threads`, `parse_model_directive`) and pre-commit. Replaced with `setup-python` + `pip install -r
requirements.txt`, which is faster and avoids the bioinformatics dependency chain entirely. The `codex-pr-lifecycle.yml`
workflow was kept on conda because it runs Codex which may need to execute pipeline code.

#### Verification

- Rebuilt `phage_env` from scratch locally with `openjdk=24.0.2` — solved and installed cleanly.
- Verified `java -version` (24.0.2), `mmseqs version` (18.8cc5c), `pharokka.py --version` (1.9.1) all functional.
- CI verification pending on next push.

### 2026-03-24: gh skill parse-logs command for CI log analysis

#### Executive summary

Added a `parse-logs` command to the `/gh` skill that documents how to fetch, filter, and extract timing from GitHub
Actions workflow logs via `gh run view --log`. This exists because the GitHub Actions web UI strips timestamps from
individual log lines — only showing sequential line numbers — making it impossible to diagnose timing of long-running
Codex implement and feedback runs from the browser. The raw logs accessed via CLI include full ISO 8601 timestamps.

#### Design decisions

**1. Skill command, not a script.** The `gh` CLI already provides all the primitives (`gh run view --log`, `--log-failed`,
`--job`). What was missing was documented recipes for combining them — filtering by PR + workflow + step, extracting
wall-clock duration, and finding stuck points. A skill command (agent documentation) is the right level of abstraction;
a wrapper script would add maintenance burden for no additional capability.

**2. All filters are combinable.** The command accepts any combination of: run ID, PR number, workflow name, job name,
step name, text pattern, status, and failed-only. Previous draft used "one of" which was too restrictive — real
debugging involves layering filters (e.g., "errors in the Codex step of PR #42's latest lifecycle run").

**3. Timing section with step-transition timeline.** Beyond simple first/last timestamp extraction, includes an awk
command that shows step transitions with timestamps — useful for seeing where time is spent across a run.

### 2026-03-24: Pre-push hook to enforce rebase on origin/main

#### Executive summary

Added a `check-rebase-on-main` pre-push hook via pre-commit that blocks `git push` when the branch is not rebased on
`origin/main`. The hook enforces two things: (1) `origin/main`'s tip is an ancestor of HEAD, and (2) no merge commits
exist between `origin/main` and HEAD (enforcing linear history). Contributors activate it once per clone with
`pre-commit install --hook-type pre-push`.

#### Design decisions

**1. Pre-commit framework rather than a standalone `.githooks/` directory.**

The repo already uses pre-commit for pre-commit stage hooks (ruff, pymarkdown, gitignore enforcement). Adding a
pre-push stage hook to the same `.pre-commit-config.yaml` keeps everything in one system. The alternative — checking
in a `.githooks/` directory and setting `core.hooksPath` — would create a parallel hook management system and conflict
with pre-commit's own `core.hooksPath` usage.

**2. Requires a separate install command: `pre-commit install --hook-type pre-push`.**

`pre-commit install` (without flags) only installs the `pre-commit` hook type. There is no single-command way to
install all hook types — each needs a separate `-t` invocation, and multiple `-t` flags in one call are not supported.
This is a pre-commit framework limitation. The install command is documented in `INSTALL.md` and `AGENTS.md`.

**3. Hook logic extracted to `scripts/check-rebase-on-main.sh`.**

The initial implementation inlined all logic as a bash one-liner in `.pre-commit-config.yaml`. This was hard to read,
test, and edit. Extracting to a standalone script referenced via `language: script` in pre-commit makes it maintainable
and directly testable (`bash scripts/check-rebase-on-main.sh`).

**4. Rejects merge commits — enforces linear history.**

The merge-base check alone passes for both `git rebase origin/main` and `git merge origin/main`. Since the policy
requires rebase (linear history), the hook additionally checks `git log --merges origin/main..HEAD` and rejects any
merge commits between origin/main and HEAD.

**5. Skips check on main branch.**

Pushing main itself (e.g., after a merge) should not be blocked. The hook exits 0 immediately when
`git rev-parse --abbrev-ref HEAD` is `main`.

**6. Fetches origin/main before checking.**

The hook runs `git fetch origin main --quiet` to ensure it checks against the latest remote state, not a stale local
ref. This adds a small network call but prevents false passes when origin/main has advanced since the last fetch.

#### PRs

- PR #193: hook implementation, AGENTS.md/INSTALL.md docs, CI workflow updates.

### 2026-03-22: CI token usage baseline — 100-run snapshot

#### Summary

First comprehensive token/cost analysis across all LLM-invoking CI workflows using the `ci-token-usage` skill.
Covers 100 most recent workflow runs (2026-03-21 to 2026-03-22).

#### Report

```text
Run ID       Workflow              Date        Status     Model    Cost       PR/Issue
-----------  --------------------  ----------  ---------  -------  ---------  ------------------
23393507574  Codex Implement Task  2026-03-22  failure    gpt-5.4  $1.23      TG04
23393126080  Codex Implement Task  2026-03-22  success    gpt-5.4  $1.36      TG03
23392903619  Codex Implement Task  2026-03-22  success    gpt-5.4  $0.96      TG02
23392501930  Codex Implement Task  2026-03-22  success    gpt-5.4  $1.46      TG01
23392240588  Codex Implement Task  2026-03-22  success    gpt-5.4  $1.29      TE03
23392054053  Codex Implement Task  2026-03-22  success    gpt-5.4  $1.17      TE02
23391859832  Codex Implement Task  2026-03-22  success    gpt-5.4  $1.51      TE01
23391685029  Codex Implement Task  2026-03-22  success    gpt-5.4  $0.66      TD03
23401448088  Codex PR Lifecycle    2026-03-22  skipped             skipped    PR#110
23401446130  Codex PR Lifecycle    2026-03-22  failure             ?          PR#111
23401099955  Codex PR Lifecycle    2026-03-22  failure             no LLM     PR#110
23400156019  Codex PR Lifecycle    2026-03-22  failure             no LLM     PR#108
23399190036  Codex PR Lifecycle    2026-03-22  failure             no LLM     PR#107
23399162621  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#107
23393491419  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#105 / Issue#104
23393109495  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#103 / Issue#102
23392886948  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#101 / Issue#100
23392484693  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#99 / Issue#98
23392445365  Codex PR Lifecycle    2026-03-22  skipped             skipped    PR#99 / Issue#98
23392445346  Codex PR Lifecycle    2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392445342  Codex PR Lifecycle    2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392408551  Codex PR Lifecycle    2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392222905  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#97 / Issue#96
23392037136  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#95 / Issue#94
23391842693  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#93 / Issue#92
23391664643  Codex PR Lifecycle    2026-03-22  success             no LLM     PR#91 / Issue#90
23401661800  Claude PR Review      2026-03-22  skipped             skipped    PR#111
23401564745  Claude PR Review      2026-03-22  skipped             skipped    PR#109
23401448139  Claude PR Review      2026-03-22  skipped             skipped    PR#110
23401446395  Claude PR Review      2026-03-22  skipped             skipped    PR#111
23401446310  Claude PR Review      2026-03-22  skipped             skipped    PR#111
23401433093  Claude PR Review      2026-03-22  skipped             skipped    PR#111
23401425335  Claude PR Review      2026-03-22  skipped             skipped    PR#110
23401377585  Claude PR Review      2026-03-22  skipped             skipped    PR#111
23401364528  Claude PR Review      2026-03-22  skipped             skipped    PR#110
23401099963  Claude PR Review      2026-03-22  skipped             skipped    PR#110
23401063163  Claude PR Review      2026-03-22  skipped             skipped    PR#110
23400672332  Claude PR Review      2026-03-22  skipped             skipped    PR#109
23400156013  Claude PR Review      2026-03-22  skipped             skipped    PR#108
23400156012  Claude PR Review      2026-03-22  skipped             skipped    PR#108
23400099350  Claude PR Review      2026-03-22  skipped             skipped    PR#108
23399190076  Claude PR Review      2026-03-22  skipped             skipped    PR#107
23399190041  Claude PR Review      2026-03-22  skipped             skipped    PR#107
23399143141  Claude PR Review      2026-03-22  skipped             skipped
23399138290  Claude PR Review      2026-03-22  success             $0.38      PR#107
23393461248  Claude PR Review      2026-03-22  skipped             skipped    TG03
23393459493  Claude PR Review      2026-03-22  cancelled           cancelled  TG03
23393458339  Claude PR Review      2026-03-22  success             $0.84      PR#105 / Issue#104
23393109763  Claude PR Review      2026-03-22  skipped             skipped    PR#103 / Issue#102
23393109483  Claude PR Review      2026-03-22  cancelled           cancelled  PR#103 / Issue#102
23393079767  Claude PR Review      2026-03-22  cancelled           cancelled  TG02
23393076800  Claude PR Review      2026-03-22  cancelled           cancelled  TG02
23393075508  Claude PR Review      2026-03-22  success             $0.69      PR#103 / Issue#102
23392886997  Claude PR Review      2026-03-22  skipped             skipped    PR#101 / Issue#100
23392886969  Claude PR Review      2026-03-22  cancelled           cancelled  PR#101 / Issue#100
23392886955  Claude PR Review      2026-03-22  cancelled           cancelled  PR#101 / Issue#100
23392834873  Claude PR Review      2026-03-22  cancelled           cancelled  TG01
23392831636  Claude PR Review      2026-03-22  cancelled           cancelled  TG01
23392830381  Claude PR Review      2026-03-22  success             $1.47      PR#101 / Issue#100
23392445476  Claude PR Review      2026-03-22  skipped             skipped    PR#99 / Issue#98
23392445394  Claude PR Review      2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392445372  Claude PR Review      2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392443487  Claude PR Review      2026-03-22  success             $0.40      PR#99 / Issue#98
23392408605  Claude PR Review      2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392408604  Claude PR Review      2026-03-22  skipped             skipped    PR#99 / Issue#98
23392408584  Claude PR Review      2026-03-22  cancelled           cancelled  PR#99 / Issue#98
23392367636  Claude PR Review      2026-03-22  cancelled           cancelled  TE03
23392365847  Claude PR Review      2026-03-22  cancelled           cancelled  TE03
23392364225  Claude PR Review      2026-03-22  success             $1.01      PR#99 / Issue#98
23392222717  Claude PR Review      2026-03-22  skipped             skipped    PR#97 / Issue#96
23392188368  Claude PR Review      2026-03-22  cancelled           cancelled  TE02
23392186591  Claude PR Review      2026-03-22  cancelled           cancelled  TE02
23392184746  Claude PR Review      2026-03-22  success             $0.87      PR#97 / Issue#96
23392037153  Claude PR Review      2026-03-22  skipped             skipped    PR#95 / Issue#94
23392003520  Claude PR Review      2026-03-22  cancelled           cancelled  TE01
23392001064  Claude PR Review      2026-03-22  cancelled           cancelled  TE01
23391999833  Claude PR Review      2026-03-22  success             $0.95      PR#95 / Issue#94
23391842750  Claude PR Review      2026-03-22  skipped             skipped    PR#93 / Issue#92
23391842735  Claude PR Review      2026-03-22  cancelled           cancelled  PR#93 / Issue#92
23391789152  Claude PR Review      2026-03-22  cancelled           cancelled  TD03
23391785627  Claude PR Review      2026-03-22  success             $1.15      PR#93 / Issue#92
23391466934  Codex Implement Task  2026-03-21  success    gpt-5.4  $0.85      TD02
23390906065  Codex Implement Task  2026-03-21  success    gpt-5.4  $1.45
23390432373  Codex Implement Task  2026-03-21  success    gpt-5.4  $1.85      TC04
23390282751  Codex Implement Task  2026-03-21  success    gpt-5.4  $0.96      TC03
23390122019  Codex Implement Task  2026-03-21  success    gpt-5.4  $0.74      TC02
23376412851  Codex Implement Task  2026-03-21  success    gpt-5.4  $1.09
23368200679  Codex Implement Task  2026-03-21  failure             no LLM     TI05
23367963808  Codex Implement Task  2026-03-21  failure    gpt-5.4  $0.83
23391593132  Codex PR Lifecycle    2026-03-21  success    gpt-5.4  $0.34      PR#91 / Issue#90
23391450667  Codex PR Lifecycle    2026-03-21  success             no LLM     PR#87 / Issue#21
23391423648  Codex PR Lifecycle    2026-03-21  skipped             skipped    PR#87 / Issue#21
23391423647  Codex PR Lifecycle    2026-03-21  cancelled           cancelled  PR#87 / Issue#21
23391345012  Codex PR Lifecycle    2026-03-21  cancelled  gpt-5.4  $0.39      PR#87 / Issue#21
23391247403  Codex PR Lifecycle    2026-03-21  success             no LLM     PR#86 / Issue#83
23391242099  Codex PR Lifecycle    2026-03-21  skipped             skipped    PR#87 / Issue#21
23391221579  Codex PR Lifecycle    2026-03-21  skipped             skipped    PR#86 / Issue#83
23391208905  Codex PR Lifecycle    2026-03-21  cancelled           cancelled  PR#87 / Issue#21
23391182628  Codex PR Lifecycle    2026-03-21  cancelled  gpt-5.4  $0.28      PR#86 / Issue#83
23390764607  Codex PR Lifecycle    2026-03-21  success             no LLM     PR#85

Total estimated cost: $26.18
  Codex:   $18.42  (18 runs, avg $1.02)
  Claude:  $7.76  (9 runs, avg $0.86)

  ⚠ Codex costs are estimates (blended 30% in / 70% out rate)
```

#### Interpretation

**Cost breakdown.** Total LLM spend across 100 runs: **$26.18** (Codex $18.42, Claude $7.76). Codex costs are
estimates using a 30/70 input/output blended rate against gpt-5.4 pricing; Claude costs are exact values reported by
`anthropics/claude-code-action`.

**Codex implementation runs** average $1.02 per run (range $0.66–$1.85). All use gpt-5.4. One failure (TI05) never
reached the LLM — correctly shows `no LLM`. One failure (run 23367963808) did consume tokens ($0.83) before failing.

**Claude review runs** average $0.86 per successful review (range $0.38–$1.47). Most expensive was PR#101/Issue#100
at $1.47. Many runs show `skipped` or `cancelled` — this is expected concurrency behavior when multiple pushes
trigger the review workflow in quick succession; only the latest run executes.

**Lifecycle runs** correctly split between `no LLM` (auto-merge jobs) and Codex-detected runs (address-feedback jobs
that invoke Codex, e.g. PR#91 at $0.34). Log-based detection handles this mixed-workflow case well.

**Concurrency waste.** PR#99/Issue#98 triggered 7 Claude review runs but only 2 completed ($0.40 + $1.01). The rest
were cancelled or skipped. Similarly, 4 lifecycle runs were cancelled for PR#87/Issue#21. This is not token waste
(cancelled runs do not consume LLM resources) but does consume GitHub Actions minutes.

#### Per-ticket drill-down

Two tickets verified with `--ticket`:

- **Issue #104 (TG03):** Codex implementation $1.36 (120,957 tok, gpt-5.4) + Claude review $0.84 (20 turns) =
  **$2.20 total**. Lifecycle run correctly shows `no LLM`.
- **Issue #98 (TE03):** Codex implementation $1.29 + 2 Claude reviews ($0.40 at 17 turns + $1.01 at 35 turns) =
  **$2.71 total**. Cancelled/skipped concurrency runs correctly show zero cost.

#### Tool used

Report generated with: `python -m lyzortx.orchestration.ci_token_usage --runs 100` and `--ticket 104` / `--ticket 98`
(see `.agents/skills/ci-token-usage/` for skill documentation and design decisions).

### 2026-03-22: Per-task model selection for Codex CI and orchestrator safety fix

#### Executive summary

Added per-task LLM model selection to the Codex orchestration pipeline. Each task in `plan.yml` now has a required
`model` field (`gpt-5.4` or `gpt-5.4-mini`) that flows through the orchestrator into CI workflows. 5 complex tasks
(SHAP, feature sweep, harmonization protocol, confidence tiers, external data integration) use `gpt-5.4`; 11
straightforward tasks (stats, visualization, parameterized loops, docs) use `gpt-5.4-mini` at ~70% lower cost. Projected
savings: ~48% on implementation runs (~$8.25 saved across the remaining 16 tasks). During rollout, discovered and fixed
a latent orchestrator bug where manually closing an issue would incorrectly mark its task as done — now gated on
GitHub's `state_reason` field.

#### Problem statement

All Codex implementation runs use `gpt-5.4` regardless of task complexity. The CI token usage baseline above shows 8
implementation runs consuming 856K tokens at an average of $1.07/run. Many tasks are straightforward (bootstrap CIs, bar
charts, parameterized training loops) and don't need the full model. `gpt-5.4-mini` costs ~70% less ($0.75/1M input,
$4.50/1M output vs $2.50/$15.00) and scores 54.4% on SWE-Bench Pro with a 400K context window — sufficient for this
repo's tasks.

#### Design decisions

**1. Model field lives in plan.yml, not in the workflow or issue template.**

The user runs a single planning conversation with subscription-model Claude Opus 4.6 to assign models to all pending
tasks at once. This avoids an additional LLM call at dispatch time and keeps model assignment as a reviewable data change
in version control. The `plan.yml` file already contains implementation-adjacent fields (`implemented_in`, `baseline`),
so `model` fits the existing schema pattern.

Alternative considered: parsing model from the issue body at dispatch time (set by the orchestrator based on heuristics).
Rejected because it requires either a second LLM call per dispatch or hand-coded heuristics that would be fragile.

**2. Model is a required field — no defaults, no fallbacks.**

`Task.model` is a required positional field on the dataclass (no default value). `load_pending_tasks()` and `run_once()`
both validate that every pending task has a non-empty model before proceeding, raising `ValueError` if any are missing.
The `codex-implement.yml` workflow also fails with `::error::` if the model directive is absent from the issue body.

This was a deliberate tightening during implementation. The initial version had `model: str = ""` with validation after
construction, but the user correctly identified that a default empty string is a silent fallback that defeats the purpose.
Making it required at the type level means you cannot construct a `Task` without specifying a model, which is enforced by
Python's `TypeError` on missing arguments to frozen dataclasses.

**3. Model travels as an HTML comment in the issue body: `<!-- model: gpt-5.4-mini -->`.**

The orchestrator emits this directive when creating issues. Both `codex-implement.yml` and `codex-pr-lifecycle.yml`
extract it using `parse_model_directive.py`, a small pure-function module that reads stdin and prints the model ID.

Why an HTML comment rather than a visible field: the directive is machine-readable metadata, not something the
implementing agent needs to see or act on. HTML comments don't render in GitHub's issue view, keeping the issue body
clean for human readers.

Why a dedicated Python module rather than inline `grep -oP`: PCRE lookbehind (`grep -P`) portability is uncertain across
CI runner images. A 12-line Python script is portable, testable, and reusable from both workflows.

**4. Lifecycle (review feedback) runs use the same model as the original implementation.**

`codex-pr-lifecycle.yml` extracts the linked issue number from the PR body (`Closes #N`), fetches that issue's body, and
extracts the model directive from it. This ensures consistency: if TG04 was implemented with `gpt-5.4`, all review
feedback rounds also use `gpt-5.4`.

Alternative considered: always use the cheaper model for feedback rounds since they're typically simpler fixes. Rejected
for now — consistency is more important until we have data on whether mini handles feedback adequately.

**5. Two-PR rollout: data first, then code.**

PR #109 (merged) added the `model` field to all 16 pending tasks in `plan.yml` as a pure data change. The existing
`load_plan()` function ignores unknown YAML fields, so this was a no-op. PR #112 wires the field through the system.
This separation allows model assignments to be reviewed and adjusted independently of code changes, and ensures a clean
rollback path if the wiring has issues.

#### Model assignments

5 tasks assigned `gpt-5.4` (complex reasoning, architectural design, domain interpretation):

<!-- pyml disable-num-lines 7 md013-->
| Task | Rationale |
|------|-----------|
| TG04 — SHAP explanations | TreeExplainer integration (new to codebase), cross-referencing ablation + model outputs, prescriptive recommendation of which feature blocks to keep |
| TG05 — Feature-subset sweep | 10 combinatorial model runs, comparison logic, locks final v1 feature config for all downstream tracks |
| TI05 — Harmonization protocol | Multi-source schema design, domain-critical decisions cascading to TI06–TI10 |
| TI07 — Confidence tiers | Subjective tier design + weighting strategy, cascades to TI08/TI09 |
| TI08 — External data integration | Conditional injection architecture, leakage prevention, fallback handling |

11 tasks assigned `gpt-5.4-mini` (stats, visualization, parameterized loops, docs):

<!-- pyml disable-num-lines 13 md013-->
| Task | Rationale |
|------|-----------|
| TF01 — Bootstrap CIs | Standard NumPy resampling, dual-slice filtering already exists in codebase |
| TF02 — v0 vs v1 comparison | Side-by-side metric table, algorithmic error bucket identification |
| TH02 — Explained recommendations | Data assembly from TG02+TG04 outputs, formatting |
| TI06 — Tier B ingestion | ID cross-referencing, lookup joins, follows TI04 patterns |
| TI09 — Sequential ablations | Parameterized loop over TG01 training, metric collection |
| TI10 — Lift tracking | GroupBy aggregation + failure mode detection |
| TJ01 — One-command regeneration | Orchestration script calling existing run_track_*.py |
| TJ02 — Environment freeze | Documentation + version pinning |
| TP01 — Digital phagogram | Plotly/Matplotlib visualization |
| TP02 — Panel coverage heatmap | Standard Seaborn heatmap |
| TP03 — Feature lift bar chart | Bar chart from existing TG03 CSV |

Each task was evaluated on four axes: novelty (new pattern vs reuse), domain criticality (does a mistake cascade?),
reasoning depth (multi-step logic vs straightforward assembly), and established patterns (can it largely copy TG01/TG03
structure?).

#### Projected cost impact

Assuming average ~107K tokens per run (observed today):

- **Status quo (all gpt-5.4):** 16 × $1.07 ≈ $17.12
- **With model selection:** (5 × $1.07) + (11 × $0.32) ≈ $5.35 + $3.52 ≈ $8.87
- **Estimated savings:** ~48% on implementation runs

This does not include lifecycle runs, which scale proportionally.

#### Orchestrator safety fix: state_reason gate for issue closure

During implementation, we closed 4 pre-existing orchestrator-task issues (#22, #25, #70, #106) that predated model
selection. This exposed a latent bug: the orchestrator's `sync_status_from_issues` treated any closed issue as
"completed," which would incorrectly mark unfinished tasks as done in `plan.yml`.

GitHub's REST API provides `state_reason` on issues: `"completed"` (closed by PR merge or explicitly marked done) vs
`"not_planned"` (manually closed without completion). The fix:

- Added `state_reason` field to the `IssueRef` dataclass.
- `list_task_issues()` now parses `state_reason` from the API response.
- `sync_status_from_issues()` only marks a task "completed" when `state_reason == "completed"`. Issues closed as
  `"not_planned"` revert to their previous status (pending or blocked), preserving the task's position in the dispatch
  queue.

This is a correctness fix independent of model selection — the bug existed before but was never triggered because issues
were only closed by PR merges (which set `state_reason: "completed"`). The model-selection migration was the first time
issues were manually closed for housekeeping.

Tests added: 5 test cases covering all `state_reason` paths (completed, not_planned, open, blocked preservation,
no-issue fallback).

#### PRs

- PR #109 (merged): data-only — added `model` field to all 16 pending tasks in `plan.yml`.
- PR #112: wired model through orchestrator, workflows, and added state_reason safety gate.

### 2026-03-24: Orchestrator task invalidation — reopening issues triggers spurious runs

#### Executive summary

Reopening 13 closed orchestrator issues (TI03-TI10, TK01-TK05) to change their close reason from "completed" to "not
planned" accidentally triggered 13 Codex implementation runs. The `issues.reopened` event fires the implementation
workflow before the issue can be re-closed. The safe workaround is to change `state_reason` via the GitHub API without
reopening: `gh api repos/OWNER/REPO/issues/NUMBER -X PATCH -f state=closed -f state_reason=not_planned`.

#### Design note

The root issue is that the orchestrator uses issue state changes as its trigger. Reopening an issue is indistinguishable
from a new dispatch signal. A more robust design would trigger on merged PRs rather than closed issues — PR merges are
unambiguous completion signals and cannot be accidentally triggered by state changes on issues.

### 2026-03-24: Claude PR review approvals are nondeterministic — formal reviews silently skipped

#### Executive summary

`claude-code-action@v1`'s built-in system prompt tells the agent "You CANNOT submit formal GitHub PR reviews" and "You
CANNOT approve pull requests (for security reasons)." This directly contradicts our custom prompt that instructs Claude
to use `mcp__github__submit_pending_pull_request_review`. The result is nondeterministic: sometimes the agent follows
our instructions and submits a formal review (PRs #217, #220, #223), sometimes it follows the built-in restriction and
writes its entire review into the sticky issue comment instead (PR #225). When this happens, the `claude-pr-review.yml`
dispatch step sees no formal review and takes no action — no auto-merge, no Codex lifecycle dispatch.

#### Root cause analysis

The conflicting instructions are in `anthropics/claude-code-action` source file
[`src/create-prompt/index.ts`](https://github.com/anthropics/claude-code-action/blob/v1/src/create-prompt/index.ts)
(the `generateDefaultPrompt()` function, "CAPABILITIES AND LIMITATIONS" section):

> ```text
> What You CANNOT Do:
> - Submit formal GitHub PR reviews
> - Approve pull requests (for security reasons)
> ```

These lines were present in the [initial commit `f66f337f`](https://github.com/anthropics/claude-code-action/commit/f66f337f)
on 2025-05-19 and have never been modified. The action provides the MCP review tools in its `allowedTools` list, but the
prompt-level "CANNOT" instruction causes the model to sometimes refuse to use them.

Evidence from PR #225 logs:
- `permission_denials_count: 3` — the agent attempted to use the review tools and was denied.
- `No buffered inline comments` — the `classify_inline_comments` feature captured nothing.
- `Unexpected review state: . No action taken.` — the dispatch step found zero formal reviews.

The contradiction is a prompt-level issue, not a tool-permission issue. The MCP tools are available and allowed; the
model simply refuses to call them when it prioritizes the "CANNOT" instruction over the custom prompt.

#### Impact on downstream workflows

When Claude fails to submit a formal review:
1. `claude-pr-review.yml` dispatch step sees empty `LATEST_REVIEW` → takes the `else` branch → no auto-merge, no Codex
   dispatch.
2. `codex-pr-lifecycle.yml` (if manually triggered) queries `review_threads.py` which only reads formal review threads
   → finds zero unresolved threads → labels PR `ready-for-human-review` and posts "Review passed with no issues" even
   though actionable feedback exists in the issue comment.

#### Workarounds investigated

- `USE_SIMPLE_PROMPT=true` env var: switches to a simplified prompt that omits the CANNOT section. Stays in tag mode.
  Risk: the simplified prompt may lose other useful scaffolding.
- Stronger override language in custom prompt: fragile — sometimes works, sometimes doesn't, as PR #225 demonstrated.
- Post-condition guard in the workflow: check whether a formal review was submitted after the action completes; if not,
  fail loudly or submit a synthetic review via `gh api`. This is a band-aid.

#### Long-term fix: LangGraph review orchestrator (PR #48)

The correct fix is to wrap the review agent in a LangGraph orchestrator that can verify tool calls were actually made
and loop the agent back with feedback if the formal review submission is missing. This moves the reliability guarantee
from prompt-level persuasion (unreliable) to programmatic verification (deterministic). PR #48 contains the
implementation plan.

#### Immediate actions taken

- Updated `/gh` skill (section 10) to document that PR feedback can appear in either review threads or issue comments,
  and that both must be checked.
- PR #225 was merged manually after human review confirmed Claude's comment-based review was accurate.

#### Future: migrate orchestrator trigger from issues to PRs

Revisit when: the current trigger model causes more incidents, or the orchestrator is being refactored for other
reasons. The change would involve updating `codex-implement.yml` to trigger on `pull_request.closed` (with
`merged == true`) instead of `issues.reopened`/`issues.opened`, and updating `sync_status_from_issues` to read from PR
merge state instead of issue close state.

#### Future: ticket dependency support in plan.yml with GitHub issue relationships

Revisit when: the orchestrator is being refactored or task scheduling becomes more complex.

Currently `plan.yml` tasks have a flat `status` field (pending, blocked, done) with no structured dependency information.
Task ordering is implicit — humans and agents infer dependencies from track structure and naming conventions (e.g., TI05
depends on TI04 because it follows in the track). This is fragile and prevents the orchestrator from automatically
determining which tasks are unblocked when a predecessor completes.

The improvement would add a `depends_on` field to each task in `plan.yml` listing predecessor task IDs. The orchestrator
would use this to:
1. Automatically set tasks to `blocked` when predecessors are not yet `done`.
2. Unblock tasks when all predecessors complete.
3. Reflect dependencies as GitHub issue relationships (sub-issues or "blocked by" links) so the dependency graph is
   visible in the GitHub UI, not just in the YAML file.

This would eliminate manual `blocked`/`pending` status management and make the dispatch order deterministic based on the
dependency graph rather than YAML ordering.

#### Future: implementation workflow should comment on the issue with a link to the action run

Revisit when: improving observability of the orchestrator pipeline.

Currently when `codex-implement.yml` picks up an issue and starts an implementation run, there is no trace on the issue
itself. A human watching the issue has no way to know that work is in progress or where to find the logs — they have to
go to the Actions tab and hunt for the right run.

The workflow should post a comment on the dispatching issue at the start of the run, e.g.:
"Implementation started — [View action run](https://github.com/OWNER/REPO/actions/runs/RUN_ID)."

This is a single `gh issue comment` step early in the job, using `${{ github.run_id }}` to construct the URL. It gives
immediate visibility into which issues are being worked on and provides a direct link to the logs for debugging when a
run fails.

#### Future: codex-pr-lifecycle should register as a check on the PR

Revisit when: improving PR observability or refactoring the lifecycle workflow.

`codex-pr-lifecycle.yml` is triggered via `workflow_dispatch`, which means it does not automatically appear in the PR's
checks tab. A human looking at the PR sees only the CI checks from `push`/`pull_request`-triggered workflows — there is
no indication that a Codex fix cycle is running, succeeded, or failed.

The workflow should use the GitHub Checks API (`gh api repos/OWNER/REPO/check-runs`) to create a check run linked to
the PR's head SHA at the start of the job, update it with progress, and mark it completed (success/failure) at the end.
This makes the lifecycle run visible alongside CI checks in the PR's "Checks" tab, so reviewers can see at a glance
whether the Codex fix cycle has completed.

The check run should be created early in the job using:
```bash
gh api repos/OWNER/REPO/check-runs -f name="Codex PR Lifecycle" \
  -f head_sha="$(git rev-parse HEAD)" -f status="in_progress" \
  -f details_url="https://github.com/OWNER/REPO/actions/runs/${{ github.run_id }}"
```
And updated to `completed` with the appropriate conclusion at the end.

### 2026-03-24: GitHub Actions now bootstrap `phage_env` with conda

#### Decision

Updated `codex-implement.yml`, `codex-pr-lifecycle.yml`, and `claude-pr-review.yml` to create `phage_env` from
`environment.yml` using the runner's preinstalled Miniconda, then invoke repo Python and pre-commit commands via
`conda run -n phage_env ...` instead of exporting env-specific paths.

Also pinned the Claude review jobs from `ubuntu-latest` to `ubuntu-24.04`. Relying on a preinstalled runner package is
too brittle on the floating `ubuntu-latest` label.

This supersedes the earlier sketch above that still assumed a separate `pip install -r requirements.txt` step after
conda bootstrapping. `environment.yml` is now the canonical CI environment definition and already includes the pip
requirements.

#### Source evidence

- Ubuntu 24.04 runner image README:
  https://github.com/actions/runner-images/blob/main/images/ubuntu/Ubuntu2404-Readme.md
  Quote: `Miniconda 26.1.1` is listed under "Package Management" and `CONDA /usr/share/miniconda` is listed under
  "Environment variables".
- GitHub Actions changelog:
  https://github.blog/changelog/2024-09-25-actions-new-images-and-ubuntu-latest-changes/
  Quote: "`ubuntu-latest` label will migrate to Ubuntu 24 over the course of the next month".

#### Rationale

- Keeps GitHub Actions aligned with local development and `.envrc`, both of which now target `conda activate
  phage_env`.
- Ensures Codex and Claude runs see the same Python toolchain and pip-installed packages as the declared environment,
  instead of a partially reconstructed `pip install -r requirements.txt` subset.
- Avoids depending on reconstructed env prefixes or manual `PATH` injection; the workflows only assume the `conda`
  command is available on the pinned runner image.

### 2026-03-24: Pharokka POC and TL01 prereqs

#### Executive summary

Track L (TL01) will add bioconda tools (Pharokka, MMseqs2, tRNAscan-SE, MinCED, ARAGORN, mash, dnaapler) to
`environment.yml`. The conda CI bootstrap above already handles `conda env create --file environment.yml`, so TL01
only needs to add the bioconda dependencies to the YAML — no workflow changes required.

#### Open item: Pharokka database caching

Pharokka databases are ~656MB. They need to be either downloaded in CI per run or cached via GitHub Actions cache.
TL01 should address this.

#### POC results (local, 2026-03-24)

- Pharokka 1.9.1 installed via pip, bioconda deps via micromamba
- Pharokka databases: 656MB, downloaded in ~66s
- Single phage annotation (LF82_P8): 2m 43s, 276 CDS annotated
- Produced: 29 tail genes, 7 lysis genes, 1 anti-restriction nuclease, 140 hypothetical proteins
- Extrapolation for 97 phages: ~4-5 hours sequential, ~1 hour at 4 threads

### 2026-03-28: Add code coverage reporting to CI via Codecov

#### Executive summary

Added `pytest-cov` and Codecov integration to the unit-tests workflow so that PR diffs show line-level coverage
annotations. This gives contributors immediate feedback on which new or changed lines lack test coverage, without
leaving the PR review UI.

#### Design decisions

**1. Codecov over alternatives.** Evaluated Codecov, Coveralls, and manual coverage-comment actions. Codecov was chosen
because it provides line-level annotations directly in PR diffs — per Codecov docs: "Codecov will [...] display results
in all future pull requests" ([source](https://docs.codecov.com/docs/pull-request-comments)). Free for open-source
repos ([pricing](https://about.codecov.io/pricing/): "Free — Open Source — $0/month").

**2. `fail_ci_if_error: false`.** The Codecov upload step is non-blocking. Coverage reporting is informational; a
transient Codecov outage should not fail an otherwise-green CI run. Per `codecov-action` README: "Specify if the CI
pipeline should fail when Codecov runs into errors during upload" ([source](https://github.com/codecov/codecov-action)).

**3. Coverage scope.** `.coveragerc` measures `lyzortx/` source, omitting `lyzortx/tests/` and
`lyzortx/research_notes/` from coverage metrics. This focuses the signal on production code.

**4. SHA-pinned action + explicit token.** The `codecov-action` is pinned to a full commit SHA (`75cd116...` = v5.5.4)
for reproducibility. An explicit `CODECOV_TOKEN` secret is used for authentication — tokenless uploads can silently
fail in some configurations per the action docs: "Required for [...] uploading coverage to Codecov"
([source](https://github.com/codecov/codecov-action#arguments)).

**5. PR comment reporting.** Added `codecov.yml` with `comment.layout: "reach,diff,flags,components"` to enable
coverage delta summaries on every PR, per Codecov's common recipes
([source](https://docs.codecov.com/docs/common-recipe-list#show-project-coverage-changes-on-the-pull-request-comment)).

#### Prerequisites

The Codecov GitHub App must be installed on the repo for annotations to appear on PRs. A `CODECOV_TOKEN` repository
secret must be configured (available from the Codecov dashboard after app installation).

### 2026-03-30: Replace personal PAT with Czarphage GitHub App for orchestration auth

#### Executive summary

Replaced the `ORCHESTRATOR_PAT` secret (a long-lived personal access token belonging to a human account) with the
**Czarphage** GitHub App (App ID 3227910) across all four orchestration workflows. Each job now mints a short-lived
installation token via `actions/create-github-app-token@v2`.

#### Motivation

The orchestration pipeline (`orchestrator.yml`, `codex-implement.yml`, `codex-pr-lifecycle.yml`, `claude-pr-review.yml`)
authenticates with a PAT to bypass GitHub's anti-recursion protection (the default `github.token` cannot trigger other
workflows). Using a personal PAT ties automation to an individual account, creates audit ambiguity (commits and API
calls appear as the PAT owner), and requires manual token rotation.

A GitHub App provides: short-lived auto-rotating tokens (1-hour expiry), a distinct `czarphage[bot]` identity in
commits and API calls, scoped permissions without repo-wide personal access, and no dependency on any individual's
account.

#### Design decisions

**1. Token generation as the first step in each job.** The `actions/create-github-app-token@v2` step runs before
checkout because `actions/checkout` also needs the token. The step requires no repo access — it authenticates directly
with GitHub's API using the app's private key.

**2. Per-step `GITHUB_TOKEN` env instead of job-level.** Job-level `env` blocks cannot reference step outputs (`${{
steps.*.outputs.* }}`). Rather than introducing a workaround, each step that needs the token declares it in its own
`env` block. This is slightly more verbose but explicit about which steps are authenticated.

**3. Git identity set to `czarphage[bot]`.** The noreply email format for GitHub Apps is
`<APP_ID>+<app-slug>[bot]@users.noreply.github.com`. Commits now show as authored by `czarphage[bot]` rather than
`github-actions[bot]`, making it clear which automation system produced them.

**4. Token lifetime vs workflow timeout.** GitHub App installation tokens expire after 1 hour. The longest workflow
timeout is `codex-implement.yml` at 45 minutes, leaving a 15-minute margin. If this proves tight, the token step can
be repeated mid-job — `actions/create-github-app-token` is idempotent and fast (~2s).

#### Changes

- `orchestrator.yml`: Added token step, moved `GITHUB_TOKEN` to per-step env, updated git identity.
- `codex-implement.yml`: Same pattern, also updated the agent prompt text (removed PAT reference).
- `codex-pr-lifecycle.yml`: Same pattern across all conditional steps.
- `claude-pr-review.yml`: Added token step for the dispatch step only (the rest uses `github.token` or Claude's OIDC).
- `.github/workflows/AGENTS.md`, `lyzortx/orchestration/README.md`: Removed `ORCHESTRATOR_PAT` references from docs.

#### Post-merge cleanup

After verifying the full orchestration loop works:
1. Delete the `ORCHESTRATOR_PAT` secret from repo settings.
2. Revoke the old PAT on the original account.

### 2026-03-30: Allow Czarphage bot-authored PR updates through Claude review workflow

#### Executive summary

Fixed `claude-pr-review.yml` so `anthropics/claude-code-action@v1` accepts PR events initiated by the repo's own
GitHub App (`czarphage` / `czarphage[bot]`). Without this allowlist entry, Claude auto-review failed on synchronize
events after Codex pushed follow-up commits.

#### Design decisions

**1. Narrow allowlist, not wildcard bots.** Set `allowed_bots: "czarphage,czarphage[bot]"` on the auto-review step
instead of `'*'`. This keeps the workflow scoped to the one automation identity that legitimately pushes PR updates in
this repo, rather than allowing arbitrary third-party bot-authored prompts.

**2. Fix only the auto-review path.** The failure occurred on `pull_request` synchronize events for orchestrator-task
PRs, so the change was limited to the `auto-review` job. The interactive `@claude` path was left unchanged because it
does not need bot-triggered invocation.

**3. Document the bot handoff explicitly.** Updated `lyzortx/orchestration/README.md` to say that Claude re-reviews
Codex follow-up pushes from the `czarphage` GitHub App. That behavior is otherwise non-obvious when reading the
workflow list in isolation.

#### Failure evidence

GitHub Actions run `23767069559`, job `69250897074` failed inside `anthropics/claude-code-action@v1` with:
"Workflow initiated by non-human actor: czarphage (type: Bot). Add bot to allowed_bots list or use '*' to allow all
bots."

### 2026-03-31: Prebake a Codex CI image for the implement/fix loop

#### Executive summary

Started replacing per-run Conda bootstrap in `codex-implement.yml` and `codex-pr-lifecycle.yml` with a prebaked GitHub
Container Registry image. The image is built from `.github/ci/Dockerfile` and includes a lightweight Python/helper
`phage_env`, a dedicated `phage_annotation_tools` env for heavy Pharokka/MMseqs2/OpenJDK-style CLIs, and three split
typing caller envs (`phylogroup_caller`, `serotype_caller`, `sequence_type_caller`).

#### Design decisions

**1. Cut over the Codex workflows first.** These two jobs were paying the biggest repeated `conda env create` cost and
are not pull-request-triggered, so they can move to a prebaked image without introducing a PR bootstrap dependency on a
new image publish.

**2. Keep `phage_env` small and move heavy CLIs out.** The default env for Codex workflows only needs Python, tests,
pre-commit, and repo helper modules. The Pharokka/MMseqs2/OpenJDK toolchain now lives in a dedicated
`phage_annotation_tools` env, and the raw typing tools live in their own caller envs. This avoids dragging the whole
bioinformatics stack into every helper solve.

**3. Refresh env manifests on job start.** Even though the image is prebaked, each run executes env refreshes from the
checked-in `environment.yml` and dedicated env manifests. That keeps dependency changes from waiting on an image rebuild
and matches the repo policy that checked-in env manifests remain the source of truth.

**4. Keep transient workflow files under `.scratch/ci/`.** Prompt and feedback files moved out of `/tmp` so the
containerized jobs and action steps share workspace-visible paths and the repo stays consistent with its own scratch
space policy.

**5. Treat the rollout as a cutover, not a compatibility promise.** The new Codex workflows require explicit
`ci-image:*` labels and the checked-in env manifests that match the selected profile. Older orchestrator issues/PRs
created before that contract are intentionally unsupported; they must be re-dispatched or rebased rather than routed
through fallback behavior.

**6. Keep specialized manifests at repo root.** The split caller/toolchain env manifests were moved alongside
`environment.yml` instead of living under `.github/` so local developers can create the same envs on bare metal. They
are CI inputs, but they are not CI-only.

#### Validation

- Verified GHCR package tags for the branch cutover images using the package API and a token with `read:packages`:
  - `base-branch-worktree-ci-image-cutover`
  - `host-typing-branch-worktree-ci-image-cutover`
  - `full-bio-branch-worktree-ci-image-cutover`
- Verified the corresponding SHA tags from the initial branch image line also exist:
  - `base-sha-7231a81bb21eddd0e34095e4613506220774fde7`
  - `host-typing-sha-7231a81bb21eddd0e34095e4613506220774fde7`
  - `full-bio-sha-7231a81bb21eddd0e34095e4613506220774fde7`

### 2026-04-02: Dead code pruning — criterion selection and first-attempt failure

#### Executive summary

Pruned ~9.2k lines of dead-end code from `lyzortx/` (33 files). The first attempt used `pydeps` transitive import
closure from the TL18 entry point as the reachability criterion, which deleted ~15k lines (86 files) including active CI
gates, the Track J regeneration pipeline, ad-hoc analysis scripts referenced by notebooks, and the
steel-thread-regression workflow. An independent review caught the over-deletion. The correct criterion turned out to be
"has no future use and no traceability value," applied case-by-case.

#### What was deleted

- **Track H** (cocktail recommendations): built on the leaked model, never updated post-leakage-fix. Would need
  rebuilding from scratch on the DEPLOY model anyway.
- **Track I** (external data ingestion): proven dead end — zero joinable external rows with the internal 404×96 panel.
  Documented in `project.md`. Download infrastructure would need rewriting for any future attempt.
- **ST v0 superseded experimental steps**: st03b (split suite), st04b (ablation suite), st06b (ranking policy
  comparison), st08 (tier A ingest ablation), plus the st03b regression check and baseline. All superseded by Track G.
- **Track A dead-end experiments**: `build_mechanistic_proxy_features`, `run_phistruct_rbp_pilot`. Findings recorded in
  lab notebooks; code not needed to reproduce current pipeline.
- **Track G analysis-only steps**: `investigate_non_leaky_candidate_features`, `run_feature_subset_sweep`. Their
  findings are locked in the lab notebook and `v1_feature_configuration.json`. Not needed to reproduce the current
  model.
- Tests for all of the above.

#### What was kept (and why the first attempt was wrong to delete it)

| Code | Why it must stay |
|------|-----------------|
| ST v0 regression checks + baselines | Active CI gate (`steel-thread-regression.yml`). Validates the data pipeline hasn't regressed. Needed when DEPLOY retrains from re-derived features. |
| `steel-thread-regression.yml` | The workflow that runs the checks above. |
| Track J (`run_track_j.py`) | One-command regeneration of all v1 outputs. Needed for DEPLOY track and any future release. |
| Track G runner (`run_track_g.py`) | Dispatches the 4 kept modeling steps. Without it, there's no CLI to run the modeling pipeline. |
| Ad-hoc analysis scripts | Referenced by lab notebook entries. Deleting them breaks traceability. Cheap to keep. |
| Track A/G runners, checks, READMEs | Infrastructure for tracks whose core code is kept. |

#### Design decision: why "pydeps import closure" was the wrong criterion

The `pydeps` transitive closure from `run_track_l.py` captures everything TL18 *imports at module level*. But
"imported by TL18" ≠ "needed by the project." The project also needs:

1. **CI validation** — regression checks that don't appear in any import chain because they're invoked via CLI
2. **Regeneration** — Track J's `run_track_j.py` orchestrates a full pipeline run but isn't imported by any step
3. **Traceability** — ad-hoc scripts referenced by notebook entries that document design decisions
4. **Future hooks** — code that is dormant now but has a documented "revisit when" trigger

The correct criterion is case-by-case: "has no future use AND no traceability value." This requires judgment that a
static import graph cannot provide.

#### Cascade updates required

Deleting steps that are imported by kept runners creates broken imports. Each deletion required updating the runner
that dispatches it:

- `run_steel_thread_v0.py` — removed st03b/st04b/st06b/st08 imports and dispatch cases
- `run_track_g.py` — removed `run_feature_subset_sweep` and `investigate_non_leaky_candidate_features`
- `run_track_j.py` — removed Track H import/dispatch (Track H was deleted) and the "recommendations" step choice

Tests for all three runners were updated to match the new dispatch tables.

#### Lesson

Automated reachability analysis is a useful starting point for identifying deletion candidates, but it must not be used
as the sole decision criterion. Always have a human or independent reviewer check the deletion list against project
needs that aren't captured in import graphs: CI gates, CLI entry points, documentation references, and planned future
work.

### 2026-04-02: CI runner sizing for DEPLOY track

#### Decision

Try the free `ubuntu-24.04` runner first for DEPLOY tasks. Fall back to the paid `ubuntu-24.04-4core` (150GB disk,
$0.008/min) only if the free runner's 14GB disk proves insufficient.

#### Context

The `full-bio` CI image is 9.5GB compressed. The free runner has 14GB disk with ~6GB usable after the runner OS. The
image already exceeds this, relying on layer sharing with the runner base to fit. Adding the 1.9GB Picard assembly
set to the image is not viable on the free runner.

Instead, DEPLOY tasks download assemblies at runtime from figshare (~7 min for the 1.9GB zip). This eats into the
90-minute task timeout but avoids image bloat. The download function skips if assemblies are already present, so local
dev only pays the cost once.

If the free runner runs out of disk during DEPLOY02 (DefenseFinder on 403 genomes produces intermediate files),
switching to `ubuntu-24.04-4core` costs ~$0.24-0.48 per task run (30-60 min typical). The full DEPLOY track would cost
~$3-5 total on paid runners. This is acceptable if needed but not worth pre-committing to.

#### Future: Switch pharokka to meta mode for batch annotation

**Trigger:** Next time pharokka annotations need regenerating, or if `run_pharokka.py` is refactored for the DEPLOY
track.

**Finding (2026-04-02):** Running pharokka in `--meta --split` mode on a single concatenated multi-FASTA of all 97 phage
genomes completes in **3 minutes 13 seconds** (192s wall time, 4 threads). The current per-phage approach
(`ProcessPoolExecutor` with 97 separate `pharokka.py` subprocess calls) takes **1-2 hours** with parallelism. The
speedup is ~20-40x because mmseqs2 database indexing and profile search happen once instead of 97 times.

The meta mode output is a single combined `pharokka_cds_final_merged_output.tsv` (3.4MB, 10,114 CDS rows) with a
`contig` column containing the phage name. The `--split` flag also produces per-phage GFF/FAA/FASTA files under
`single_gffs/` etc. Downstream code currently reads per-phage `_cds_final_merged_output.tsv` files; it would need
either a trivial split-by-contig post-processing step or a parser update to accept the combined file.

The combined TSV compresses to 343KB. It could be committed directly or baked into the CI image with negligible
size impact. The existing 194 committed files in `data/annotations/pharokka/` (4.5MB total) could be replaced by
this single file, though at 4.5MB the current approach is also not a real problem.

**What to do if triggered:** Replace the `ProcessPoolExecutor` loop in `run_pharokka.py` with a single
`pharokka.py --meta --split` call on a concatenated FASTA, then split the combined TSV by the `contig` column to
produce per-phage files compatible with existing downstream parsers.

### 2026-04-04 22:26 UTC: PR lifecycle should scan visible PR feedback, not unresolved-thread state

#### Executive summary

The lifecycle workflow was still carrying an older assumption from when it also owned the merge decision: it counted
unresolved review threads to decide whether a PR was clean. That is now the wrong abstraction. `claude-pr-review.yml`
either merges the PR itself or dispatches `codex-pr-lifecycle.yml`, so the lifecycle job only needs to answer one
question: "is there any visible feedback on the PR that Codex should read before saying nothing to do?"

#### Trigger

GitHub Actions run `23988591888` on PR #331 answered that question incorrectly. The job printed `Unresolved threads: 0,
Claude comments with feedback: 0`, then posted "Review passed with no issues" even though the PR already contained a
top-level Claude review comment and a formal `COMMENTED` review body with actionable feedback.

#### Design decision

Simplified the lifecycle contract:

- `claude-pr-review.yml` remains the merge gate.
- `codex-pr-lifecycle.yml` now treats visible PR feedback surfaces as its source of truth.
- The workflow reads top-level PR comments, inline review comments, and non-empty review bodies from non-`czarphage`
  authors.
- It should not re-derive PR cleanliness from unresolved-thread state, because that confuses "merge eligibility" with
  "feedback exists somewhere on the PR."

This change keeps the lifecycle aligned with the current responsibility split: Claude decides approve/merge vs dispatch;
Codex only addresses whatever feedback is already visible on the PR.

### 2026-04-05 20:25 UTC: AUTORESEARCH RunPod workflow moved cloud spend behind a dedicated environment gate

#### Executive summary

Added `.github/workflows/autoresearch-runpod.yml` as a manual-only workflow for AUTORESEARCH GPU search. The workflow
stages a frozen bundle, uses the dedicated `runpod-autoresearch` GitHub environment with `RUNPOD_API_KEY`, provisions
one fixed RunPod pod, uploads a candidate artifact bundle, and tears the pod down. This keeps paid cloud credentials
and provisioning logic out of the generic Codex issue/PR automation path.

### 2026-04-08 22:00 UTC: Added `/sleeponit` skill and knowledge model infrastructure

#### Executive summary

New skill `/sleeponit` and supporting Python infrastructure for consolidating lab notebook knowledge into a structured
YAML knowledge model. Follows the `plan.yml` → `PLAN.md` pattern: `knowledge.yml` → `knowledge_parser.py` →
`render_knowledge.py` → `KNOWLEDGE.md`. The rendered output is wired into Claude's context via `lyzortx/CLAUDE.md`.

#### Components

- `lyzortx/orchestration/knowledge_parser.py` — frozen dataclasses (`KnowledgeUnit`, `KnowledgeTheme`,
  `KnowledgeModel`, `KnowledgeDiff`), YAML loader, validator, and diff function.
- `lyzortx/orchestration/render_knowledge.py` — renders knowledge.yml to markdown with pymarkdown lint fixing. CLI:
  `python -m lyzortx.orchestration.render_knowledge [--knowledge-path PATH] [--output PATH] [--stdout]`.
- `.agents/skills/sleeponit/SKILL.md` — 3-phase skill: extract (read notebooks), curate (user reviews draft), emit
  (write YAML + render markdown + wire context).
- `lyzortx/tests/test_knowledge_parser.py` — 15 unit tests covering parser, all validation rules, diff detection, and
  render output.

#### Design decisions

- YAML source of truth (not markdown) — enables Python validation, diffing, and programmatic manipulation.
- Knowledge units carry: `id`, `statement`, `sources`, `status`, `confidence`, `context`, `relates_to`.
- `relates_to` provides lightweight cross-references between units without requiring a full graph database.
- Merge-conflict strategy: YAML is append-mostly; after merging, re-run the renderer to regenerate KNOWLEDGE.md.

### 2026-04-19 19:25 CEST: Added `review-ml-pr` skill; self-review subagent now required to run it on ML PRs

#### Executive summary

New skill `.agents/skills/review-ml-pr/` captures the verification protocol for ML PRs that produce prediction or
metric artifacts. The reviewer subagent spawned per `lyzortx/orchestration/AGENTS.md` now MUST invoke this skill as
part of priority (2) scientific correctness whenever the PR writes under `lyzortx/generated_outputs/` or cites AUC /
Brier / bootstrap CI numbers. A read-only review of an ML results PR is no longer acceptable. The skill was distilled
from an audit of the CH05 PR (#437) that surfaced per-stratum calibration and stratification issues a read-review had
missed. Companion edit: `case-by-case` SKILL.md de-tracked and de-ranked (removed SPANDEX/autoresearch track
references, purged nDCG / top-3 from the described output since `ranking-metrics-retired` is active under CHISEL).

#### What the skill enforces

Eleven checks the reviewer subagent must work through on every qualifying PR:

1. Recompute headline AUC / Brier from predictions CSV, match to 6 decimal places
2. Verify fold disjointness empirically on per-row predictions (not from reading fold-assignment code)
3. NaN / out-of-[0,1] / duplicate-key scan on predictions
4. Check every "invariant" the PR prose asserts — is it enforced as an assert or a warn+drop?
5. Decompose aggregate metrics per source / family / fold; report spread, not just mean
6. Reliability diagram is mandatory whenever Brier is claimed (deciles: mean P vs actual rate)
7. Never conflate AUC (discrimination) with Brier (calibration) in findings
8. Direction check on every explanatory hypothesis — wrong direction rejects the hypothesis
9. Stratification semantics vs name — flag catch-all buckets misnamed as named strata
10. Bootstrap output completeness — `bootstrap_samples_used` alongside `requested`
11. Decompose PR deltas against prior baseline into scientific change vs confounder

#### Why now

The CH05 audit (worktree `../audit-ch05`) recomputed the headline AUCs and found they matched exactly, but also
surfaced: BASEL reliability is 30–55 pp off in mid-P buckets while AUC parity holds, per-family AUC spread is
0.73–0.96 while the headline is 0.88, "ICTV family" stratification puts 40% of the panel in catch-all buckets, and
the earlier "BASEL encoded at wrong log_dilution" hypothesis is direction-wrong. None of those were in the PR body or
the self-review. Codifying the protocol as a skill means the next ML PR picks up the same checks automatically rather
than depending on a human reviewer repeating the audit by hand.

#### Scope of changes

- `.agents/skills/review-ml-pr/SKILL.md` — new skill, ~10 KB, checklist + snapshot discipline + related-skills
  cross-reference to `case-by-case`.
- `.agents/skills/case-by-case/SKILL.md` — de-tracked (removed SPANDEX / SX10-13 / autoresearch references), de-ranked
  (removed nDCG / top-3 / MLC as primary outputs), added a Pre-CHISEL note flagging that `compare_predictions.py`
  still computes retired metrics internally and needs a follow-up refactor.
- `lyzortx/orchestration/AGENTS.md` — Self-Review via Subagent section amended to make the skill invocation
  MANDATORY (not advisory) for ML artifact PRs.

#### Follow-up tracked

- `compare_predictions.py` needs an AUC + Brier + reliability-diff rewrite to match the case-by-case SKILL.md
  description. The SKILL.md currently carries a Pre-CHISEL note acknowledging the lag.
