---
name: review-ml-pr
description: >
  Verify (not just read) ML PRs that produce predictions or metrics artifacts. Recomputes headline
  numbers from scratch, checks fold disjointness empirically, decomposes aggregates into per-stratum
  metrics, builds reliability diagrams when Brier is claimed, and catches "direction-wrong"
  hypotheses the author offered as explanations. Use whenever a PR publishes a predictions CSV, a
  bootstrap metrics JSON, or a new baseline number. Triggers on "review this PR", "verify these
  results", "audit this run", "self-review after push", or any PR that emits a generated_outputs
  artifact.
user-invocable: true
---

# Review ML PR

Reading the diff, reading the acceptance criteria, and running `pytest` is a **correctness check**,
not a **verification**. The numbers in the PR body can come from honest code and still mislead:
headline AUCs hide per-stratum spread, Brier numbers hide reliability shape, "invariants"
documented in the prose can be warnings in the code. This skill is the verification pass.

Run it on every PR that publishes prediction or metric artifacts, after your own self-review diff
scan, before approving.

## When to use

- The PR writes files under `lyzortx/generated_outputs/` (predictions CSV, metrics JSON)
- The PR body cites AUC / Brier or any bootstrap CI (CHISEL scorecard — see `ranking-metrics-retired`)
- The PR adds or updates a knowledge unit with a headline number
- The reviewer needs to decide approve/reject on a claim that depends on a computation

If the PR cites ranking metrics (top-3, nDCG, mAP) as a primary scorecard, that is itself a review finding:
`ranking-metrics-retired` is active under CHISEL, and ranking belongs to the product layer. Flag it and ask the
author to move to AUC + Brier before verifying the numbers.

Skip this skill for pure infra, docs, or workflow PRs with no artifacts to verify.

## Working copy

If the implementer is actively editing, **snapshot the artifacts before auditing**. Reading a
half-rewritten predictions CSV produces nonsense. Either:

- Create a worktree pinned to the PR HEAD (`git worktree add ../audit-<ticket> <branch>`), or
- `cp` the generated outputs into `.scratch/<ticket>_audit_snapshot/` inside your worktree

Never audit against the live `lyzortx/generated_outputs/` path when a rerun might be in flight.

## Verification checklist

Work through this in order. Later checks depend on earlier ones being clean.

### 1. Recompute headline metrics from the predictions CSV

Load the per-pair predictions, compute the aggregate AUC and Brier yourself, compare to the
reported numbers to ≥6 decimal places. Any delta larger than float rounding means the summary JSON
was produced by different code / different inputs than the predictions CSV, which is a bug to stop
and investigate.

Do the same for every cross-source / cross-axis subset the PR reports. "Reported" numbers and
"recomputable from predictions" numbers must match.

### 2. Fold disjointness — empirical, not from code reading

Check that the holdout unit (bacteria for bacteria-axis, phages for phage-axis, etc.) is **disjoint
across folds** and that **every pair lands in exactly one holdout fold**. Do this on the per-row
predictions (which carry `fold_id`), not by reading the fold-assignment function. A well-written
fold-assignment function can still be called on the wrong input.

### 3. Basic hygiene on predictions

Single pass over the predictions CSV:
- NaN predictions count (expect 0)
- Prediction range (expect `[0, 1]` for probabilities)
- Duplicate pair-level keys (expect 0)
- Source / stratum column NaN count (expect 0)

Any of these non-zero is a bug that invalidates downstream metrics — don't continue until resolved.

### 4. Invariant assertions claimed vs enforced

For every invariant the PR docs assert (e.g. "all 25 BASEL ECOR strains are contained in the 369
Guelin panel", "148 = 96 + 52 with no name overlap", "all N pairs evaluated at the same
concentration"), check both:

- **Does the code raise on violation, or only warn + drop?** If only a warning, and the invariant
  is load-bearing for the reported numbers, the PR is one silent drop away from wrong results.
  Flag it and require a hard assert.
- **Does the run log prove the invariant held?** A warning that never fires and a warning that
  silently fired once look the same in the artifact. If the PR asserts the invariant and the code
  uses a warning fallback, ask for positive evidence (explicit `assert`, explicit count in log).

### 5. Decompose every aggregate into strata the domain supports

An aggregate metric on a mixed panel is a weighted average that hides stratum-specific behavior.
For any AUC or Brier the PR reports, recompute the same metric:

- Per source (if the panel is unified across data sources)
- Per family (or per whatever taxonomic / structural stratum the domain uses)
- Per fold
- Per prediction bucket (for reliability — see next item)

Report the **spread**, not just the mean. A 0.88 aggregate AUC with a 0.73–0.96 per-family range
is a different finding than a 0.88 aggregate with 0.86–0.90 spread. The PR almost never reports the
spread; that's what this step is for.

### 6. Reliability diagram is mandatory whenever Brier is claimed

A single Brier number averages calibration across the probability range. It cannot distinguish
"model well-calibrated everywhere" from "model ~right on 5% and 95% predictions, wildly off in the
middle." Those are clinically different.

For any PR that reports Brier, bucket predictions into deciles of predicted P and report each
bucket's (mean predicted P, actual rate, delta). If any bucket shows |delta| > 0.2 on >50 rows,
the calibration story is substantially more complicated than the Brier number suggests, and the
knowledge unit should name it.

### 7. AUC vs Brier in every finding

Never let one stand in for the other. AUC measures ranking / discrimination. Brier measures
calibration. A model with AUC 0.88 and Brier 0.19 on a cohort is well-ranked but poorly calibrated
— that's a deployment-relevant distinction (ranking-based recommendations fine; threshold-based
decisions broken).

Reject wording like "generalizes as well" that leans on AUC parity while ignoring Brier divergence,
or vice versa.

### 8. Direction check on explanatory hypotheses

For every hypothesis the PR offers to explain an observation ("X is caused by mechanism Y"),
confirm the **direction** of the observed effect matches what Y would predict. Wrong direction
rejects the hypothesis regardless of how nice the narrative sounds.

Example: "BASEL is over-predicted because its concentration is encoded lower than actual" — but
the model learned higher-concentration → more-lysis, so lower encoding should cause *under*-prediction,
not over-prediction. Direction wrong, hypothesis rejected.

### 9. Stratification semantics vs name

When a PR says "stratified by X" (e.g. "ICTV family"), confirm the actual strata are X, not
`X + catch-all buckets`. If "UNKNOWN" or "Other" is a stratum with >10% of the panel, the
stratification is enforcing balance across a bucket, not across true strata. Name it accurately in
the knowledge unit — don't call it "X-stratified" when 40% of the panel is in a pseudo-X bucket.

### 10. Bootstrap output completeness

If the PR reports bootstrap CIs, check that the artifact includes `bootstrap_samples_used`
alongside `bootstrap_samples_requested`. Degenerate resamples (all-positive or all-negative units)
get silently skipped on small cohorts; without the used count, wide CIs and thin CIs look the same.

### 11. Comparison to prior baseline

If the PR supersedes a prior baseline, check the change decomposition:
- How much of the delta is the scientific change the PR ships (new feature, new loss, new split)?
- How much is a confounder (bug fix, fold reshuffling, label change, filter change)?

If the change is bundled, the knowledge unit should decompose it. Bundled "one number moved" claims
are hard to validate and harder to inherit.

## Related skills

- `case-by-case` — when the PR compares two prediction sets, use `case-by-case` for per-bacterium
  delta / permutation / named-case spotlight. `review-ml-pr` focuses on single-artifact
  verification; `case-by-case` focuses on pair-of-artifacts comparison. Use both when applicable.

## How to report findings

Split findings into two lists:

- **Verified clean**: the specific invariants and recomputations that passed. This is the evidence
  the PR is trustworthy on the axes checked. Report as a short table, not prose.
- **Surfaced issues**: numbered, each citing (a) what the PR claims, (b) what the artifact shows,
  (c) what change in the PR / knowledge / code would resolve it.

End with an explicit approve / reject call. Unverified scientific claims cannot be approved — if a
check couldn't be completed (artifact missing, snapshot stale, a subroutine errored), say so and
block on the missing verification rather than defaulting to approve.

## Artifact discipline

Keep audit scripts under the PR's audit worktree's `.scratch/` (gitignored, don't commit). Findings
go into the PR review comment / notebook entry, not into committed files. The audit script is
disposable; the finding is the deliverable.

## Rescoping the self-review subagent

Treat this skill as the verify step in the lyzortx/orchestration/AGENTS.md "Self-Review via
Subagent" protocol. The subagent's priority (2) step (scientific / biological / logical
correctness) should be read as **verify correctness from artifacts**, not **review correctness
from the diff**. Read-review misses exactly the class of issues this skill catches.

## Limitations

- Does not replace per-ticket acceptance criteria. Acceptance criteria define what the PR should
  deliver; this skill checks whether the delivered artifacts back the claimed numbers.
- Cannot detect bias introduced by a feature-engineering choice upstream of the predictions being
  audited. If the feature cache is contaminated, the predictions can be self-consistently wrong —
  that's a separate audit at the feature-pipeline layer.
- The stratification checks depend on the domain knowing what strata matter. Consult the knowledge
  model before picking stratification dimensions; don't invent them.
