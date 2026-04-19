---
name: case-by-case
description: >
  Compare two phage–host prediction files case-by-case (per-bacterium AUC + Brier deltas, lysis-rate decile
  stratification, narrow-host sign test, permutation significance, named-case spotlight). Use whenever a ticket
  produces a new predictions CSV to be compared against a baseline, to audit whether aggregate metric deltas
  reflect real signal or just noise. Triggers on "case-by-case", "compare predictions", "per-bacterium delta",
  "audit this result", "is this noise", "case study", or when the user wants to check whether a ticket's metric
  delta comes from real localized wins or a uniform tiny shift.
user-invocable: true
argument-hint: "<baseline_predictions.csv> <candidate_predictions.csv> [--focus BACT1,BACT2] [--out report.md]"
---

# Case-by-Case Prediction Comparison

Whenever a ticket produces a new predictions CSV to compare against a baseline, run this skill to decide whether the
aggregate metric deltas (usually small) reflect real localized signal or are statistical noise. Complements
`review-ml-pr` (which verifies a single artifact) when the question is "did this new run improve over the prior
run, and is that improvement real or noise?"

## When to use

- A ticket produces a new predictions CSV alongside bootstrap CIs
- Aggregate AUC or Brier deltas vs the reference baseline are small (< ~2 pp)
- You need to decide: clean null, sub-threshold signal, or cherry-picked outlier?
- The acceptance criteria include a named case (a specific bacterium, phage, or phage family) and you need to check
  whether its improvement generalizes

## How to invoke

```
python .agents/skills/case-by-case/compare_predictions.py \
    <baseline_predictions.csv> <candidate_predictions.csv> \
    [--label-col label_row_binary] [--entity-col bacteria] [--pred-col predicted_probability] \
    [--focus BAC1,BAC2,BAC3] \
    [--narrow-rate 0.15] [--broad-rate 0.30] \
    [--permutations 200] \
    [--out report.md]
```

Both CSVs must share a `pair_id` column (or an equivalent unique key the two files can be merged on). Override the
column defaults via the flags above when the prediction schema differs. Report prints to stdout and optionally to
`report.md`.

## What it reports

1. **Aggregate metrics comparison** — per-bacterium AUC and Brier under each prediction file (mean / median / spread)
2. **Win / loss / tie counts** — how many bacteria improve, degrade, tie at the configured delta threshold
3. **Lysis-rate decile stratification** — where in the narrow ↔ broad spectrum the delta is concentrated
4. **Sign test on narrow hosts** — is the narrow-host effect statistically significant?
5. **Permutation test on aggregate delta** — swap predictions per pair 50/50 N times; is the observed delta typical of
   random swaps?
6. **Named-case spotlight** — for each `--focus` bacterium: per-pair prediction shifts for all its positives; peer
   comparison (bacteria within ±3 pp of its lysis rate) to flag outlier wins

## Writing up the result

After running, copy the `# Headline` and `# Decision` sections into the track notebook entry. The skill emphasizes
honest readings over optimistic ones:

- A delta that passes the acceptance gate + survives permutation test + has aligned decile pattern = **real signal**
- A delta that fails permutation test but has a clean decile pattern = **sub-threshold directional signal** — note it,
  don't adopt
- A delta that's dominated by one outlier bacterium sitting well above its peer group at the same lysis rate =
  **cherry-picked outlier, not reproducible rescue** — say so explicitly

## Implementation notes

- The permutation test pairs each row's baseline/candidate predictions and swaps them randomly. This is the most
  conservative possible null ("your candidate is the baseline plus noise"); if the observed delta is typical under this
  null, the candidate adds nothing.
- Decile boundaries use `pd.qcut` on the full bacteria set; with ~300 bacteria, each decile has ~30, enough for a
  decile-level mean to be interpretable.
- Narrow-host sign test uses the default 15% lysis-rate cutoff. Override with `--narrow-rate` if the project's
  definition of narrow differs.
- Narrow/broad cutoffs are rough heuristics; always sanity-check by eyeballing the per-decile table before drawing
  conclusions.

### Pre-CHISEL note — script may still reference retired ranking outputs

`compare_predictions.py` was originally authored against the pre-CHISEL scorecard (nDCG, top-3) and still computes
those internally. Under CHISEL only AUC and Brier are on the scorecard (see `ranking-metrics-retired`). Treat any
nDCG / top-3 numbers the script prints as legacy diagnostics — do not cite them in notebook entries or knowledge
units. A follow-up is tracked to convert the script to an AUC + Brier + reliability-diff comparison end-to-end.

## Limitations

- Assumes the two prediction files cover the same (pair_id) set. If sets differ, only the intersection is analyzed; a
  warning is printed.
- The permutation test is about the delta being distinguishable from noise — it does NOT account for multiple
  comparisons across tickets. Multiple-comparison correction is a track-level concern, not a skill-level one.
