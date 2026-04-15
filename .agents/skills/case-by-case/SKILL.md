---
name: case-by-case
description: >
  Compare two phage–host prediction files case-by-case (per-bacterium nDCG deltas, lysis-rate decile stratification,
  narrow-host sign test, permutation significance, named-case spotlight, top-3 hit rate). Use after every SPANDEX /
  autoresearch ticket that produces a new predictions CSV, to audit whether aggregate metric deltas reflect real signal
  or just noise. Triggers on "case-by-case", "compare predictions", "per-bacterium delta", "audit this result",
  "is this noise", "NILS53 analysis", "case study", or when the user wants to check whether a ticket's metric delta
  comes from real localized wins or a uniform tiny shift.
user-invocable: true
argument-hint: "<baseline_predictions.csv> <candidate_predictions.csv> [--focus BACT1,BACT2] [--out report.md]"
---

# Case-by-Case Prediction Comparison

After every autoresearch or SPANDEX ticket that produces a new predictions CSV, run this skill to decide whether the
aggregate metric deltas (usually small) reflect real localized signal or are statistical noise. The same diagnostic
protocol was used in SX11, SX12, and SX13 — bake it into a reusable tool so we don't rebuild it every time.

## When to use

- A ticket produces a new predictions CSV alongside bootstrap CIs
- Aggregate AUC/nDCG deltas vs the reference baseline are small (< ~2 pp)
- You need to decide: clean null, sub-threshold signal, or cherry-picked outlier?
- The acceptance criteria include a named case (NILS53, Dhillonvirus narrow-host, ECOR-XX, etc.) and you need to check
  whether its improvement generalizes

## How to invoke

```
python .agents/skills/case-by-case/compare_predictions.py \
    <baseline_predictions.csv> <candidate_predictions.csv> \
    [--label-col mlc_score] [--entity-col bacteria] [--pred-col predicted_probability] \
    [--focus NILS53,ECOR-19,EDL933] \
    [--narrow-rate 0.15] [--broad-rate 0.30] \
    [--permutations 200] \
    [--out report.md]
```

Both CSVs must share a `pair_id` column (or an equivalent unique key the two files can be merged on). The default
column names match SX10/SX11/SX12/SX13 prediction outputs. Report prints to stdout and optionally to `report.md`.

## What it reports

1. **Aggregate metrics comparison** — per-bacterium nDCG under each prediction file (mean / median / spread)
2. **Win / loss / tie counts** — how many bacteria improve, degrade, tie at ≥1 pp threshold
3. **Lysis-rate decile stratification** — where in the narrow ↔ broad spectrum the delta is concentrated
4. **Sign test on narrow hosts** — is the narrow-host effect statistically significant?
5. **Permutation test on aggregate delta** — swap predictions per pair 50/50 N times; is the observed delta typical of
   random swaps?
6. **Named-case spotlight** — for each `--focus` bacterium: per-pair rank changes for all its positives; peer
   comparison (bacteria within ±3 pp of its lysis rate) to flag outlier wins
7. **Top-3 hit rate** — number of bacteria with at least one top-3 prediction that is a true positive

## Writing up the result

After running, copy the `# Headline` and `# Decision` sections into the track notebook entry. The skill emphasizes
honest readings over optimistic ones:

- A delta that passes the acceptance gate + survives permutation test + has aligned decile pattern = **real signal**
- A delta that fails permutation test but has a clean decile pattern = **sub-threshold directional signal** — note it,
  don't adopt
- A delta that's dominated by one outlier bacterium (e.g., NILS53 at 76th percentile of its peers) = **cherry-picked
  outlier, not reproducible rescue** — say so explicitly

## Implementation notes

- The permutation test pairs each row's baseline/candidate predictions and swaps them randomly. This is the most
  conservative possible null ("your candidate is the baseline plus noise"); if the observed delta is typical under this
  null, the candidate adds nothing.
- Decile boundaries use `pd.qcut` on the full bacteria set; with ~300 bacteria, each decile has ~30, enough for a
  decile-level mean to be interpretable.
- Narrow-host sign test uses the default 15% lysis-rate cutoff. Override with `--narrow-rate` if the project's
  definition of narrow differs.
- The script treats `mlc_score` as the relevance for nDCG and `mlc_score > 0` as the binary positive label. Override
  with `--label-col` / `--binary-label-col` if your prediction files use different schemas.
- For runtime, 200 permutations × 356 bacteria × nDCG computation = ~30 seconds. Don't crank permutations past 500
  unless you have a specific reason — the null distribution stabilizes quickly.

## Limitations

- Assumes the two prediction files cover the same (pair_id) set. If sets differ, only the intersection is analyzed; a
  warning is printed.
- nDCG comparison ignores bacteria with no positives and bacteria with only 1 pair (sklearn restriction).
- The permutation test is about the delta being distinguishable from noise — it does NOT account for multiple
  comparisons across tickets. That's a track-level concern (SX14 stratified eval).
- Narrow/broad cutoffs are rough heuristics; always sanity-check by eyeballing the per-decile table before drawing
  conclusions.
