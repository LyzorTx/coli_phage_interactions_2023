# Track CHISEL — Recap

## Goal

Rebase the phage-host lysis prediction pipeline from SPANDEX's pair-level `any_lysis`
training + MLC-graded nDCG/mAP scorecard onto a per-row binary training frame with
concentration as a numeric feature, an AUC+Brier scorecard, and disclose-everything
cross-source (Guelin ↔ BASEL) generalization semantics. In two sentences: SPANDEX
reported within-panel numbers inflated by cv_group leakage and a label-smoothed rollup
that hid BASEL's deployability gap; CHISEL re-trains honestly, names the gap, and moves
one of its two root causes (phage-side TL17 bias) into a panel-independent feature slot.

## Headline outcomes

**2026-04-21 CH13 close-out:** Arm 3 Moriniere per-receptor-class k-mer fractions
migrated to canonical `phage_projection` slot in both the autoresearch Guelin-only
cache and the unified-panel slot root. TL17 BLAST family-presence vectors retired to
`_tl17/` sensitivity fallbacks. CH04 + CH08 re-run under the new canonical; CH05/CH07/CH09
artifacts from CH10/CH11 were already under Arm 3 (via override) and reproduce
bit-for-bit under the post-migration default, so no re-run needed. **CHISEL track
closed** — all canonical numbers below are the final state.

Everything below is on the 369×96 Guelin panel (unified 148-phage panel for cross-
source numbers). Canonical = per-row binary labels, `pair_concentration__log10_pfu_ml`
feature (absolute log₁₀ pfu/ml), all-pairs only (no per-phage blending), NO neat-only
filter, Arm 3 Moriniere per-receptor-fraction `phage_projection` slot (13 dims),
TL17 BLAST families (33 dims) retired to sensitivity fallback.

| Metric | Frame | Number | 95% CI |
|---|---|---|---|
| CH04 Guelin bacteria-axis AUC | pre-filter + Arm 3 (CH13) | **0.8094** | [0.7956, 0.8226] |
| CH04 Guelin bacteria-axis Brier | pre-filter + Arm 3 (CH13) | **0.1749** | [0.1679, 0.1825] |
| CH05 unified bacteria-axis AUC | pre-filter + Arm 3 (CH11) | **0.8079** | [0.7934, 0.8223] |
| CH05 unified phage-axis AUC | pre-filter + Arm 3 (CH11) | **0.8870** | [0.8658, 0.9055] |
| CH05 BASEL bacteria-axis AUC (subset) | pre-filter + Arm 3 (CH11) | **0.7392** | 7.1 pp below Guelin (was 10.2 pp under post-filter TL17) |
| CH05 BASEL phage-axis AUC (subset) | pre-filter + Arm 3 (CH11) | **0.8952** | exceeds Guelin by 0.8 pp (was −1.0 pp) |
| CH07 both-axis AUC | pre-filter + Arm 3 (CH10) | **0.7634** | [0.7581, 0.7689] |
| CH07 both-axis Brier | pre-filter + Arm 3 (CH10) | **0.1902** | [0.1874, 0.1927] |
| CH08 SX12 Δ AUC (phage_moriniere_kmer) | pre-filter + Arm 3 (CH13) | **+0.0058** | [+0.0020, +0.0094], disjoint — lift survives Arm 3 |
| CH08 SX13 Δ AUC (host_omp_kmer) | pre-filter + Arm 3 (CH13) | +0.0002 | [−0.0011, +0.0014], null |
| CH09 Guelin LOOF ECE (bact / phage) | pre-filter + Arm 3 (CH11 refit) | **0.0057 / 0.0052** | target < 0.02 ✓ |
| CH09 BASEL ECE closure (bact / phage) | pre-filter + Arm 3 (CH11 refit) | **64.3% / 44.6%** | residual feature-space deficit, not threshold |

The load-bearing cold-start number (both-axis AUC on simultaneously unseen
bacterium × phage) is **0.7634** [0.7581, 0.7689] under the pre-filter canonical +
Arm 3 slot (CH10). Post-filter CH07 reported 0.7749; the −1.15 pp gap vs the
post-filter number is consistent with the label-shift mechanism documented in the
CH10 entry — the filter trivialised ~12.6 % of pair eval labels 1→0, producing a
label-population-easier AUC. Cross-source: Guelin 0.7654 [0.7600, 0.7710] vs BASEL
0.7193 [0.6822, 0.7580]; per-cell AUC mean 0.7663, std 0.0438, min 0.6378, max
0.8989 across 100 cells (10 bacteria folds × 10 phage folds).

## What changed in the canonical pipeline

Against the SPANDEX starting point:

+ **Training label**: any_lysis (pair-level rollup) → per-row binary, score ∈ {0, 1},
  `score == "n"` dropped as missing (see `label-policy-binary`). Each (bacterium, phage,
  log_dilution, replicate) observation is a training row.
+ **Concentration**: implicit in the rollup → explicit numeric feature
  `pair_concentration__log10_pfu_ml` with Guelin steps {4.7, 6.7, 7.7, 8.7} and BASEL
  constant 9.0 (absolute log₁₀ pfu/ml; conservative lower bound on Maffei 2021/2025
  >10⁹ pfu/ml working titer).
+ **Per-phage blending (AX02)**: dominant SPANDEX architectural gain +2 pp AUC →
  **retired track-wide** (`per-phage-retired-under-chisel`). Not deployable for unseen
  phages, and the SPANDEX +2 pp was partly leakage and partly per-phage-head artifacts
  under per-row training.
+ **Fold hashing**: name-hashed (leaked 45/48 multi-bacterium cv_groups across folds) →
  cv_group-hashed (see `cv-group-leakage-fixed`). All subsequent CHISEL results sit
  downstream of this fix.
+ **Scorecard**: nDCG + mAP + top-k + AUC + Brier → **AUC + Brier only**. Ranking
  metrics are a product-layer concern; a biological model predicts `P(lysis | host, phage,
  concentration)` and downstream code turns calibrated probabilities into rankings
  (see `ranking-metrics-retired`).
+ **Calibration layer**: none → CH09 isotonic calibrator fitted on Guelin training-fold
  predictions, persisted as `ch09_calibrator.pkl`, closes Guelin ECE from 0.13 to
  0.007/0.006 on both axes (see `chisel-unified-kfold-baseline`).
+ **Phage-side feature slot**: Guelin-derived TL17 BLAST projection (zero-vector for
  13/52 BASEL phages) → Moriniere per-receptor k-mer fractions (CH06 Arm 3,
  panel-independent 13-dim, BASEL zero-vec phage-axis +4.36 pp; see
  `moriniere-receptor-fractions-validated`). Canonical migration deferred to the
  follow-up CH13 ticket (currently a side-materialized artifact at
  `.scratch/basel/feature_slots_arm3/`).
+ **Both-axis cold-start**: never previously measured → CH07 reports 0.7749 on 100 cells
  with pair-level bootstrap.

## Dead ends and null arms

What was tested and found not to lift (one-liners for future tracks that might otherwise
re-litigate):

+ **CH06 Arm 1 (OOD shrinkage toward base rate)**: null. Feature-space-level
  out-of-distribution-detector + shrinkage on phage projection features did not rescue
  BASEL calibration on either axis.
+ **CH06 Arm 2 (MMseqs2 pairwise proteome similarity, PCA-32)**: null on BASEL non-zero-
  projection phages (cannibalized RBP-specific signal), null on Guelin. Partially rescued
  the zero-projection subset but at the cost of the non-zero-projection group.
+ **CH06 Arm 4 (tail-protein-restricted TL17 BLAST)**: null. Strict subset of baseline
  TL17 hits — smaller matched region, no new information, no lift on any subset.
+ **SX11 ordinal losses**: not retired in CHISEL, but out of scope (chased MLC-graded
  potency which doesn't exist in the CHISEL label frame). Already null under SPANDEX.
+ **CH09 Arm 2 (cross-panel calibrator transfer)**: Guelin-fitted isotonic applied to
  BASEL closes 79.5% of bacteria-axis ECE and 53.2% of phage-axis ECE, but residual BASEL
  ECE 0.044/0.111 is 6-17× Guelin's calibrated ECE — TL17-bias is a feature-level
  problem, not a threshold-level one.
+ **CH09 Arm 3 (label-threshold sensitivity) — REVERTED by CH10 (2026-04-21).**
  Dropping Guelin neat-only positives gives headline +1.3 pp AUC, −3.2 pp Brier, +0.7 pp
  ECE — but a 2026-04-21 post-hoc decomposition found the filter also flips 4,428 pair
  eval labels 1→0 (12.6% of eval set), because evaluation pulls the label from the
  pair's max-concentration observation and the removed neat-only positive leaves a 0
  replicate standing. On matched (pre-filter) labels the filter-trained model has
  −1.47 pp AUC (regression) but −0.98 pp Brier (genuine improvement survives).
  CH10 rolled the filter back to an opt-in sensitivity analysis after four reviewer
  objections (wrong proxy, testing-completeness bias, BASEL inconsistency, paper
  does not drop) — see the CH10 section in `track_CHISEL.md` and the revised
  `chisel-baseline` unit for the full decomposition + objections.

And what surprised us on re-audit:

+ **CH08 SX12 (Moriniere 815 phage 5-mers, top-100 variance pre-filter)**: **non-null
  across three successive baselines**, lift shrinks as baseline sharpens.
  Post-filter TL17 (CH08 original): +1.16 pp AUC [+0.82, +1.51]. Pre-filter TL17
  (CH12): **+0.72 pp AUC [+0.36, +1.05]** — ~38% of original was label-shift
  artifact from the CH06-followup filter. Pre-filter Arm 3 canonical (CH13):
  **+0.58 pp AUC [+0.20, +0.94]** — another ~19% subsumed by Arm 3's 13-dim
  aggregation of the same kmers. The remaining +0.58 pp is genuinely
  incremental — Arm 3 and the raw 815-kmer slot coexist, and `kmer-receptor-
  expansion-neutral` stays reopened under CHISEL.
+ **CH08 SX13 (host OMP 5546 5-mers, top-100 variance pre-filter)**: **AUC null
  confirmed across all three baselines**. Post-filter TL17: +0.17 pp [+0.03,
  +0.31] (marginal-positive). Pre-filter TL17 (CH12): +0.02 pp [−0.13, +0.17].
  Pre-filter Arm 3 (CH13): +0.02 pp [−0.11, +0.14]. The marginal post-filter
  positive was filter-driven label-shift. One new CH13 wrinkle: Brier Δ goes
  from null under prior baselines to **−0.09 pp [−0.15, −0.02], disjoint below
  zero** under Arm 3 canonical — a tiny calibration improvement as host-OMP
  kmers pick up slack that TL17's 33-dim phage_projection used to cover.
  `host-omp-variation-unpredictive` stays dead-end as a discrimination
  direction; Brier-only signal is below the deployment-materiality threshold.

## Open follow-ups

Each is a concrete, schedulable item — not a vague aspiration. (All CH10-CH13
scope complete; CH13 close-out is in the headline narrative above, not here.)

1. **Slot registry allowlist hoist.** `ch04_parallel.FEATURE_COLUMN_PREFIXES` is a
   hardcoded allowlist; a slot attached without its prefix registered is silently dropped
   from model features (caught during CH08). Move the prefix declaration into the slot
   artifact itself so attaching auto-registers.
2. **CH08 top-K variance ablation.** Re-run SX12 at K ∈ {50, 200, 400} to check whether
   the +0.58 pp CH13 Arm-3-canonical lift is K-stable. If K-sensitive, demote the claim.
3. **SX13 Brier-only signal under Arm 3 — lineage-confound check.** CH13 surfaced a
   Guelin bact-axis ΔBrier −0.09 pp (CI [−0.15, −0.02], disjoint below zero) with AUC
   null. Permute host-OMP kmer values within phylogroup and re-run SX13; if the Brier
   signal survives, it is OMP-specific; if it vanishes, it is lineage-confounded.
4. **Bootstrap `ar02_schema_manifest_v1.json` into `prepare.py`.** CH13's Arm 3 slot
   migration currently lives in `ad_hoc_analysis_code/ch13_arm3_migration/` and rewrites
   the top-level cache schema after the fact. Rewire `prepare.py` to emit Arm 3 as the
   default `phage_projection` slot so fresh clones don't need the two-step bootstrap.
5. **BASEL bacteria-axis deficit (7.1 pp residual under Arm 3 canonical).** CH11
   narrowed the deficit from 10.2 pp (post-filter TL17) to 7.1 pp (pre-filter Arm 3).
   The residual is TL17-bias-free by construction, so the remaining gap is
   panel-mismatch in the host feature encoding or in the 39 non-zero-projection BASEL
   phages. Panel expansion (new Guelin-like spot-test data) is the dominant lever —
   see `panel-size-ceiling`.
6. **Cell-level parallelism for CH07.** 100-cell both-axis CV took ~4 h; a 3-4×
   speedup is available if `fit_seeds` is refactored to run multiple cells in parallel
   with memory budgeting. Not blocking but useful for future cheap both-axis re-runs.
7. **Reliability diagnostic at top P deciles.** Raw CH04 predictions overpredict lysis
   by 18-31 pp at predicted-probability deciles 6-9 (ECE 0.122, max|gap| 0.317) — this
   shape is approximately unchanged across TL17 and Arm 3 baselines and is the reason
   CH09 isotonic is load-bearing before any deployment decision. A focused diagnostic
   (per-bacterium reliability trace, per-Straboviridae trace) would help localise
   whether the upper-P bucket is dominated by a specific failure mode.

## Artifact pointers

Canonical generated-outputs directories (under `lyzortx/generated_outputs/`):

+ `ch02_cv_group_fix/` — CV fold-hashing fix + SX10 revalidation.
+ `ch03_row_expansion/` — per-row training matrix + any_lysis regression check.
+ `ch04_chisel_baseline/ch04_aggregate_metrics.json` — CH13 Arm 3 canonical,
  bacteria-axis AUC 0.8094 [0.7956, 0.8226]. Pre-CH13 TL17 canonical preserved at
  `ch04_chisel_baseline_tl17/`.
+ `ch05_unified_kfold/ch05_combined_summary.json` — unified Guelin+BASEL two-axis
  baseline.
+ `ch06_arm1_ood_shrinkage/`, `ch06_arm2_mmseqs_proteome/`,
  `ch06_arm3_moriniere_receptor/`, `ch06_arm4_tail_restricted_tl17/` — per-arm metrics
  JSON, prediction CSVs, cross-source breakdowns.
+ `ch07_both_axis_holdout/ch07_aggregate.json`, `.../ch07_cell_metrics.csv`,
  `.../ch07_cell_distribution.png`.
+ `ch08_wave2_reaudit/ch08_summary.csv`, `.../ch08_sx12_delta.json`,
  `.../ch08_sx13_delta.json`.
+ `ch09_calibration_layer/ch09_calibrator.pkl`, `.../ch09_calibration_report.json`,
  `.../ch09_label_threshold_sensitivity.json`.

Scripts that reproduce headline numbers:

+ `lyzortx/pipeline/autoresearch/ch04_eval.py` — CH04 pre-filter canonical (CH10);
  pass `--drop-high-titer-only-positives` to reproduce the deprecated post-filter
  numbers for sensitivity analysis.
+ `lyzortx/pipeline/autoresearch/ch05_eval.py` — unified Guelin+BASEL k-fold.
+ `lyzortx/pipeline/autoresearch/ch06_arm3_moriniere_receptor.py` — Moriniere per-
  receptor k-mer-fraction slot materializer + eval driver.
+ `lyzortx/pipeline/autoresearch/ch07_both_axis_holdout.py` — 100-cell both-axis CV.
+ `lyzortx/pipeline/autoresearch/ch08_wave2_reaudit.py` — SX12 + SX13 re-audit with
  top-100 variance pre-filter + paired bacterium-level bootstrap.
+ `lyzortx/pipeline/autoresearch/ch09_calibration_layer.py`,
  `lyzortx/pipeline/autoresearch/ch09_arm3_analysis.py` — isotonic calibrator + label-
  threshold sensitivity.

Knowledge units updated or added: `label-policy-binary`, `cv-group-leakage-fixed`,
`chisel-baseline`, `chisel-unified-kfold-baseline`, `per-phage-retired-under-chisel`,
`ranking-metrics-retired`, `moriniere-receptor-fractions-validated`, and the new
`chisel-both-axis-holdout` (pending, per CH07 PR #450 notebook entry).

Detail on every individual ticket lives in `track_CHISEL.md` — this recap is the digest.
