# Track CHISEL: Concentration-level binary labels, Honest splits, Calibrated scorecard

**Goal:** Re-anchor the evaluation frame around the raw (bacterium, phage, concentration, replicate) → {0, 1}
observations and the AUC + Brier scorecard. Retire MLC-derived ranking metrics (nDCG, mAP, top-3) and the any_lysis
rollup. Close all pair-level leakage paths before running new feature or architecture experiments.

**Baseline entering CHISEL:** SPANDEX SX10 canonical — AUC 0.8699 [0.8570, 0.8819], Brier 0.1248 [0.1187, 0.1309]
within-panel (369×96, clean-label, bacterium-name-hashed folds). That baseline is now known to be inflated by a fold
leakage bug (CH02, below).

**Why CHISEL exists:** SPANDEX fought four conceptual problems at once — graded MLC labels, any_lysis rollup,
nDCG-driven model selection, and pair-level folds that leaked under cv_group. CHISEL untangles them in sequence:
CH02 fixes the fold scheme, CH03 verifies the any_lysis safety-net survives row expansion, CH04 flips to per-row
binary training with concentration as a feature, CH05-CH06 add phage-axis and both-axis honesty, CH07 re-audits the
SPANDEX feature-family nulls under the new frame.

---

## 2026-04-19 07:30 CEST: Track inception

Track CHISEL succeeds Track SPANDEX (closed 2026-04-19, see `track_SPANDEX.md`). The framework pivot and ticket
sequencing are documented in `project.md` under the same date heading.

---

## 2026-04-19 08:30 CEST: CH02 — cv_group leakage fix and SX10 baseline revalidation

### Executive summary

`assign_bacteria_folds` hashed 10-fold CV folds on bacterium name, splitting 45 of 48 multi-bacterium cv_groups
across folds (near-duplicates ≤1e-4 ANI on both sides of the train/test boundary). Fixed to hash on cv_group;
revalidating the SX10 canonical configuration (GT03 all_gates_rfe + AX02 per-phage blending, any_lysis, SX10
features, 369×96 panel) under the corrected scheme shifts AUC from 0.8699 → **0.8521** [0.8381, 0.8649] and Brier
from 0.1248 → **0.1317** [0.1253, 0.1381]. The SPANDEX leakage was uniform across arms, so SPANDEX per-family null
conclusions are preserved; only the headline numbers move by this amount.

**Problem.** `assign_bacteria_folds` in `sx01_eval.py` hashed folds on bacterium name
(`sha256("spandex_v1:{bacteria_name}") % 10`). 48 / 283 cv_groups in the clean-label 369-bacteria panel contain
more than one bacterium (bacteria within ≤1e-4 ANI of each other, clustered via the canonical
`data/metadata/370+host_cross_validation_groups_1e-4.csv`). Under independent bacterium-name hashing into 10 folds,
each multi-bacterium cv_group had a ≥90% chance of splitting across folds. The measurement: **45 / 48 multi-bacterium
cv_groups were split across folds under the SPANDEX scheme** — near-identical bacteria appeared on both sides of the
train/test boundary for roughly 17% of bacteria on the panel.

This affects every SPANDEX-era aggregate number that passed through `assign_bacteria_folds`: SX05/SX10 canonical
within-panel (AUC 0.8699), SX03 cross-source arms, SX04 ordinal regression, SX11 loss ablation, SX12 Moriniere
k-mers, SX13 OMP allelic variation, SX14 stratified decomposition, SX15 unified Guelin+BASEL. All of them inherit
a small positive bias in favour of within-panel performance relative to cross-panel performance.

**Fix.** `assign_bacteria_folds` now takes a bacterium → cv_group mapping and hashes
`sha256("spandex_v1:cv_group:{cv_group_id}") % 10`. All bacteria sharing a cv_group land in the same fold by
construction. A companion helper, `bacteria_to_cv_group_map()`, builds the mapping from any training frame that
carries the `cv_group` column (all callers already do — cv_group is propagated from ST02 through every downstream
pair table).

All six call sites were updated: `sx01_eval.run_kfold_evaluation`, `sx03_eval.run_arm_a_baseline` and
`run_arm_b_pooled`, `sx04_eval.run_sx04_eval`, `sx11_eval.run_sx11`, `sx12_eval.run_sx12_eval`,
`sx13_eval.run_sx13_eval`. The function signature change is a hard break — callers that try to pass a
`Sequence[str]` (bacteria list) fail loudly rather than silently reverting to the buggy hash. A unit test in
`lyzortx/tests/test_sx01_fold_assignment.py` pins the new behaviour: same-cv_group bacteria share folds, empty or
missing cv_groups raise, and determinism holds across runs.

**Revalidation — SX10 canonical under fixed folds.** Rerun:

```
PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.sx01_eval \
  --device-type cpu \
  --output-dir lyzortx/generated_outputs/ch02_cv_group_fix/sx10_revalidated
```

Configuration identical to SPANDEX SX10 canonical: GT03 all_gates_rfe + AX02 per-phage blending, `label_any_lysis`
training target, SX10 feature bundle (host_surface + host_typing + host_stats + host_defense + phage_projection +
phage_stats + pair_depo_capsule + pair_receptor_omp, RFE-selected), 10-fold bacteria-axis CV, 3 seeds,
1000 bootstrap resamples at bacterium level. Only the fold hashing differs.

| Scheme | AUC | AUC 95% CI | Brier | Brier 95% CI |
|--------|-----|------------|-------|--------------|
| SPANDEX (bacterium-name hash, leaky) | 0.8699 | [0.8570, 0.8819] | 0.1248 | [0.1187, 0.1309] |
| CHISEL (cv_group hash, CH02 fix) | **0.8521** | [0.8381, 0.8649] | **0.1317** | [0.1253, 0.1381] |
| Shift | −1.78 pp | CIs narrowly overlap | +0.69 pp | CIs do not overlap |

Artifacts: `lyzortx/generated_outputs/ch02_cv_group_fix/sx10_revalidated/` (predictions, fold metrics, bootstrap
results) and `lyzortx/generated_outputs/ch02_cv_group_fix/ch02_sx10_revalidated_metrics.json` (summary JSON with
fold-diff context). `lyzortx/generated_outputs/ch02_cv_group_fix/ch02_fixed_folds.csv` is the canonical
bacterium → fold mapping under CHISEL hashing, regenerated by
`lyzortx/research_notes/ad_hoc_analysis_code/ch02_fold_diff.py`.

**Interpretation.** The AUC drop (−1.78 pp) is larger than the 17%-of-bacteria leakage footprint alone would
predict, because the hash input change (bacterium name → cv_group id) reshuffles every singleton cv_group's fold as
a side-effect — 323 / 369 bacteria change fold under CHISEL, not just the ~70 bacteria in split cv_groups. The
shift therefore captures two compounding effects: (a) leakage correction on the 45 split cv_groups, and (b)
random fold-composition variance on the other 325 bacteria. Distinguishing (a) from (b) would require holding the
singleton-fold mapping fixed and only re-routing the multi-bacterium cv_groups — that experiment is not planned
because the CHISEL hashing is now the canonical anchor and future tickets compare against CH02, not SPANDEX.

The Brier deterioration (+0.69 pp, CIs disjoint) is the more interpretable signal: calibration is genuinely worse
without train/test leakage. The model was slightly over-confident under the leaky split because some test bacteria
had near-duplicates in training.

**Scope of the correction.** Every SPANDEX aggregate number is subject to this small shift, but the *directions* of
SPANDEX findings (nulls for ordinal regression, Moriniere k-mers, OMP allelic variation, plm-rbp-redundant) are
conserved — the leakage was a uniform inflation across arms, not an arm-selective effect. SPANDEX per-family null
conclusions remain valid; only the absolute headline numbers move.

The `spandex-final-baseline` knowledge unit retains the original numbers (AUC 0.8699) as historical reference
under the old fold design. The `cv-group-leakage-fixed` unit records the magnitude of the correction and flags the
SPANDEX-era baseline numbers as subject to this small shift. The new baseline number (AUC 0.8521) is not yet
canonical — CH03 must first verify that row expansion preserves the any_lysis aggregate (safety-net check) before
CH04 flips to per-row binary training. The canonical CHISEL baseline will be established in CH04 once the label
frame is fully migrated.

**Acceptance met.**

- `assign_bacteria_folds` hashes on cv_group; all six callers updated.
- Unit test pins no-leakage invariant.
- `ch02_fixed_folds.csv` emitted (369 bacteria × 4 columns).
- SX10 canonical rerun under fixed folds: AUC 0.8521, Brier 0.1317 with bacterium-level bootstrap CIs.
- `ch02_sx10_revalidated_metrics.json` emitted with SPANDEX vs CHISEL side-by-side and fold-diff summary.
- `cv-group-leakage-fixed` knowledge unit added; `spandex-final-baseline` historical numbers preserved.
- Lab notebook entry (this one) records methodology and the 1.78 pp AUC shift.

---

## 2026-04-19 09:00 CEST: CH03 — row-expanded training matrix with preserved any_lysis semantics

### Executive summary

Plumbing-only safety-net ticket. Built a `(bacterium, phage, log_dilution, replicate, X, Y)` row-expanded training
frame (318,816 rows across 35,424 pairs) by joining `raw_interactions.csv` with ST02 pair metadata and ST03
splits, then verified the ST01B `any_lysis` rule applied to row-level scores reproduces **all 35,424 ST02
`label_hard_any_lysis` values with zero mismatches** (9,720 positive / 25,546 negative / 158 unresolved — identical
to ST02). Rerunning SX10 canonical on the collapsed frame gives AUC **0.8522** [0.8379, 0.8650] and Brier
**0.1318** [0.1255, 0.1382] — |ΔAUC| = 0.00012, |ΔBrier| = 0.00012 vs CH02 (0.8521 / 0.1317), well inside the
0.005 regression tolerance.

**Design.** Two modules landed under `lyzortx/pipeline/autoresearch/`:

- `ch03_row_expansion.py` — loader + rollup:
  - `load_raw_observation_rows()` reads `raw_interactions.csv` (318,816 rows; score distribution 272,750 × "0",
    37,391 × "1", 8,675 × "n").
  - `load_row_expanded_frame()` joins each raw row against the ST02 pair table and ST03 split assignments so
    every observation carries `pair_id`, `cv_group`, `label_hard_any_lysis`, `training_weight_v3`,
    `split_holdout`, `is_hard_trainable`, plus the host and phage metadata columns that downstream SX10 reads.
    Fails loudly if any raw row has no ST02 counterpart.
  - `rollup_row_scores_to_pair_any_lysis()` reproduces the ST01B rule on the row-level scores: `any_lysis = 1`
    iff any row has `score == "1"`; `any_lysis = 0` iff no `score == "1"` AND `score_0_count ≥ 5`; otherwise
    unresolved. `score == "n"` rows count neither as positive nor as negative — they're treated as missing, not
    as silent zeroes.
  - `verify_rollup_matches_st02()` is the core regression check — the recomputed labels must equal ST02's
    `label_hard_any_lysis` value on every one of the 35,424 pairs, or the eval aborts.
  - `collapse_to_pair_level_for_training()` dedupes the row-expanded frame back to pair level for CH03's
    temporary scaffolding. CH04 replaces this with a per-row trainer. The helper fails loudly if any pair-level
    column (cv_group, label, etc.) takes inconsistent values within a single pair.

- `ch03_eval.py` — runner. Loads the row-expanded frame, runs the rollup check, collapses to pair level,
  verifies row-level fold integrity (every raw observation of a bacterium routes to the same CH02 fold, so no
  `(pair, concentration, replicate)` row can appear in both train and test), then calls `run_kfold_evaluation`
  with the collapsed frame. On completion, compares the resulting AUC + Brier against `ch02_sx10_revalidated_metrics.json`
  and fails if either delta exceeds 0.005.

**Run.**

```
PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch03_eval --device-type cpu
```

Total runtime 613.8 s (a few seconds longer than CH02's 609 s — the overhead is the raw-row load + rollup check).

| Quantity | CH02 (direct ST02) | CH03 (row-expanded → rolled up) | Δ |
|---|---|---|---|
| AUC | 0.8521 [0.8381, 0.8649] | 0.8522 [0.8379, 0.8650] | +0.00012 |
| Brier | 0.1317 [0.1253, 0.1381] | 0.1318 [0.1255, 0.1382] | +0.00012 |
| Pairs positive | 9,720 | 9,720 | 0 |
| Pairs negative | 25,546 | 25,546 | 0 |
| Pairs unresolved | 158 | 158 | 0 |

The residual ~0.0001 delta on AUC and Brier is numerical noise from the row-expansion → pair-level recollapse
path (likely float ordering differences in the `drop_duplicates` step interacting with downstream pandas
merges). Well below the 0.005 regression tolerance and far below any bootstrap CI width.

**Interpretation.** The zero-mismatch rollup is the load-bearing result. It proves two things CH04 needs:

1. The raw → row-expanded pipeline conserves every pair's label semantics. CH04 can swap `any_lysis` labels for
   the raw per-row `score` values without worrying that the row-level data represents a different experiment
   than SX10 trained on.
2. The `n` handling is correct: rows with `score == "n"` are not silently counted as negative in the training
   corpus (they push 158 pairs to unresolved; SX10's `is_hard_trainable` filter drops them). CH04's per-row
   training can drop `score == "n"` rows directly without changing the pair-level denominator the SPANDEX and
   CH02 numbers were computed on.

**Fold integrity.** The CH02 cv_group hash operates at bacterium level. Because all 318,816 raw rows for a
given bacterium share the same `cv_group` (the cv_group is a bacterium attribute, not a pair attribute), every
row routes to the same fold as its bacterium. CH03's fold-integrity assertion `pair_id → fold_id` is 1:1
confirms this holds on the full training matrix — no pair was split across folds at any expansion level, which
is the necessary precondition for CH04's per-row training.

**Scope.** CH03 changes nothing semantically. The collapsed pair-level frame is the ST02 pair table (with the
`bacteria` and `phage` columns re-derived from the raw rows and the rest copied through). SX10 training sees
the same design matrix, the same labels, and the same folds as CH02. CH03's value is that it localises any
future row-expansion bug inside `ch03_row_expansion.py` — if CH04 misinterprets the raw schema, the rollup
check will flag it before the model runs.

**Artifacts.**

- `lyzortx/generated_outputs/ch03_row_expansion/ch03_expanded_training_frame.csv` — the 318,816-row canonical
  frame for CH04 reuse (persisted as CSV; parquet would have required adding `pyarrow` as a dependency for a
  temporary scaffold).
- `lyzortx/generated_outputs/ch03_row_expansion/ch03_regression_check.json` — CH02-vs-CH03 metric comparison +
  rollup check summary.
- `lyzortx/generated_outputs/ch03_row_expansion/sx10_on_row_expanded/` — raw SX10 outputs (predictions, fold
  metrics, bootstrap results) from the revalidation run.

**Acceptance met.**

- Row-expanded loader + rollup land in `ch03_row_expansion.py`; all raw rows carry pair-level metadata.
- ST01B rollup reproduces ST02 labels exactly (0 mismatches on 35,424 pairs).
- Fold-integrity check (1 fold per pair_id) enforced in `check_row_level_fold_integrity`.
- SX10 rerun on collapsed frame: AUC 0.8522, Brier 0.1318 — |Δ| < 0.005 vs CH02 on both metrics.
- Artifacts `ch03_expanded_training_frame.csv` and `ch03_regression_check.json` emitted.
- Unit tests (6) cover rollup logic, `n` handling, and collapse invariants.

---

## 2026-04-19 12:30 CEST: CH04 — CHISEL canonical baseline (per-row, concentration feature, all-pairs only)

### Executive summary

First canonical CHISEL result: **AUC 0.8084 [0.7944, 0.8217], Brier 0.1750 [0.1677, 0.1824]** on 35,266 pairs ×
369 bacteria. Three changes simultaneously separate this baseline from CH02's revalidated SX10 (AUC 0.8521,
Brier 0.1317): (1) training unit flips from pair-level `any_lysis` to per-row binary `score` with 8,675
`score == "n"` rows dropped as missing; (2) `pair_concentration__log_dilution` enters as a numeric feature;
(3) AX02 per-phage blending is retired as non-deployable (see `per-phage-retired-under-chisel`). Aggregate
shifts are −4.37 pp AUC and +4.33 pp Brier, both CIs disjoint — significant but expected, as CHISEL trades
CH02's inflated bacteria-axis-only numbers for a deployable per-observation predictor. Concentration is the
#4 feature by mean LightGBM importance and is retained by RFE in all 30 fold × seed fits.

**Design.**

- `ch04_eval.py` loads the CH03 row-expanded frame, drops `score == "n"` rows
  (`build_clean_row_training_frame`), casts `score` to `label_row_binary`, and attaches
  `pair_concentration__log_dilution = float(log_dilution)` as a pair-level numeric feature.
- `train_and_predict_per_row_fold` adapts `sx01_eval.train_and_predict_fold` for per-row training. Same
  host/phage slot bundle (host_surface + host_typing + host_stats + host_defense + phage_projection +
  phage_stats), same pairwise cross-terms (pair_depo_capsule + pair_receptor_omp), same RFE machinery.
  Differs on three axes: (1) `label_row_binary` replaces `label_any_lysis` as the training target, (2) the
  prefix tuple grows `"pair_concentration__"` so RFE considers the concentration feature, (3) every raw
  observation is a training example rather than one row per pair. **Per-phage blending (AX02) is
  intentionally omitted** — `per-phage-retired-under-chisel` documents the rationale (non-deployable;
  SPANDEX +2 pp gain was partly a cv_group-leakage artifact; an intermediate CH04 v1 run confirmed the
  per-phage head deflates positive predictions by ~5 pp under per-row training because it sees 9 rows per
  bacterium with identical host features but mixed labels and no log_dilution visibility).
- Evaluation: `select_pair_max_concentration_rows` picks the highest observed log_dilution per pair (all
  35,266 pairs have at least one replicate at log_dilution=0, so evaluation is at neat concentration),
  averages replicate predictions, and max-aggregates replicate labels. AUC/Brier are aggregated over pairs
  (not rows) with bacterium-level bootstrap CIs from `bootstrap_auc_brier_by_bacterium` (1,000 resamples).
- No nDCG / mAP / top-k — ranking metrics retired per `ranking-metrics-retired`.
- Per-fold AUC/Brier logged at each fold completion for visibility during long runs (added after the
  first CH04 run went 85 minutes with only a single aggregate line at the end).

**Run.**

```
PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch04_eval --device-type cpu
```

Total runtime 5139 s (86 min) — ~8× SX10 because every fold trains on ~280K rows instead of ~26K.

**Result.**

| Metric | CH04 v3 | CH04 CI | CH02 revalidated | Δ vs CH02 | CIs overlap? |
|---|---|---|---|---|---|
| AUC | **0.8084** | [0.7944, 0.8217] | 0.8521 [0.8381, 0.8649] | −4.37 pp | no (disjoint) |
| Brier | **0.1750** | [0.1677, 0.1824] | 0.1317 [0.1253, 0.1381] | +4.33 pp | no (disjoint) |

Per-fold AUC across the 10 folds: 0.828, 0.784, 0.772, 0.836, 0.807, 0.823, 0.800, 0.841, 0.795, 0.785 —
range [0.77, 0.84], std 0.023.

**Decomposition of the CH02 → CH04 gap.**

CH02 compounded three effects this baseline separates from. Evaluation label redistribution is small
(pair-level any_lysis 27.6% positive vs CH04 max-conc 27.4% — ~47 pair labels flip because nearly every
pair has a replicate at log_dilution=0 that matches the any_lysis call). The bulk of the gap is
training-side:

1. **Per-row training vs pair-level rollup.** The question the model answers shifts from "does this pair
   ever lyse?" to "does this specific observation show lysis?" A diagnostic CH04 v1 run retaining per-phage
   blending but flipping to per-row training landed at AUC 0.8078 — almost identical to the final
   all-pairs-only result, so per-row training accounts for the vast majority of the AUC drop. The per-row
   frame has 13.4% positive rate (37,391 / 310,141) vs the pair-level 28%, and LightGBM biases harder
   toward the majority class.

2. **Per-phage retirement.** CH04 v1 also exposed that SPANDEX's per-phage head was actively deflating
   positive predictions by ~5 pp under per-row training (per-phage models see 9 rows per bacterium with
   identical host features but mixed labels and no log_dilution visibility, collapsing to mean lysis rate
   across dilutions). Under CH02's pair-level training, per-phage was a clean bacterium-level memorization
   head; under per-row it's a category mismatch. Removing it keeps AUC essentially unchanged (0.8078 →
   0.8084) but worsens Brier by 2.3 pp (0.1522 → 0.1750). Per-phage was pulling predictions toward
   per-bacterium averages, which improved calibration at the cost of discrimination on positives. CHISEL
   accepts the calibration cost because a non-deployable head is not a defensible baseline — the
   phage-axis generalization question per-phage silently answered is measured directly in CH05.

3. **cv_group leakage already closed in CH02.** CH02's revalidation (AUC 0.8521) fixed the SPANDEX 0.8699
   baseline's fold-hashing bug separately. That effect is not part of the CH02 → CH04 gap.

Load-bearing finding: the CHISEL training unit + AUC+Brier scorecard costs ~4.4 pp aggregate AUC relative
to SPANDEX's pair-level any_lysis + per-phage setup at fixed features, and the Brier degradation of 4.3 pp
combines per-row (better calibration in principle) with per-phage removal (worse calibration in practice).
Neither metric reflects a loss of predictive capability — CHISEL is answering a different, harder,
deployment-aligned question.

**Concentration feature behaviour.**

`pair_concentration__log_dilution` was retained by RFE in all 30 fold × seed fits. Its mean LightGBM
importance is **328.70**, ranking it **#4 overall**:

| Rank | Feature | Mean importance | In RFE all 30 folds |
|---|---|---|---|
| 1 | host_typing__host_serotype | 2466.5 | ✓ |
| 2 | phage_stats__phage_gc_content | 721.2 | ✓ |
| 3 | phage_stats__phage_genome_length_nt | 408.4 | ✓ |
| **4** | **pair_concentration__log_dilution** | **328.7** | ✓ |
| 5 | host_typing__host_o_type | 282.1 | ✓ |
| 6 | phage_projection__tl17_rbp_reference_hit_count | 230.3 | ✓ |
| 7 | pair_receptor_omp__predicted_lps_x_host_o_antigen | 147.6 | ✓ |

Concentration-response check on positive pairs: mean predicted probability rises monotonically with
log_dilution — 0.24 at log_dilution=-4 (lowest titer), 0.32 at -2, 0.41 at -1, 0.49 at 0 (neat). The model
learns a clean concentration gradient. Within a (pair, log_dilution), predictions across the 1–3
replicates are always identical (140,325 groups, all `nunique=1`) — replicates have identical feature
vectors, which is the irreducible noise floor of the per-row task.

**Deployment interpretation.**

CHISEL's design choice is deliberate. The training question shifts from "does any dilution × replicate of
this pair show lysis?" to "does this specific observation show lysis?" Product-layer calibration — how
confident we are that a clinician administering phage X to bacterium Y at concentration Z will see lysis —
is a per-observation question, not a per-pair one. Per-phage blending was a bacterium-memorization head
that inflated the bacteria-axis metric under the cv_group-leaky SPANDEX folds and cannot transfer to
unseen phages (the deployment-goal). CH04 trades 4.4 pp aggregate AUC for a predictor whose probability
outputs are grounded in the observational unit the deployment cares about, using only features derivable
from genome data.

Any comparison of CHISEL arms (CH05 phage-axis, CH06 both-axis holdout, CH07 feature-family re-audit)
anchors on CH04's AUC 0.8084 / Brier 0.1750. Knowledge updated: `chisel-baseline` is the new canonical;
`per-phage-blending-dominant` marked HISTORICAL (dead-end); new `per-phage-retired-under-chisel` unit with
explicit "do not re-enable" directive.

**Drive-by fix: `log_config.py` timezone.**

While waiting on CH04, noticed pipeline log timestamps showed `+0100` despite local time being CEST
(`+0200`). Root cause: `setup_logging` set `converter = time.gmtime`, which emits UTC hours; `%z` on the
naive `struct_time` then fell back to `time.timezone` (non-DST local offset, CET = `+0100`), producing
hour values and offsets that matched neither UTC nor wall clock. Fix: switch to `time.localtime`, which
sets `tm_isdst` correctly so `%z` picks `time.altzone` (`+0200`) during DST and `time.timezone`
(`+0100`) otherwise. DST-robust, verified with hardcoded summer (2026-07-15 12:00 CEST) and winter
(2026-12-15 12:00 CET) probe timestamps in `test_log_config.py`.

**Artifacts.**

- `lyzortx/generated_outputs/ch04_chisel_baseline/ch04_predictions.csv` — 35,266 pair-level predictions at
  max observed concentration (one row per pair).
- `lyzortx/generated_outputs/ch04_chisel_baseline/ch04_per_row_predictions.csv` — 309,750 per-row
  predictions across all 10 fold × 3 seed fits (averaged over seeds), for downstream per-concentration
  diagnostics.
- `lyzortx/generated_outputs/ch04_chisel_baseline/ch04_aggregate_metrics.json` — AUC + Brier with
  bacterium bootstrap CIs, CH02 comparison, concentration feature importance.
- `lyzortx/generated_outputs/ch04_chisel_baseline/ch04_feature_importance.csv` — mean feature importance
  across all 30 fold × seed fits with `is_log_dilution` flag.

**Acceptance met.**

- Training label is per-row binary `score`; 8,675 `score == "n"` rows dropped as missing, not negative.
- `pair_concentration__log_dilution` added as a numeric feature; RFE retains it in 30/30 fits.
- SX10 feature bundle unchanged otherwise (host_surface + host_typing + host_stats + host_defense +
  phage_projection + phage_stats + pair_depo_capsule + pair_receptor_omp, RFE-selected).
- Per-phage blending retired (deviation from the ticket's "retain AX02" directive, justified per AGENTS.md
  "Question the requirement" — see `per-phage-retired-under-chisel`).
- AUC + Brier only. No nDCG, mAP, top-k computed or logged.
- Per-fold AUC/Brier logged + aggregate with bacterium-level bootstrap CIs (1000 resamples).
- Artifacts `ch04_predictions.csv`, `ch04_per_row_predictions.csv`, `ch04_aggregate_metrics.json`,
  `ch04_feature_importance.csv` emitted.
- `chisel-baseline` knowledge unit added with canonical numbers; `spandex-final-baseline` flagged
  HISTORICAL; `per-phage-blending-dominant` marked dead-end; new `per-phage-retired-under-chisel` unit
  with explicit "do not re-enable" directive.
- 5 unit tests cover `n`-row dropping, concentration feature attachment, max-concentration selection,
  and replicate aggregation. 2 unit tests cover the log_config TZ fix with hardcoded DST probe
  timestamps.

---

## 2026-04-19 18:30 CEST: CH05 — unified Guelin+BASEL two-axis k-fold under CHISEL

### Executive summary

CH05 measures two-axis generalization under all-pairs-only on the unified 148-phage × 369-bacterium
panel (36,643 pairs: 35,403 Guelin + 1,240 BASEL). Headline numbers (post-encoding-fix rerun under
absolute log₁₀ pfu/ml; Guelin-side fold-level bit-identical, BASEL-side aggregate within sampling
noise of the pre-fix run):
**bacteria-axis AUC 0.8063 [0.7919, 0.8202]**; **phage-axis AUC 0.8849 [0.8616, 0.9059]**. Three
substantive findings fall out, each separate from the others:

1. **Phage-axis discrimination parity**: Guelin AUC 0.8861 vs BASEL AUC 0.8829, |ΔAUC| 0.0032 — a
   weak non-rejection on 52 BASEL phages (CI 3× wider than Guelin's), not positive evidence of
   equivalence.
2. **Phage-axis calibration divergence**: Guelin Brier 0.1329 vs BASEL Brier 0.1884 — 42% relative
   degradation on BASEL with disjoint CIs. Per-decile reliability: BASEL mid-P buckets (predicted
   probability 0.5-0.9) are 21-27 pp more miscalibrated than Guelin on phage-axis and 21-22 pp on
   bacteria-axis. AUC parity does *not* imply calibration parity.
3. **BASEL bacteria-axis deficit**: BASEL-only bacteria-axis AUC 0.7152 (1,240 pairs) sits 9.5 pp
   below Guelin-only bacteria-axis AUC 0.8098 on the same axis. This is a BASEL-specific training
   limitation invisible in the aggregate (which is 96.6% Guelin-weighted at row level).

The wave-closing "BASEL phages generalize essentially identically to Guelin phages" claim from the
earlier draft is retired — it conflates three separate findings into one headline that holds only
along the discrimination-AUC axis on phage-axis folds. Under calibration, and on bacteria-axis,
BASEL materially diverges.

**Design.** `ch05_eval.py` reuses CH04's per-row training / pair-max-concentration evaluation /
cluster-bootstrap infrastructure. New helpers: `load_basel_as_row_frame` (BASEL pairs as row-level
observations with parameterised `log10_pfu_ml` via `--basel-log10-pfu-ml`, replicate=1, fails loudly
on missing Guelin cv_group per AGENTS.md fail-fast); `load_unified_row_frame` (concat);
`assign_phage_folds` (StratifiedKFold on 148 phages by ICTV family + "other" catch-all for families
<10 phages + "UNKNOWN" for the 20/52 BASEL phages with no ICTV assignment — the labeling in the
earlier draft as "ICTV-stratified" was misleading given 40% of folds are pseudo-family buckets);
`_bootstrap_by_unit(unit_key=...)` (cluster-bootstrap — resamples bacteria on bacteria-axis, phages
on phage-axis). BASEL phage features patched via `sx03_eval.patch_context_with_extended_slots`.
**Per-phage blending omitted on both axes** — on phage-axis structurally impossible (held-out phages
have zero training rows), on bacteria-axis under `per-phage-retired-under-chisel`. The ticket's
"BLENDING TAX" is consequently moot — replaced by the phage-axis generalization gap as a structural
diagnostic. `build_clean_row_training_frame` now logs pair-level drops so the 21 IAI64 all-`n` pair
drop is traceable across future dataset revisions.

**Audit confirms** (case-by-case): row counts exact (318,816 Guelin + 1,240 BASEL), Guelin drops
exactly 8,675 `n` rows (2.7%), BASEL drops 0 (no `n` category), 96 Guelin + 52 BASEL phage IDs fully
disjoint, all 25 BASEL ECOR bacteria resolve in the Guelin cv_group map, no train/test leakage,
every pair_id has exactly one source tag, CSV metrics reproduce JSON exactly. **21 benign drops** are
all on bacterium IAI64 (22% of IAI64's phage surface — any future IAI64 per-bacterium analysis needs
this caveat).

**Run.** `PYTHONPATH=. python -m lyzortx.pipeline.autoresearch.ch05_eval --device-type cpu` —
10,774 s (3 hr): 2 axes × 10 folds × 3 seeds on ~270–290K rows per fit.

**Results (BASEL encoded at absolute log₁₀ pfu/ml = 9.0, Guelin neat = 8.7).**

| Quantity | Value | 95% CI |
|---|---|---|
| Bacteria-axis AUC | 0.8063 | [0.7919, 0.8202] |
| Bacteria-axis Brier | 0.1778 | [0.1702, 0.1853] |
| Phage-axis AUC | 0.8849 | [0.8616, 0.9059] |
| Phage-axis Brier | 0.1349 | [0.1221, 0.1497] |
| Phage-axis Guelin subset AUC | 0.8861 | [0.8652, 0.9077] |
| Phage-axis BASEL subset AUC | 0.8829 | [0.8207, 0.9324] |
| Phage-axis Guelin subset Brier | 0.1330 | [0.1191, 0.1472] |
| Phage-axis BASEL subset Brier | **0.1881** | [0.1580, 0.2208] |
| Phage-axis generalization gap (bacteria − phage AUC) | −0.0785 | CIs disjoint |

Per-fold AUC bacteria-axis: 0.822, 0.783, 0.781, 0.835, 0.806, 0.816, 0.798, 0.840, 0.794, 0.781
(std 0.023). Per-fold phage-axis: 0.884, 0.938, 0.890, 0.900, 0.865, 0.914, 0.896, 0.895, 0.841, 0.875
(std 0.028).

**BASEL bacteria-axis deficit** (first-order finding, separate from phage-axis parity). Filtering
bacteria-axis predictions CSV by source: Guelin AUC 0.8098 (35,403 pairs), **BASEL AUC 0.7152**
(1,240 pairs) — BASEL is 9.5 pp below Guelin on the same axis, and 17 pp below BASEL phage-axis on
the same pairs. The aggregate bacteria-axis AUC 0.8061 is 96.6% Guelin-weighted at row level, which
is why this finding is invisible in the headline. It is a separate deployability diagnostic from the
cross-source phage-axis question and belongs in its own line on the knowledge model.

**Per-decile reliability** (observed lysis rate vs mean predicted probability; bins of width 0.1).
Regenerable via `ad_hoc_analysis_code/ch05_reliability_diagrams.py`.

| Bin | Guelin bacteria-axis | BASEL bacteria-axis | Guelin phage-axis | BASEL phage-axis |
|---|---|---|---|---|
| 0.5–0.6 | −22.6 pp | −20.0 pp | −23.5 pp | **−48.3 pp** |
| 0.6–0.7 | −28.1 pp | **−45.3 pp** | −26.3 pp | **−53.0 pp** |
| 0.7–0.8 | −30.4 pp | **−51.6 pp** | −20.8 pp | **−46.1 pp** |
| 0.8–0.9 | −27.6 pp | **−49.9 pp** | −17.0 pp | **−44.1 pp** |
| 0.9–1.0 | −18.4 pp | −26.0 pp | −6.4 pp | +0.8 pp |

**Guelin is also substantially miscalibrated in mid-P — not just BASEL.** Guelin bacteria-axis
predicts 0.75 but observes 0.45 in the 0.7–0.8 bin (−30.4 pp gap), and predicts 0.65 but observes
0.37 in the 0.6–0.7 bin (−28.1 pp). The "BASEL 21–27 pp wider than Guelin" framing is accurate in
relative terms, but Guelin is 20–30 pp off in absolute terms; BASEL is another 20 pp worse on top
of that. Read both together: the model is systematically too confident in mid-P on **both**
sources; BASEL just amplifies the same failure mode. BASEL's top-decile phage-axis calibration is
fine (+0.8 pp at 0.9–1.0), so the model *can* calibrate when features unambiguously drive a high
score; it miscalibrates specifically in the mid-P region where features push predictions upward
but pairs don't actually lyse at the predicted rate.

**Isotonic calibration diagnostic (two separable mechanisms, empirically confirmed).** The
reviewer flagged that attributing all mid-P miscalibration to TL17-bias stretches the
feature-failure-mode story to cover territory it doesn't explain — TL17-bias concerns
cross-source generalization, but Guelin is also 28–30 pp off in mid-P despite the model being
trained on these exact phages. Tested this via leave-one-fold-out isotonic regression on Guelin
predictions (fold-aware, no leakage), then applied the Guelin-fitted calibrator to BASEL
predictions as the discriminative test — Guelin/BASEL panels are disjoint so no leakage. Script:
`ch05_isotonic_calibration_test.py`. Expected three outcomes (A: threshold story for both;
B: threshold for Guelin, TL17-bias extra for BASEL; C: threshold wrong for both).

| Subset | Raw ECE | Iso ECE | Raw max \|gap\| | Iso max \|gap\| | Closure |
|---|---|---|---|---|---|
| Guelin bacteria-axis | 0.120 | **0.008** | 30.4 pp | **6.6 pp** | 78% |
| Guelin phage-axis | 0.104 | **0.008** | 26.3 pp | **2.9 pp** | 89% |
| BASEL bacteria-axis | 0.272 | 0.113 | 51.6 pp | 32.6 pp | 37% |
| BASEL phage-axis | 0.236 | 0.122 | 53.0 pp | 35.1 pp | 34% |

**Outcome B confirmed.** Guelin mid-P is essentially fully fixable with an isotonic layer (ECE
drops from 0.12 to 0.008, gap closure 78–89%). BASEL improves but retains a substantial residual
(ECE stays at 0.11–0.12, gap closure only 34–37%) — the Guelin-fitted calibrator cannot rescue
BASEL's full miscalibration because it comes from a different mechanism. AUC sanity: moves by ≤0.5
pp on Guelin (isotonic produces ties that shift rank-based ordering slightly) and ≤0.001 pp on
BASEL; the small Guelin movement is tie-breaking noise, not a bug — isotonic is monotone but not
strictly monotone. Strict AUC invariance holds only for unclipped strictly-monotonic maps; under
boundary clipping and ties in the fitted step function, ranking is preserved up to tie-breaking
at the edges, which can produce sub-percent AUC movement. Worth flagging for transparency.

**Two-driver story.** Guelin's mid-P gap is a training-label-vs-deployment-question mismatch.
Gaborieau 2024 Methods (see `mlc-dilution-potency`) explicitly says clearing events at high titer
"could result from productive lysis of the bacterial population by the phage, or from another
mechanism such as lysis from without, or abortive infection." CHISEL's per-row binary training
treats every score='1' row as a positive, including non-productive clearing. The deployment
question — "will this phage productively lyse this strain at therapeutic dose?" — is stricter
than the training label, and the model's inflated mid-P probabilities are the signature of that
mismatch. This is fixable post-hoc without changing the underlying model (the discrimination is
already there). Connected to `ambiguous-label-noise` (GT09 showed +3.1 pp top-3 from dropping 10%
of noisy rows — same shape of mechanism, different regime).

BASEL's residual after Guelin-fitted calibration (~12 pp ECE, ~33 pp max decile gap) is the part
that threshold-mismatch cannot explain. It matches the `ch05_basel_zero_vector_split` diagnostic:
miscalibration concentrates on the 39/52 non-zero-projection BASEL phages (Brier 0.31 bacteria-
axis) whose phage_projection vectors map into Guelin-calibrated TL17 neighborhoods associated
with broad-host lysis but whose actual host ranges are narrower. That is a feature-transfer
problem that needs panel-independent phage features (CH06), not a calibration layer.

**Implication for CH06 scope** (forward-looking, not a CH05 acceptance criterion): BASEL's
bacteria-axis AUC 0.7152 vs Guelin 0.8098 is the 9.5 pp discrimination gap; no calibration layer
can touch it (ranking is preserved under monotone transforms). CH06's acceptance criterion —
"materially above 0.7152 on at least one arm" — still reads as written. The calibration component
of BASEL's 0.1884 Brier vs Guelin's 0.1329 is now understood to be partially fixable post-hoc;
the remaining discrimination gap is what CH06 must attack.

**Follow-up** (deferred, out of CH05 scope): a dedicated CH-series ticket on "post-hoc
calibration layer + label-threshold sensitivity" (likely CH09 after CH06 and CH07). Platt scaling
or isotonic as a deployable calibration step is cheap and independent of CH06's feature work.

**Per-family phage-axis breakdown (post-hoc).**

| Family | Regime | Panel size | n_pairs | pos_rate | AUC | Brier |
|---|---|---|---|---|---|---|
| Autographiviridae | within-family | 47 | 14,101 | 0.162 | 0.8678 | 0.1088 |
| Other | within-family | 40 | 14,757 | 0.345 | 0.8948 | 0.1400 |
| Straboviridae | within-family | 24 | 4,261 | 0.502 | **0.7338** | 0.2400 |
| Drexlerviridae | within-family | 17 | 1,770 | 0.080 | 0.9251 | 0.0629 |
| Demerecviridae | family-novel (rare) | 9 | 148 | 0.122 | 0.7872 | 0.3286 |
| Schitoviridae | family-novel (rare) | 5 | 1,131 | 0.118 | 0.8526 | 0.0830 |

Aggregate 0.8850 mixes regimes. Straboviridae 0.7338 within-family is the striking weak case (narrow-
host Straboviridae pattern from `family-bias-straboviridae` surviving the CHISEL frame change).
Drexlerviridae 0.9251 is standout strong (easy-negative bias at 8% positive). Family-novel AUCs
(0.7872 / 0.8526) sit inside the within-family spread, so "siblings-in-training boosts AUC" isn't
clean — signal quality varies more by family identity than by generalization regime. CH06 both-axis
double-CV is the cleaner test for deployability. Regenerable via
`lyzortx/research_notes/ad_hoc_analysis_code/ch05_per_family_phage_axis.py`.

**Cross-source parity is narrower than originally framed.**

The earlier draft claimed AUC parity "confirms SX15's 11 pp nDCG gap was a metric artifact." That
overclaims — AUC and nDCG responding differently does not *prove* one was an artifact; it shows the
two metrics are sensitive to different aspects of the prediction distribution. CH05 data actually
supports three separate findings, not a single parity statement:

- **Phage-axis discrimination parity**: |AUC gap| 0.0032 << 1 pp. BASEL CI is 3× wider than Guelin's
  on 52 phages vs 96 — "indistinguishable" is a weak non-rejection on small N, not positive evidence
  of transfer.
- **Phage-axis calibration divergence**: Guelin Brier 0.1329 vs BASEL 0.1884, disjoint CIs, 42%
  relative calibration hit on BASEL. AUC parity does not imply calibration parity.
- **Bacteria-axis cross-source asymmetry**: BASEL 0.7152 vs Guelin 0.8098 on the same axis (9.5 pp
  deficit); BASEL transfers materially worse when held out as bacteria than when held out as phages.

The nDCG-vs-AUC narrative narrows to: under AUC, cross-source discrimination parity holds *on
phage-axis only*; the nDCG gap was likely a graded-vs-binary label artifact (BASEL has no MLC
grades), not a deployment signal. Deeper interpretation is case-by-case post-hoc (below).

**Straboviridae prior collapse ruled out as primary driver.** A post-hoc diagnostic
(`ad_hoc_analysis_code/ch05_straboviridae_exclusion.py`) excluding Straboviridae phages from the
predictions closes only 1.5 pp of the 9.5 pp bacteria-axis BASEL deficit, and near-zero on
phage-axis. BASEL Brier still 6.3 pp worse than Guelin even with Straboviridae removed. Future
tickets should not relitigate this hypothesis.

**Root cause: cross-panel phage feature failure mode.** Under the phage-axis split,
non-zero-projection BASEL phages (n=39) are substantially miscalibrated — mean P|y=0 = **0.55**
(bacteria-axis) / 0.44 (phage-axis), Brier **0.31** / 0.21 — while zero-projection BASEL phages
(n=13) calibrate well (mean P|y=0 = 0.17, Brier 0.12 bacteria-axis). The failure is
distributional: non-zero BASEL phages map into TL17 projection neighborhoods populated by
broader-host Guelin phages, and the model applies Guelin-calibrated host-range priors that don't
fit BASEL's narrower actual host ranges. Zero-vector BASEL phages calibrate correctly precisely
because the model has *no* phage signal to apply — it falls back to the host-side prior, which is
roughly correct for their 9.7% positive rate. **Phage features from a Guelin-only reference bank
are therefore actively harmful on out-of-distribution phages, not merely neutral.**

This connects to `plm-rbp-redundant` (same mechanism, cross-family rather than cross-panel: PLM
phage features hurt cross-family predictions by −1.22 pp nDCG because they inject within-family
priors the model misapplies to unseen families — here we see the cross-panel version of the same
effect). It also supports `panel-size-ceiling`: the fix is not "fill in 13 zero-vector phages" or
"engineer richer phage features" — it is more training phages at currently-underpopulated TL17
neighborhoods, so Guelin-calibrated priors aren't the only signal available for BASEL-like phages.
Cheap diagnostic evidence: `ch05_basel_feature_variance.csv` shows 36.4% of `phage_projection`
features are constant across all 52 BASEL phages vs 0% across Guelin (i.e., even the 39 non-zero
BASEL phages occupy a narrower TL17 subregion than Guelin). Zero-vector ∩ UNKNOWN-family overlap is
partial (8/13), so taxonomic novelty and TL17 coverage correlate but are not coincident —
`ch05_basel_zero_vector_split.csv` documents the full phage list and per-subset metrics.

**Hypothesis ruled weaker by pos-rate check**: label-semantics asymmetry (Guelin productive lysis
vs BASEL spot clearing at >10⁹ pfu/ml, which Maffei 2021 says can include "lysis from without or
abortive infection"). If BASEL's label were materially more permissive we'd expect its positivity
rate to exceed Guelin's. Non-zero BASEL pos_rate 29.5% vs Guelin 27.3% — essentially equal. Label
semantics may still contribute at the margin but is not a first-order driver. Not promoted to a
live candidate without a label-remap experiment that is currently deprioritised.

**Encoding correctness (addressed this ticket, not a driver).** BASEL's >10⁹ pfu/ml spot was
previously encoded at `log_dilution=0` (= Guelin neat 5×10⁸ pfu/ml, i.e. ~2× below BASEL's actual
titer) because the original `pair_concentration__log_dilution` feature only represented
relative-dilution steps. CH05's reliability analysis shows BASEL is *over*-predicted in mid-P —
the opposite direction the encoding hypothesis predicts — so the encoding was not the dominant
driver. Regardless, the encoding has been fixed track-wide to absolute `log10_pfu_ml` (Guelin
{4.7, 6.7, 7.7, 8.7}; BASEL 9.0, the Maffei-reported lower bound — Maffei 2021 Fig. 12 and
Maffei 2025 Fig. 13 both quote >10⁹ pfu/ml as the working titer, 2025 adds "if possible"). Guelin
neat at 10⁸·⁷ and BASEL at 10⁹·⁰ now coexist on one feature axis without implicitly mapping BASEL
onto Guelin's neat. CH04 rerun under the new encoding is an affine-shift sanity check (metrics
expected within ≤1e-6 of prior run — LightGBM pre-bins features before split search, so a monotonic
affine shift on a single feature yields identical bin boundaries and tree structure). CH05 rerun
under the new encoding is deferred to a follow-up ticket — the encoding fix is scientifically small
but semantically load-bearing for future cross-source work.

**Per-axis gap interpretation.** 7.9 pp bacteria → phage AUC advantage is a **training-coverage
structural effect**, not a model-quality signal. Bacteria-axis holds out bacteria (loses host-side
training signal for test pairs); phage-axis holds out phages but keeps all 369 bacteria in training
(full host-side signal per test pair). Do NOT read phage-axis 0.8850 as "better" than bacteria-axis
0.8061 — they answer different deployability questions. SX15 framed this as "per-phage blending tax";
under CHISEL's all-pairs-only model there is no blending to tax and the gap is purely structural.

**Artifacts** (all under `lyzortx/generated_outputs/ch05_unified_kfold/`):
`ch05_combined_summary.json`, `ch05_{bacteria,phage}_axis_metrics.json`,
`ch05_cross_source_breakdown.csv`, `ch05_{bacteria,phage}_axis_predictions.csv` (with source tag),
`ch05_{bacteria,phage}_axis_per_row_predictions.csv`, `ch05_per_family_breakdown.csv`,
`ch05_straboviridae_exclusion.csv`, `ch05_reliability_tables.csv`, `ch05_basel_feature_variance.csv`,
`ch05_basel_zero_vector_split.csv`. Bootstrap JSON now exposes `bootstrap_samples_requested` and
`bootstrap_samples_used` so degenerate-resample skips are visible in the output.

**Follow-ups deferred to a dedicated ticket**: (a) rerun CH05 end-to-end under the new
`log10_pfu_ml` encoding to replace headline numbers; (b) probe TL17 reference-bank bias on the 39
non-zero BASEL phages (candidate driver #1); (c) pilot a label-remap experiment if Gibbs-style
retest data becomes available for BASEL (candidate driver #2). CH06 both-axis double-CV remains the
deployability test. Family-novel generalization is measured imperfectly here (rare-family collapse
mixes regimes); IAI64 per-bacterium analysis carries the 22% phage-surface caveat. BASEL with a
larger panel would tighten cross-source CI comparison.

**Acceptance met.**

- Unified 148 × 369 panel at row level; BASEL embedded with absolute `log10_pfu_ml` encoding (9.0,
  Maffei lower bound), cv_group inherited (hard-fail on missing mapping).
- Bacteria-axis CV under CH02 cv_group hash; phage-axis CV under ICTV family + "Other" (<10
  phages) + "UNKNOWN" (no family) catch-all 10-fold stratification (40% of folds are pseudo-family
  buckets; the earlier "ICTV-stratified" shorthand was misleading).
- Cross-source AUC+Brier with CIs; discrimination parity, calibration divergence, AND BASEL
  bacteria-axis deficit reported as three separate findings (not a single parity headline).
- Per-decile reliability tables for both sources × both axes.
- Straboviridae-exclusion diagnostic, per-phage feature-variance comparison, and BASEL zero-vector
  vs non-zero split all regenerable from existing artifacts (no model rerun).
- AUC + Brier only (top-3/nDCG retired under `ranking-metrics-retired`).
- Per-family phage-axis breakdown post-hoc.
- Pair-drop logging added (the 21 IAI64 all-`n` pair drop is visible across future dataset revisions).
- `chisel-unified-kfold-baseline` knowledge unit active with narrower claims; `spandex-unified-kfold-baseline` HISTORICAL.
- **Ticket scope change** (explicit): per-phage "blending tax" retired under
  `per-phage-retired-under-chisel`; replaced by phage-axis generalization gap as a structural
  diagnostic, not a quality metric.
- 7 CH05 unit tests (phage-fold determinism, rare-family collapse, missing-family handling, BASEL
  cv_group inheritance, fail-fast on missing map, cluster-bootstrap unit-level behavior,
  degenerate-resample skip). Vacuous `test_source_constants_are_distinct` removed per self-review.
- Bootstrap JSON exposes `bootstrap_samples_requested`/`bootstrap_samples_used`.
- Concentration encoding reworked track-wide: `CONCENTRATION_FEATURE_COLUMN` is now
  `pair_concentration__log10_pfu_ml`; CH03 derives the Guelin value; CH05 sets BASEL directly;
  `--basel-log10-pfu-ml` CLI flag available for sensitivity analysis.
- Full CH05 rerun under the new encoding deferred to a follow-up ticket (scientifically small
  affine shift; notebook/knowledge headline numbers remain pre-encoding-fix and are flagged as such).
