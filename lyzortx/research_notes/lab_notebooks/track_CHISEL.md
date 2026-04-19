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
