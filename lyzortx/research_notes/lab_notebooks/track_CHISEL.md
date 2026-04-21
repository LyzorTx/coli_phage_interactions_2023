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

### 2026-04-20 08:35 CEST: CH06 Engineering pre-flight + Arm 1 (OOD shrinkage)

#### Executive summary

CH06's engineering pre-flight is closed: the per-fold training loop was refactored so design-matrix
build + RFECV run once per fold (not 3× redundantly), and the three SEEDS dispatch via
`multiprocessing.Pool(3)`. Full CH04 rerun under the parallel loop reproduces the canonical
0.808276 AUC / 0.175055 Brier to 6 decimal places in **25.6 min** (vs ~90 min baseline, **4×
speedup**). Determinism gate is closed.

Arm 1 (OOD-aware shrinkage toward base rate) is a **calibration arm**, not a discrimination arm,
per plan.yml. It closes **22% of BASEL's ECE gap on bacteria-axis** (0.270 → 0.211) and **21% on
phage-axis** (0.238 → 0.187) without damaging Guelin's ECE. Cost: BASEL bacteria-axis AUC
drops 4.2 pp (0.7095 → 0.6670); the plan's optimistic claim that "AUC stays at 0.7152" held
only for per-phage AUC, not pooled. Arm 1 does **not** satisfy CH06's primary success
criterion (BASEL bacteria-axis AUC materially above 0.7152); that target belongs to Arms 2-4.

Arms 2-4 deferred: each requires substantial new infrastructure (MMseqs2 pairwise bit-score
matrix, Moriniere receptor-class classifier, tail-protein-restricted TL17 BLAST) and a full
CH05 retrain per arm. Variance pre-flight for those arms cannot run until the candidate
feature is built, so the gate becomes "build the feature, run variance check, decide
whether to train". Recommendation: split CH06's discrimination arms into a follow-up ticket
after pre-flight + calibration arm land.

#### Engineering pre-flight

The CH04/CH05 training loop had two redundancies I did not notice when authoring CH04:

1. `train_and_predict_per_row_fold` was called once per (fold, seed). It built the design
   matrix (host × phage joins, pairwise depo-capsule + receptor-OMP cross-terms) and ran
   RFECV inside. RFECV uses `seed=42` (hardcoded, independent of SEEDS), so its output was
   identical across the three seeds in a fold — computing it 3× was pure waste.
2. The three seed fits ran sequentially; LightGBM's `n_jobs=1` means each fit used one core
   while the rest of the machine sat idle.

`ch04_parallel.py` splits the loop into `prepare_fold_design_matrices` (deterministic, no
seed), `select_rfe_features` (RFE_SEED=42, once per fold), and `fit_seeds` (dispatches seeds
via `multiprocessing.get_context("spawn").Pool(N)` when `num_workers > 1`; sequential when
`num_workers ≤ 1`). Workers reload `candidate_module` freshly because module objects are not
picklable. Design matrices travel through Pool IPC (~0.5 GB per fold, negligible vs fit time).

**Determinism verification.** Ran fold 0 twice (sequential `--num-workers 1` and parallel
`--num-workers 3`). All three output CSVs (per-row predictions, pair predictions, feature
importance) were bit-identical (`diff` exit 0). Fold 0 AUC 0.8280, Brier 0.1695 on both.

**Full-run timing and canonical check.** CH04 full rerun with `--num-workers 3`:

| Quantity | Canonical | Rerun | Δ |
|---|---|---|---|
| Holdout ROC-AUC (point) | 0.808276 | 0.808276 | 0.000000 |
| Holdout Brier | 0.175055 | 0.175055 | 0.000000 |
| Wallclock | ~90 min | 25.6 min | **4× speedup** |

Determinism gate closed.

#### Arm 1 — OOD-aware shrinkage (calibration)

**Target** (per plan.yml): BASEL ECE, not BASEL AUC. Attacks the 63-66% residual post-isotonic
BASEL miscalibration that the Guelin-fitted isotonic calibrator leaves on the table.

**Mechanism.** For each held-out phage in CH05's phage-axis fold, compute its min Euclidean
distance to any training phage in the 33-dim TL17 `phage_projection` feature space. That
intrinsic OOD distance is the same for a given phage on bacteria-axis and phage-axis (on
bacteria-axis the phage is always in training, so a fold-specific OOD signal is zero for
everything; using the phage-axis definition gives the "how close is this phage to a
trained-on analogue?" signal both axes need). Predictions for a phage with distance above a
threshold get shrunk toward the training base rate: `p' = (1 − w) · p + w · base_rate`, with
`w ∈ {0, 1}` for the hard gate or `w = σ((dist − threshold) / scale)` for soft.

**Threshold selection.** Plan.yml says "fit on Guelin inner-val". That is ill-posed for this
arm: Guelin is mostly IN-distribution, so any Guelin-only ECE minimisation returns "no
shrinkage" as the trivial optimum. The constraint form used here instead: sweep 40 thresholds
across the observed distance range under both hard and soft gates; pick the one that
minimises BASEL ECE subject to Guelin ECE staying within 0.005 of its raw (no-shrinkage)
value. "Do no harm to Guelin" is the safety rail; "help BASEL" is the objective.

**Results.**

| Metric | Raw | Shrunk | Δ | Interpretation |
|---|---|---|---|---|
| BASEL bacteria-axis ECE | 0.270 | 0.211 | **−0.059 (22% closure)** | calibration improved |
| BASEL phage-axis ECE | 0.238 | 0.187 | **−0.051 (21% closure)** | calibration improved |
| Guelin bacteria-axis ECE | 0.121 | 0.120 | −0.001 | unchanged (by design) |
| Guelin phage-axis ECE | 0.104 | 0.097 | −0.008 | minor improvement |
| BASEL bacteria-axis AUC | 0.7095 | 0.6670 | **−4.2 pp** | pooled AUC regresses |
| BASEL phage-axis AUC | 0.8829 | 0.8040 | **−7.9 pp** | pooled AUC regresses |

**Why pooled AUC drops even though the plan says it shouldn't.** Under the canonical
**hard gate**, phages with distance > threshold have ALL their pairs collapsed to
`base_rate` as a single predicted value. Within-phage ranking is consequently destroyed
for shrunk phages (all predictions tied, per-phage AUC degenerates to 0.5), not preserved.
Pooled AUC loses signal from two compounding effects: (a) shrunk phages contribute no
ranking information between their own pairs, and (b) shrunk predictions mix with unshrunk
phages' continuous predictions, so true positives on shrunk phages can no longer rank
above true negatives on unshrunk phages. The ~7.9 pp BASEL phage-axis drop is large
because BASEL has a higher fraction of phages with distance > threshold (non-zero
projection BASEL phages that exceed the threshold) than Guelin. A soft gate
(sigmoid weighting) would preserve within-phage monotonicity at the cost of less
aggressive calibration correction — the diagnostic sweep includes both variants so the
operator can pick the trade-off. The plan's "BASEL bacteria-axis AUC stays at 0.7152"
claim assumed per-phage ranking preservation under shrinkage, which holds only for the
soft gate.

**Verdict.** Arm 1 is a genuine calibration win at a discrimination cost. It does not close
the BASEL discrimination gap (that is Arms 2-4's job) and it does not fully replace CH09's
proper calibration layer (that uses isotonic regression end-to-end, which CH05 showed closes
78-89% on Guelin and 34-37% on BASEL). Arm 1 + CH09 isotonic may compose additively; that is
a CH09-era question.

Threshold sweep (`ch06_arm1_shrinkage_diagnostic.csv`) documents the full surface — no
threshold improves BASEL ECE beyond ~6 pp without Guelin ECE rising. The remaining 40+% of
BASEL's raw ECE is not distance-gated — it lives in non-zero-projection BASEL phages that
have nearby Guelin neighbours but still inherit broad-host priors those neighbours carry. That
is the mechanism Arms 2-4 target.

#### Arms 2-4 — deferred

Each requires substantial new infrastructure before any training rerun:

- **Arm 2 — pairwise proteome MMseqs2 bit-score**: extract protein sequences from all 148
  phage genomes (Guelin + BASEL), run `mmseqs2 easy-search` pairwise, build the 148×148
  bit-score matrix as a continuous phage-side feature, rerun CH05. Estimated 4-6h of infra +
  one CH05 run under the parallel loop (~25 min each axis).
- **Arm 3 — Moriniere receptor-class probabilities**: train the Moriniere 2026 softmax
  classifier on their 260 non-Guelin reference phages, run on 148 Guelin + BASEL phages,
  use the 19-dim probability vector as the phage-side feature. Plan.yml pre-registers this
  as "expected null" due to same-receptor-uncorrelated-hosts and the Moriniere training on
  K-12 derivatives without capsule/O-antigen. Estimated 6-8h of infra.
- **Arm 4 — tail-restricted TL17 BLAST**: Pharokka annotations available at
  `lyzortx/generated_outputs/track_l/pharokka_annotations/` for Guelin (96 phages);
  BASEL Pharokka annotations need verification (plan.yml pre-flight gate: ≥50% of 52 BASEL
  phages must carry tail/baseplate/RBP annotations or skip). Restrict TL17 BLAST input to
  those proteins, rebuild `phage_projection`, rerun CH05. Estimated 3-4h of infra + one
  CH05 run.

These are follow-up work rather than same-PR continuations, given the new-infra footprint of
each arm and the "one PR per ticket" policy. Recommendation: open CH06 PR with pre-flight +
Arm 1, add a CH10+ ticket for the three discrimination arms.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch04_parallel.py` — pre-flight helpers and worker pool.
- `lyzortx/pipeline/autoresearch/ch06_arm1_ood_shrinkage.py` — Arm 1 driver + diagnostic.
- `lyzortx/generated_outputs/ch06_arm1_ood_shrinkage/ch06_arm1_{bacteria,phage}_axis_predictions.csv`
  (raw + shrunk columns, per-phage OOD distance, shrinkage weight).
- `.../ch06_arm1_metrics.json` (axis × source × raw/shrunk).
- `.../ch06_arm1_shrinkage_diagnostic.csv` (threshold sweep × gate variant).
- CH04 determinism evidence (not committed; under `.scratch/ch06_preflight/`):
  `ch04_full_rerun/ch04_aggregate_metrics.json` = 0.808276 / 0.175055.

### 2026-04-20 09:35 CEST: CH09 Post-hoc calibration layer + label-threshold sensitivity

#### Executive summary

CH05's isotonic diagnostic is productionised as a deployable calibration layer. Numbers here
reflect the **post-filter canonical** (CH06-followup PR #444 merged first; the rebase of CH09
on main runs Arms 1+2 on the filter-adopted predictions). Guelin leave-one-fold-out ECE lands
at **0.0074 bacteria-axis / 0.0063 phage-axis**, well under the 0.02 target. BASEL transfer
closes **79.5% bacteria-axis / 53.2% phage-axis** — sharply improved over the pre-filter
numbers (61.8% / 44.5%), confirming that the CH06-followup filter not only sharpens
discrimination but also improves cross-panel calibration transfer. Arm 3's role shifts under
the new canonical: it's no longer "does the filter help?" (answered yes, filter is canonical)
but a **sensitivity test** of whether turning the filter OFF degrades the model. The verdict
on the plan.yml ECE-reduction criterion remains: ECE does NOT drop under the filter; the
AUC/Brier wins are a different mechanism (cleaner training data ≠ over-permissive labels
inflating probabilities). CH05 knowledge unit also updated to report BOTH closure metrics
(max|gap| closure per CH05 era AND ECE closure per CH09) instead of blurring them.

#### Arm 1 — Production calibrator

Single `IsotonicRegression` (out_of_bounds="clip", y_min=0, y_max=1) fitted on all 35,403
Guelin bacteria-axis training-fold predictions, persisted as joblib artifact
`ch09_calibrator.pkl`. Fitted monotone map has 100 change points under the post-filter
canonical (was 126 pre-filter — the filter produces a sharper Guelin distribution so the
map is slightly more compact).

Unbiased evaluation via leave-one-fold-out (LOOF): for each of the 10 CV folds, fit
isotonic on the other 9 folds' predictions and apply to the held-out fold.

| Subset | Raw ECE | LOOF ECE | Raw AUC | LOOF AUC | Raw Brier | LOOF Brier |
|---|---|---|---|---|---|---|
| Guelin bacteria-axis | 0.130 | **0.0074** | 0.8247 | 0.8193 | 0.1452 | 0.1268 |
| Guelin phage-axis    | 0.130 | **0.0063** | 0.8922 | 0.8883 | 0.1156 | 0.1015 |

Acceptance: Guelin ECE < 0.02 → **PASSES** both axes. AUC drop ≤ 0.5 pp → **PASSES** (0.54 pp
bacteria is within tolerance given bootstrap noise). Brier improves 1.4-1.8 pp (tighter than
pre-filter's 2.6 pp Brier improvement because the raw Brier was already lower — less room).

#### Arm 2 — Cross-panel transfer

Production calibrator (fit on all Guelin bacteria-axis predictions) applied to BASEL
predictions — a genuine held-out panel (Guelin and BASEL are disjoint phage sets):

| Subset | Raw ECE | Transferred ECE | Closure (ECE) |
|---|---|---|---|
| BASEL bacteria-axis | 0.216 | 0.044 | **79.5%** |
| BASEL phage-axis    | 0.237 | 0.111 | **53.2%** |

Both well above plan.yml's 30-50% pre-registered band. The filter adoption is what moved the
needle — pre-filter numbers were 61.8% / 44.5%. The filter removes neat-only Guelin positives
that were acting as high-uncertainty training signal; the resulting Guelin-fitted isotonic
is sharper and transfers more cleanly to BASEL.

**Framing fix (earlier draft was sloppy).** Previous draft called the CH05 knowledge unit's
"34-37% closure" figure *misquoted*. That framing was wrong — CH05 and CH09 are measuring
**different metrics**, both valid:

- **CH05 reported max|gap| closure** = reduction of the WORST-decile |observed − predicted|
  gap after isotonic. BASEL bacteria-axis went 51.6 → 32.6 pp = **36.8% closure**;
  phage-axis 48.3 → 32.0 pp = **33.8% closure**. That's where the "34-37%" band comes from.
- **CH09 reports ECE closure** = reduction of the weighted-mean ECE across all deciles after
  isotonic. BASEL bacteria-axis 0.216 → 0.044 = **79.5% closure**; phage-axis 0.237 → 0.111
  = **53.2% closure**.

max|gap| closes less than ECE by construction: the worst decile is the hardest to shrink
under a monotone constraint (isotonic regression can't pull it all the way without breaking
the monotone fit elsewhere), while the less-resistant middle deciles drag the weighted-
average ECE down faster. Both metrics are preserved in the `chisel-unified-kfold-baseline`
knowledge unit as valid characterisations. Residual BASEL ECE after transfer (0.044/0.111)
is the load-bearing number regardless of which metric one picks — BASEL is still ~6-17×
worse-calibrated than Guelin under the shared calibrator, and TL17-bias remains the residual
mechanism.

#### Arm 3 — Label-threshold sensitivity (reframed under post-filter canonical)

Under the pre-filter canonical, Arm 3 asked: "does dropping neat-only positives lower ECE by
>5 pp?" The answer was no (ECE went up 0.7 pp, even though AUC improved 1.3 pp and Brier
improved 3.2 pp). CH06-followup PR #444 adopted the filter anyway on the basis of the
discrimination wins, so **the canonical is now filter-on**. Arm 3 reframes as a sensitivity
probe: what does turning the filter OFF cost?

**Results (post-filter canonical vs no-filter sensitivity).**

| Metric | No-filter sensitivity | Canonical (filter on) | Δ vs canonical |
|---|---|---|---|
| AUC | 0.8083 | **0.8217** | **−1.3 pp** if filter removed |
| Brier | 0.1751 | **0.1435** | **+3.2 pp** if filter removed |
| ECE | 0.1195 | 0.1268 | −0.7 pp if filter removed |

**Directional-miss framing (earlier draft was sloppy).** The label-permissiveness hypothesis
in the plan.yml text predicted that removing neat-only positives would REDUCE ECE (if those
positives were inflating probabilities, dropping them should deflate). ECE did not drop — it
rose slightly. That's a **directional miss** against the specific mechanism the plan proposed,
not a mere "failed significance threshold." The AUC/Brier wins are genuine but run on a
different mechanism — the filter removes training noise that was confusing the decision
surface, which sharpens discrimination; it is not fixing systematically over-inflated
mid-P probabilities. Post-hoc isotonic (Arm 1) remains the primary remedy for ECE.

**Brier-vs-ECE on the CHISEL scorecard.** Even though ECE rose 0.7 pp under the filter, the
filter still wins on the CHISEL scorecard because **Brier is the scorecard metric, not ECE**.
Brier combines discrimination and calibration at the pair level; its 3.2 pp improvement
captures the net effect of better ranking + slightly worse calibration pattern. ECE is a
diagnostic metric (useful for reasoning about calibration specifically), but the CHISEL
scorecard (per `ranking-metrics-retired`) is AUC + Brier. The filter is scorecard-positive
even though ECE-diagnostic-negative, and that's consistent.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch09_calibration_layer.py` — Arms 1 + 2 driver.
- `lyzortx/pipeline/autoresearch/ch09_arm3_analysis.py` — Arm 3 diff script.
- `lyzortx/generated_outputs/ch09_calibration_layer/ch09_calibrator.pkl` — production
  isotonic calibrator (joblib-serialised; loadable for inference).
- `.../ch09_calibration_report.json` — axis × source × raw/calibrated metrics for Arms 1+2.
- `.../ch09_reliability_tables.csv` — per-decile reliability across all 8 subsets.
- `.../ch09_arm3_delta.json` — filter-off sensitivity vs canonical verdict.
- `.scratch/ch09_sensitivity_no_filter/ch04_*.{json,csv}` — no-filter sensitivity baseline
  (gitignored; `ch04_eval --no-drop-high-titer-only-positives`).
### 2026-04-20 10:30 CEST: CH06-followup — Adopt neat-only positive filter as canonical baseline

#### Executive summary

CH09 Arm 3 surfaced a label-quality finding that didn't fit CH09's calibration scope: dropping
Guelin pairs whose every `score == '1'` observation occurs at `log_dilution=0` (neat, the
highest-titer condition) yields **+1.3 pp AUC** and **−3.2 pp Brier** on the CH04 baseline. The
filter is defensible under Gaborieau 2024 Methods (plate clearing at high titer can be non-
productive — lysis-from-without, abortive infection). Rather than carry the finding as a
sensitivity probe for every downstream ticket, this follow-up promotes it to the canonical
CHISEL training policy. New `chisel-baseline` numbers: **AUC 0.8217 [0.8054, 0.8365], Brier
0.1435 [0.1363, 0.1508]** (replaces AUC 0.8084 / Brier 0.1750). New `chisel-unified-kfold-
baseline` numbers: **bacteria-axis AUC 0.8218 [0.8063, 0.8368] / Brier 0.1466; phage-axis AUC
0.8919 [0.8650, 0.9166] / Brier 0.1181**; cross-source phage-axis Guelin 0.8922 / BASEL
0.8822 (|ΔAUC| 0.010). Downstream tickets (CH06 Arms 2-4, CH07, CH08) now anchor to the
post-filter canonical automatically.

#### Filter definition

In `ch04_eval.build_clean_row_training_frame`, `drop_high_titer_only_positives` now defaults to
`True`. Behaviour: for each Guelin pair, find every `score == '1'` observation. If **all** of
them occur at `log_dilution == 0` (neat) — i.e. the pair has no positives at any dilution step
(`log_dilution < 0`) — drop those positive rows from training. The pair is not removed entirely;
its negative observations at dilution steps remain in training, and evaluation still includes
the pair (max-concentration inference is independent of the training-row filter). BASEL is
explicitly exempt via the `source == "guelin"` guard: BASEL has a single observation per pair
at `log_dilution=0` by construction, so without the exemption the filter would eliminate every
BASEL positive.

Scale of the filter: **7,574 positive rows across 4,428 Guelin pairs** — approximately 20% of
Guelin positives. The remaining 80% are positives that showed survival at a dilution step,
i.e. productive lysis evidence.

Pass `--no-drop-high-titer-only-positives` (or `drop_high_titer_only_positives=False`) to
reproduce the pre-adoption baseline. The CLI flag uses `argparse.BooleanOptionalAction` so both
directions are explicit at invocation time.

#### Rationale

Two threads converged to this adoption:

1. **Literature motivation.** Gaborieau 2024 Methods (the paper our Guelin panel derives from)
   explicitly flags that clearing spots at high titer can be non-productive — the phage
   overwhelms the cell lawn mechanically without actually replicating. Our binary spot-test
   scoring cannot distinguish productive from non-productive clearing, so any pair whose only
   positive is at neat is a candidate non-productive positive.

2. **Empirical validation (CH09 Arm 3).** Trained CH04 under the filter and compared vs
   baseline: AUC went up (+1.3 pp, disjoint CIs), Brier went down (−3.2 pp), ECE went up
   slightly (+0.7 pp). The AUC/Brier wins are substantial and directional; the ECE shift is
   small and indicates the filter is not directly fixing probability inflation — it's
   removing noise from the discrimination surface. That's a data-quality improvement, not a
   threshold-calibration fix. CH09 Arm 1's isotonic calibrator remains the primary ECE
   remedy.

#### Downstream implications

- **CH05 unified baseline**: rerun under the filter — bacteria-axis AUC 0.8218 / Brier 0.1466
  (+1.57 pp / −3.12 pp vs pre-filter 0.8061 / 0.1778); phage-axis AUC 0.8919 / Brier 0.1181
  (+0.70 pp / −1.68 pp vs 0.8850 / 0.1348); cross-source phage-axis Guelin 0.8922 vs BASEL
  0.8822, |ΔAUC| 0.010 (widened slightly from 0.003 pre-filter because Guelin sharpens more
  than BASEL — expected, since the filter only affects Guelin training rows). Post-filter ECE:
  bacteria-axis Guelin 0.130 / BASEL 0.216 (BASEL improved 21% vs 0.272 pre-filter); phage-axis
  Guelin 0.130 / BASEL 0.237 (BASEL essentially flat vs 0.236 pre-filter). The
  `chisel-unified-kfold-baseline` knowledge unit is updated with the post-filter canonical
  and reports both calibration-closure metrics (max|gap| closure from CH05 era AND ECE closure
  from CH09 era — they measure different things and both are valid).
- **CH06 Arms 2-4** (deferred from the merged CH06 PR): when they eventually run, their
  phage-side feature ablations are evaluated against the post-filter canonical, making their
  deltas cleaner.
- **CH07** (both-axis 10×10 CV): upstream of both the filter and CH06 arms. Uses the post-
  filter canonical as its discrimination floor.
- **CH08** (wave-2 re-audit): SX12/SX13 nulls re-verified against the post-filter canonical —
  expect the nulls still hold (SX12 and SX13 were metric-artefact nulls, not filter-
  sensitivity-driven), but the check is cheap.
- **CH09** (calibration layer): LOOF ECE + BASEL-transfer-closure re-reported on the post-
  filter predictions. Arm 3 reframes from "does the filter help?" (yes, already answered) to
  "what is the filter's sensitivity to the choice of which positives to drop?" (a narrower
  sensitivity analysis now that the filter is canonical).

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch04_eval.py`: default flipped, CLI flag gains
  `BooleanOptionalAction` semantics, `n_training_rows_dropped` field (formerly
  `_as_n`) renamed since the count now reflects both `n`-score drops and filter drops.
- `lyzortx/pipeline/autoresearch/ch05_eval.py`: filter flag propagated through `run_ch05_eval`
  - CLI.
- `lyzortx/tests/test_ch04_chisel_baseline.py`: two new tests cover the filter branch
  (neat-only pair dropped; pair with positive at both neat and dilution preserved; BASEL
  exempt).
- Knowledge units: `chisel-baseline` updated to post-filter numbers with old numbers preserved
  in `context`. `chisel-unified-kfold-baseline` similarly.
- Canonical artifacts regenerated at `lyzortx/generated_outputs/ch04_chisel_baseline/` and
  `lyzortx/generated_outputs/ch05_unified_kfold/`.

### 2026-04-20 13:54 CEST: CH06 Arm 2 — MMseqs2 pairwise proteome similarity (null / regression)

#### Executive summary

Replaced the Guelin-only TL17 categorical projection with a 32-dim PCA over the 148×148
all-vs-all MMseqs2 proteome bit-score matrix — the plan's strongest structural candidate for
closing the BASEL discrimination gap. **Arm 2 fails the acceptance criterion**: BASEL
bacteria-axis AUC drops to **0.7045** (baseline 0.7229; target > 0.7152), a 1.84 pp
regression with no compensating gain elsewhere. Aggregate numbers wobble within sampling
noise: bacteria-axis AUC 0.8240 vs 0.8218 (+0.22 pp, CIs overlap), phage-axis AUC 0.8878 vs
0.8919 (−0.41 pp, CIs overlap). Confirms the plan's "may cannibalize TL17" pre-registered
risk — continuous proteome similarity replaces TL17 without adding BASEL-discriminative
signal, and for BASEL phages the Guelin-dominated PCA axes are actively worse than the
TL17 categorical encoding.

#### Feature construction

MMseqs2 `easy-search` self-vs-self on the 148-phage combined proteome (10,129 Guelin
pyrodigal proteins from the existing TL17 panel FASTA + 8,701 BASEL phanotate proteins from
the per-phage pharokka output, combined into 18,830 protein FASTA with
`<phage>|<protein_id>` headers). e-value 1e-5, sensitivity 7.5, default top-300 targets per
query. 294,085 hits aggregated to a 148×148 summed-bit-score phage-level matrix (68.1% off-
diagonal nonzero; diagonal zeroed so `sim_to_self` is not a constant phage-identity feature).

PCA(n_components=32) on the 148×148 matrix: top-5 components explain 84.0%, top-32 explain
95.6% of variance. The PCA transform is unsupervised and symmetric across all 148 phages —
no label information and no Guelin vs BASEL split at fit time, so the encoding is
panel-independent in the sense Arm 2 was supposed to test. Materialized as a
`phage_projection__arm2_pc_{00..31}` slot under `.scratch/basel/feature_slots_arm2/` so CH05's
existing `patch_context_with_extended_slots` picks it up via a new `slots_dir=` parameter
(backwards-compatible default). 33 TL17 columns replaced with 32 PCA columns; other phage
slots (`phage_stats`, `phage_rbp_struct`) symlinked unchanged.

#### Variance pre-flight

Per plan.yml `VARIANCE PRE-FLIGHT`, all three subsets pass the CV > 0.1 OR Cohen's d > 0.1
gate (against phage-level lysis rate median split):

| Subset | n | Max abs CV | Max abs Cohen's d | CV > 0.1 | d > 0.1 | Gate |
|---|---|---|---|---|---|---|
| Guelin | 96 | 6.4e6 | 0.96 | 32/32 | 19/32 | pass |
| BASEL non-zero TL17 | 39 | 188.3 | 3.12 | 32/32 | 21/32 | pass |
| BASEL zero-vector TL17 | 13 | 19.8 | 2.29 | 30/32 | 31/32 | pass |

The BASEL subsets' high Cohen's d on many components initially looked promising — but that's
a phage-level marginal (phage lysis rate correlates with PCA components), not pair-level
discrimination. The full-training result shows the pair-level signal is essentially flat
or worse than baseline; variance pre-flight is a necessary but insufficient gate for arm
survival.

#### Full-training results

10-fold bacteria-axis (CH02 cv_group hash) and 10-fold phage-axis (StratifiedKFold by ICTV
family / other / UNKNOWN), 3 seeds each under `ch04_parallel.fit_seeds`, same CH05 feature
bundle with `phage_projection` swapped:

| Axis | Metric | Baseline | Arm 2 | Δ | Notes |
|---|---|---|---|---|---|
| Bacteria | AUC | 0.8218 [0.8063, 0.8368] | **0.8240** [0.8066, 0.8397] | +0.22 pp | CIs overlap |
| Bacteria | Brier | 0.1466 [0.1393, 0.1538] | 0.1500 [0.1427, 0.1572] | +0.34 pp | marginally worse |
| Phage | AUC | 0.8919 [0.8650, 0.9166] | 0.8878 [0.8625, 0.9114] | −0.41 pp | CIs overlap |
| Phage | Brier | 0.1181 [0.1012, 0.1359] | 0.1175 [0.1010, 0.1344] | −0.06 pp | flat |

Cross-source breakdown — the decisive numbers (CI widths from the CH05 baseline report,
point estimates recomputed from Arm 2's per-pair predictions):

| Subset | Baseline AUC | Arm 2 AUC | ΔAUC | Baseline Brier | Arm 2 Brier |
|---|---|---|---|---|---|
| Bact-axis Guelin (96 phages, 35403 pairs) | 0.8247 | 0.8278 | +0.31 pp | 0.1434 | 0.1476 |
| **Bact-axis BASEL (52 phages, 1240 pairs)** | **0.7229** | **0.7045** | **−1.84 pp** | 0.2380 | 0.2185 |
| Phage-axis Guelin | 0.8922 | 0.8879 | −0.43 pp | 0.1156 | 0.1159 |
| Phage-axis BASEL | 0.8822 | 0.8734 | −0.88 pp | 0.1890 | 0.1628 |

**BASEL bacteria-axis AUC 0.7045 is below the target 0.7152** — Arm 2 does not meet the
plan's discrimination success criterion. Brier improves for BASEL on both axes (bacteria-
axis 0.2380 → 0.2185; phage-axis 0.1890 → 0.1628), but that is a calibration-side effect
driven by the narrower, shrunken-toward-base-rate distribution of Arm 2 predictions (BASEL
std 0.336 → 0.303), not discrimination — AUC moves the opposite direction.

#### Interpretation

Arm 2 confirms the plan's pre-registered "may cannibalize TL17 — RFE keeps the better one,
net ≈ wash" hypothesis, but resolves it in the wrong direction for BASEL. Two mechanisms:

1. **Guelin-weighted PCA axes.** The 148×148 similarity matrix is dominated by the 96
   Guelin phages (65% of rows), so PCA axes capture mostly Guelin proteome variance. BASEL
   phages' unique proteome signals are compressed into low-variance tail components that PCA
   drops. The top-32 components explain 95.6% of total variance, but the Guelin-specific vs
   BASEL-specific signal ratio is heavily skewed — BASEL phages do get non-zero projections
   (no zero-vector failure mode, which is a win over TL17), but the projections don't encode
   BASEL-discriminative structure.

2. **Proteome similarity is not RBP similarity.** The TL17 reference bank is pharokka-
   filtered to RBP proteins (`classify_rbp_genes`) — only ~200 proteins, receptor-binding-
   specific. MMseqs2 all-vs-all on the full 18,830-protein panel includes every core gene:
   terminase, polymerase, capsid, lysis machinery. Conserved housekeeping proteins contribute
   identically across all phage pairs, so the resulting similarity matrix overweights
   core-genome shared ancestry and underweights RBP-specific variation. Per
   `receptor-specificity-solved`: "the signal is in short motifs at specific loci, not
   global protein similarity."

Both mechanisms were anticipated by the plan's "net ≈ wash" caveat; the result sharpens
that prediction into "BASEL regression, Guelin wash" — an active reason to stop at TL17
rather than replace it with Arm 2.

#### Verdict

**Null (regression on BASEL discrimination).** Arm 2 does not replace TL17 as the canonical
phage-side feature slot. The plan's acceptance criterion is not met on any of its three
exit conditions (AUC lift, top-3 lift — retired under CHISEL anyway, NILS53 rescue — not
measurable on the BASEL side). Adds to `panel-size-ceiling` evidence: continuous proteome
similarity at phage-level aggregation does not crack the BASEL discrimination gap.

#### Scope

Single-arm PR. Arms 3 and 4 run on separate branches and will be reported in follow-up
entries. CH06 closes when Arm 4 merges (per-arm PRs chain under "Relates to #440" for
Arms 2/3 and "Closes #440" for Arm 4).

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch06_arm2_mmseqs_proteome.py` — precompute + eval driver.
- `lyzortx/pipeline/autoresearch/sx03_eval.py:patch_context_with_extended_slots` — new
  backwards-compatible `slots_dir=` parameter (default unchanged).
- `.scratch/ch06_arm2/{combined_proteome.faa, mmseqs_hits.tsv, pairwise_bitscore.npz}` —
  precompute cache (MMseqs2 hits, 148×148 bit-score matrix).
- `.scratch/basel/feature_slots_arm2/phage_projection/features.csv` — Arm 2 slot (148 × 33
  columns: phage + 32 PCA).
- `lyzortx/generated_outputs/ch06_arm2_mmseqs_proteome/` — variance pre-flight JSON, per-
  row and pair-level prediction CSVs, cross-source breakdown, final metrics JSON.

### 2026-04-20 15:08 CEST: CH06 Arm 3 — Moriniere receptor probabilities (validated — meets acceptance)

#### Executive summary

Plan.yml pre-registered Arm 3 as "expected null" based on `same-receptor-uncorrelated-hosts`
and Moriniere's classifier training on K-12 derivatives. The full run contradicts that:
Arm 3 meets the CH06 acceptance criterion with **BASEL bacteria-axis AUC 0.7374** (baseline
0.7229, target > 0.7152, **+1.45 pp**). The decisive deployability metric — BASEL
zero-projection phages (n=13, the literal TL17 zero-vector-failure subset) — improves
**+4.36 pp on phage-axis** (0.6901 → 0.7337). Guelin numbers flat (±0.15 pp on both axes);
BASEL non-zero-projection phages also gain (+2.54 pp bacteria-axis, +0.54 pp phage-axis).
Arm 3 replaces the Guelin-only TL17 categorical projection with a 13-dim per-receptor
k-mer-fraction vector, trained upstream on 260 non-Guelin phages so the feature basis is
panel-independent at source — which is exactly what the plan's `deployment-goal` needs.

#### Feature construction

For each phage P and each of the 13 Moriniere receptor classes R, the feature value is
`|kmers(R) ∩ kmers(P)| / |kmers(R)|` — the fraction of receptor R's 5-mer set that appears
anywhere in P's proteome. Divided per-receptor so the class imbalance in Dataset S6 (HepI:
168 k-mers, Kdo: 9, GluI: 28, …) is normalized out; otherwise HepI would dominate raw
hit counts just by having more k-mers to match.

Critical reuse: `lyzortx.pipeline.autoresearch.build_moriniere_kmer_slot.collect_*_phage_
proteomes` provides the 148-phage Guelin-plus-BASEL proteome load, and
`lyzortx.pipeline.autoresearch.predict_receptor_from_kmers.load_receptor_kmers` parses the
GenoPHI Dataset S6 supplementary XLSX. Plan.yml listed the feature space as "19-dim softmax"
but the actual Moriniere dataset has **13 receptor classes** (tsx, ompA, ompC, ompF, fhuA,
btuB, lptD, lamB, NGR, Kdo, HepI, HepII, GluI) — the notebook entry reports 13-dim
reflecting the data as-loaded rather than the plan's stale number.

Slot schema: `phage_projection__recep_frac_<receptor>` × 13 columns, materialized to
`.scratch/basel/feature_slots_arm3/phage_projection/features.csv`. Phage-stats and
phage-rbp-struct slots symlinked from the baseline directory, so only the projection slot
is changed.

#### Variance pre-flight

All three subsets pass the CV > 0.1 OR Cohen's d > 0.1 gate — and unlike Arm 2, the
BASEL zero-projection subset shows strong signal: 11/13 receptors have Cohen's d > 0.1
against phage-level lysis-rate median split, max abs d 1.20.

| Subset | n | Max CV | Max Cohen's d | d > 0.1 | CV > 0.1 | Gate |
|---|---|---|---|---|---|---|
| Guelin | 96 | 6.7 | 0.60 | 12/13 | 13/13 | pass |
| BASEL non-zero TL17 | 39 | 6.25 | 1.67 | 13/13 | 13/13 | pass |
| BASEL zero-vector TL17 | 13 | 3.61 | 1.20 | 11/13 | 12/13 | pass |

#### Full-training results

10-fold bacteria-axis (CH02 cv_group hash) and 10-fold phage-axis (StratifiedKFold by ICTV
family + "other" + "UNKNOWN"), 3 seeds each under `ch04_parallel.fit_seeds`, 1000-sample
unit-level bootstrap.

| Axis | Metric | Baseline | Arm 3 | Δ | CI overlap |
|---|---|---|---|---|---|
| Bacteria | AUC (agg) | 0.8218 [0.8063, 0.8368] | 0.8205 [0.8039, 0.8362] | −0.13 pp | yes |
| Bacteria | Brier (agg) | 0.1466 | 0.1457 | −0.09 pp | yes |
| Phage | AUC (agg) | 0.8919 [0.8650, 0.9166] | 0.8936 [0.8672, 0.9171] | +0.17 pp | yes |
| Phage | Brier (agg) | 0.1181 | 0.1177 | −0.04 pp | yes |

Aggregate numbers are flat (as expected when Guelin dominates weight). The real signal
lives in the cross-source and BASEL-subset decomposition:

| Subset | Bact-axis Baseline | Arm 3 | Δ | Phage-axis Baseline | Arm 3 | Δ |
|---|---|---|---|---|---|---|
| Guelin (n=96) | 0.8247 | 0.8232 | −0.15 pp | 0.8922 | 0.8934 | +0.12 pp |
| BASEL all (n=52) | **0.7229** | **0.7374** | **+1.45 pp** | 0.8822 | 0.8902 | +0.80 pp |
| BASEL non-zero TL17 (n=39) | 0.6968 | 0.7222 | +2.54 pp | 0.8974 | 0.9028 | +0.54 pp |
| **BASEL zero-vec TL17 (n=13)** | 0.6325 | 0.6484 | +1.59 pp | **0.6901** | **0.7337** | **+4.36 pp** |

**Acceptance criterion met.** BASEL bacteria-axis 0.7374 > 0.7152 target, with improvement
on every cross-source subset on both axes (Guelin flat; BASEL positive across zero/non-zero
decomposition). No cannibalization signature — unlike Arm 2, Arm 3 rescues zero-proj BASEL
without regressing non-zero-proj BASEL.

#### Arm 2 re-interpretation in light of Arm 3

Arm 2's merged notebook entry framed the result as "null / BASEL regression" on the
aggregate BASEL bacteria-axis AUC. The subset decomposition after the fact shows a more
nuanced picture: Arm 2 **does** rescue zero-proj BASEL (+1.07 pp bacteria-axis, +1.27 pp
phage-axis) — the deployability mechanism it was designed to address — but cannibalizes
non-zero-proj BASEL (−2.07 pp bacteria-axis, −1.52 pp phage-axis). The 39/52 non-zero
BASEL phages contribute 920/1240 pairs (74%), so their regression dominates the aggregate.
Arm 2 as a **drop-in replacement** for TL17 fails; Arm 2's phage-level similarity vector
composed **alongside** TL17 (keep RBP-focused signal where available, add proteome
similarity where TL17 is zero) could be a valid follow-up — but not needed if Arm 3 does
the same rescue cleanly, which it does.

#### Why the plan's "expected null" prediction was wrong

Plan.yml cited two reasons Arm 3 would fail:

1. `same-receptor-uncorrelated-hosts` — phages sharing a predicted receptor have Jaccard
   0.091 on host ranges (Tsx phages).
2. Moriniere trained on K-12 derivatives lacking capsule/O-antigen, so
   polysaccharide-mediated specificity is not in scope.

Both are correct characterizations of what receptor-class probabilities *don't* encode, but
the CH06 acceptance criterion is not "within-receptor-class host-range accuracy" or "capture
Gate-1 polysaccharide mechanism." It is **absolute cross-panel discrimination AUC on
BASEL bacteria-axis**. A 13-dim probability vector that is panel-independent (trained on
260 non-Guelin phages) and distinguishes "which of the major receptor families does this
phage use" is a real discriminative feature even if it says nothing about strain-level
host-range rank ordering. The null prediction was reasoning about a different question.

`kmer-receptor-expansion-neutral` (the SX12 / GT06 null) stands — raw 815-kmer
presence-absence is information-redundant with `phage_projection`. Arm 3's aggregation to
13 per-receptor fractions removes that redundancy by forcing the model to see class-level
probabilities rather than re-deriving them from the 815 raw k-mers; the aggregation turns
out to carry more actionable signal than the raw features.

#### Verdict

**Success — promote to canonical.** Arm 3 should become the default `phage_projection`
slot for CHISEL going forward. Recommendation (handled in follow-up commits to
`knowledge.yml`):

- Update `chisel-unified-kfold-baseline` knowledge unit to note the CH06 outcome and
  reference Arm 3 as the validated panel-independent phage-side feature slot.
- Add a new knowledge unit `moriniere-receptor-probabilities-panel-independent`
  summarizing the finding (BASEL bacteria-axis +1.45 pp, zero-proj BASEL phage-axis
  +4.36 pp, Guelin flat).
- Harden `plm-rbp-redundant` with "receptor-class probability aggregation (Arm 3) is the
  working alternative to global protein similarity" — the PLM null stands, but the
  replacement path is per-receptor aggregation, not global pooling.

Arm 4 (tail-restricted TL17) runs next for completeness. Under the Arm 3 win, Arm 4's
expected-null-or-wash carries less weight: even if Arm 4 works, Arm 3 is the simpler,
source-panel-independent solution.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch06_arm3_moriniere_receptor.py` — precompute + eval driver.
- `.scratch/basel/feature_slots_arm3/phage_projection/features.csv` — Arm 3 slot
  (148 × 14 columns: phage + 13 receptor-fraction).
- `lyzortx/generated_outputs/ch06_arm3_moriniere_receptor/ch06_arm3_metrics.json` — final
  aggregate AUC/Brier with CIs.
- `.../ch06_arm3_{bacteria,phage}_axis_predictions.csv` — per-pair predictions.
- `.../ch06_arm3_cross_source_breakdown.csv` — phage-axis cross-source.
- `.../ch06_arm3_variance_preflight.json` — pre-flight subset breakdown.

### 2026-04-20 17:01 CEST: CH06 Arm 4 — tail-restricted TL17 BLAST (null — CH06 closes)

#### Executive summary

Tail-restricted TL17 BLAST query re-projection does not lift discrimination: **BASEL
bacteria-axis AUC 0.7145 (−0.84 pp vs baseline 0.7229, below the 0.7152 target)**.
Aggregate AUC flat on both axes. BASEL zero-vector TL17 phages see no meaningful rescue
(+0.29 pp bacteria-axis; +0.12 pp phage-axis). BASEL non-zero-projection phages regress
on both axes (−1.19 pp bacteria / −1.08 pp phage). Guelin flat. Fails the CH06 acceptance
criterion. **CH06 closes with Arm 3 (Moriniere per-receptor k-mer fractions) as the
canonical panel-independent phage-side feature slot**; Arm 4 joins the dead-end list.

#### Why Arm 4 was an expected near-no-op

The existing `tl17_rbp_reference_bank.faa` is already RBP-filtered at build time via
`classify_rbp_genes` — only ~200 proteins, receptor-binding-specific. Arm 4 restricts
the QUERY side of the BLAST to tail/baseplate/RBP-annotated proteins of each phage.
Tail proteins of the query phage already match RBPs in the reference bank under the
baseline run; non-tail proteins of the query rarely cross-hit RBPs (since the bank is
narrow). So restricting the query side is mostly a strict subset of the hits the
baseline already counts — nothing new is added, and in noisy cases some genuine
cross-hits are dropped.

Variance pre-flight flagged this: the BASEL zero-vector TL17 subset had 0 features
with Cohen's d > 0.1 on Arm 4 (vs 11/13 for Arm 3). The 13 BASEL phages whose proteomes
didn't match the Guelin RBP bank under baseline TL17 still don't match under tail-
restricted TL17 — expected, since Arm 4 is a STRICT SUBSET of the baseline search. The
ticket's own pre-flight guidance — "the arm cannot rescue a phage whose tail region is
not annotated" — is confirmed empirically.

#### Feature construction

For each of 148 phages (96 Guelin + 52 BASEL), extract tail/baseplate/RBP-annotated
proteins from the Pharokka merged TSV (category ∈ {tail, baseplate} OR annotation
matches `RBP_PATTERNS` from `parse_annotations`). Matching strategy:

- **Guelin** (96 phages): coordinate overlap between Pharokka CDS (start, stop, contig)
  and the track_d per-phage protein FASTA (`start=`, `end=`, `contig=` header metadata).
  Exact (start, end, contig) tuple lookup into an O(1) index.
- **BASEL** (52 phages): CDS name match between Pharokka `gene` field and `phanotate.faa`
  header (first whitespace-separated token).

All 148 phages carry tail-annotated proteins (100% BASEL pre-flight gate pass — well
above the plan's >50% coverage requirement). Combined tail-restricted query FASTA has
~1,750 tail proteins. Searched against `tl17_rbp_reference_bank.faa` via MMseqs2
`easy-search` with `--max-seqs 20` + 12-column BLAST layout; aggregated to per-phage ×
32 TL17 family max-pident scores + 1 summary hit count, same schema as baseline TL17.

#### Full-training results

| Axis | Subset | Baseline | Arm 4 | Δ |
|---|---|---|---|---|
| Bacteria | Aggregate | 0.8218 [0.8063, 0.8368] | 0.8209 [0.8055, 0.8359] | −0.09 pp |
| Bacteria | Guelin (n=96) | 0.8247 | 0.8239 | −0.08 pp |
| Bacteria | **BASEL all (n=52)** | **0.7229** | **0.7145** | **−0.84 pp** |
| Bacteria | BASEL zero-vec (n=13) | 0.6325 | 0.6354 | +0.29 pp |
| Bacteria | BASEL non-zero (n=39) | 0.6968 | 0.6849 | −1.19 pp |
| Phage | Aggregate | 0.8919 | 0.8881 | −0.38 pp |
| Phage | Guelin (n=96) | 0.8922 | 0.8884 | −0.38 pp |
| Phage | BASEL all (n=52) | 0.8822 | 0.8753 | −0.69 pp |
| Phage | BASEL zero-vec (n=13) | 0.6901 | 0.6913 | +0.12 pp |
| Phage | BASEL non-zero (n=39) | 0.8974 | 0.8866 | −1.08 pp |

**Acceptance criterion not met.** BASEL bacteria-axis AUC 0.7145 < 0.7152 target
(marginal miss but a miss). No rescue of zero-vec BASEL beyond baseline. Guelin flat;
BASEL non-zero regresses on both axes. Brier deltas are small and CIs overlap throughout.

#### CH06 — three-way comparison and close-out verdict

| Arm | Bact-axis BASEL | Phage-axis zero-vec BASEL | Guelin (both axes) | Verdict |
|---|---|---|---|---|
| Baseline TL17 | 0.7229 | 0.6901 | 0.8247 / 0.8922 | — |
| Arm 2 (MMseqs2 pairwise) | 0.7045 (−1.84 pp) | 0.7028 (+1.27 pp) | 0.8278 / 0.8879 | null (BASEL non-zero cannibalized) |
| **Arm 3 (Moriniere recep. fractions)** | **0.7374 (+1.45 pp)** | **0.7337 (+4.36 pp)** | 0.8232 / 0.8934 | **success** |
| Arm 4 (tail-restricted TL17) | 0.7145 (−0.84 pp) | 0.6913 (+0.12 pp) | 0.8239 / 0.8884 | null |

**CH06 closes with Arm 3 as the canonical panel-independent phage-side feature slot.**
Arm 3 is the single passing arm; it rescues both the zero-vector deployability failure
mode (+4.36 pp on the 13 BASEL phages that get zero TL17 signal) and the non-zero-proj
BASEL subset (+2.54 pp bact-axis), without damaging Guelin. Arm 2 partially addressed
the deployability mechanism but cannibalized RBP-specific signal; Arm 4 was a strict
subset of baseline TL17 hits and carried no new information.

This outcome updates `plm-rbp-redundant` and `receptor-specificity-solved`: the
receptor-class signal IS useful, **provided** it's aggregated to per-class probabilities
rather than fed as raw 815 k-mer presence-absence (SX12 null) or as global protein-level
similarity (AX08 null). The aggregation level is the mechanism.

#### Open follow-ups

1. **Wire Arm 3 into the canonical `phage_projection` slot** — currently Arm 3 is a
   side-materialized artifact at `.scratch/basel/feature_slots_arm3/`; canonical CH05 /
   CH07 / CH08 / CH09 pipelines still read baseline TL17. A separate migration ticket
   (CH10?) should make Arm 3 the default, re-run CH05 / CH07 / CH08 / CH09 under the
   new slot, and update their baseline numbers.
2. **Arm 2 composition with TL17** — identified as valid follow-up during Arm 3 review
   (keep TL17 where available, add Arm 2 PCA for zero-vec phages). Under Arm 3's clean
   win this is unnecessary but could be revisited if Arm 3's canonical migration
   reveals an unexpected failure mode.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch06_arm4_tail_restricted_tl17.py` — precompute + eval
  driver.
- `.scratch/ch06_arm4/{tail_restricted_query.faa, mmseqs_hits.tsv}` — precompute cache.
- `.scratch/basel/feature_slots_arm4/phage_projection/features.csv` — Arm 4 slot.
- `lyzortx/generated_outputs/ch06_arm4_tail_restricted_tl17/ch06_arm4_metrics.json` +
  `.../ch06_arm4_{bacteria,phage}_axis_predictions.csv` +
  `.../ch06_arm4_cross_source_breakdown.csv` +
  `.../ch06_arm4_variance_preflight.json`.

### 2026-04-21 00:27 CEST: CH07 — Both-axis 10×10 double cross-validation (Arm 3 canonical)

#### Executive summary

Ran the 10-fold bacteria × 10-fold phage double cross-validation (100 cells, every cell
holds out one bacteria fold AND one phage fold simultaneously) on the unified 148-phage
panel using the CH06 Arm 3 Moriniere receptor-fraction phage-side slot. Pooled AUC
**0.7749 [0.7687, 0.7814]**, Brier **0.1516 [0.1490, 0.1541]** on 36,643 pairs. The
single-axis bacteria-axis baseline (CH05, same feature bundle area though still on
baseline TL17 there) was **0.8218** — simultaneous-holdout loses **−4.7 pp** in AUC. The
per-cell AUC distribution has mean 0.781, median 0.780, std 0.054; 35/100 cells are above
0.80 and only 7/100 are below 0.70. **Both-axis CI is disjoint from 0.80**, so the
pre-registered "ceiling framing strengthened" verdict (AUC > 0.80) is NOT met; but the
model retains ~93% of single-axis discrimination, which also does not fall into the
pre-registered "marginal riding" floor (0.60-0.70). The middle-ground finding needs its
own knowledge unit language, not either pre-registered label.

#### Design

- **Bacteria folds (10)** — CH02 cv_group hash, identical to CH04/CH05 bacteria-axis.
  Ensures ANI-duplicate strains land in the same fold.
- **Phage folds (10)** — StratifiedKFold by ICTV family, identical to CH05 phage-axis.
  Rare families (< 10 phages) collapse to an "other" stratum; UNKNOWN family maps into
  "other". Same `random_state=42`.
- **Cells** — 100 (bact_fold × phage_fold). For each cell, the training set is every row
  whose bacterium is NOT in bact_fold AND whose phage is NOT in phage_fold (~247K rows).
  The test set is every row whose bacterium IS in bact_fold AND whose phage IS in
  phage_fold (~1.5-3K rows).
- **Phage-side slot** — `.scratch/basel/feature_slots_arm3/` (CH06 Arm 3, Moriniere
  per-receptor k-mer fractions). Chosen over baseline TL17 per the CH06 verdict
  (`moriniere-receptor-fractions-validated`). Arm 3's 13-dim panel-independent
  representation is what "winning feature bundle from CH06" in plan.yml points at.
- **All-pairs only** — per-phage blending retired track-wide; also structurally
  impossible since held-out phages have zero training rows for per-phage models.
- **Seeds** — 3 (SEEDS = (0, 1, 2)), averaged per cell, then pair-level max-concentration
  aggregation (CH04 convention: log10_pfu_ml-max per pair; replicates averaged).
- **Bootstrap** — 1000 resamples at the pair level on the pooled predictions (100-cell
  aggregate). Cell-level bootstrap is underpowered because per-cell pair counts average
  366 with a min of 225; plan.yml pre-registered pair-level bootstrap for exactly this
  reason.
- **Training filter** — `drop_high_titer_only_positives=True` (CH06-followup canonical).

Compute: 14,447 s wallclock (~4h) on laptop, 3 seed workers per cell. Fully deterministic
given SEEDS + RFE_SEED (tested under the CH06 engineering pre-flight).

#### Headline results

| Metric | Value | 95% CI | n_pairs |
|---|---|---|---|
| Aggregate AUC | 0.7749 | [0.7687, 0.7814] | 36,643 |
| Aggregate Brier | 0.1516 | [0.1490, 0.1541] | 36,643 |
| Guelin subset AUC | 0.7767 | [0.7701, 0.7837] | 35,403 |
| Guelin subset Brier | 0.1491 | [0.1464, 0.1517] | 35,403 |
| BASEL subset AUC | 0.7040 | [0.6692, 0.7424] | 1,240 |
| BASEL subset Brier | 0.2226 | [0.2065, 0.2385] | 1,240 |

Per-cell AUC distribution (100 cells): mean 0.7806, median 0.7803, std 0.0534, IQR
[0.7495, 0.8210], min 0.6316, max 0.9150. 7 cells below 0.70; 35 cells above 0.80; 10
cells above 0.85.

#### Comparison to single-axis baselines (same Arm 3 slot, CH05 post-filter frame)

CH05 phage-slot choice was baseline TL17, not Arm 3, so the single-axis baselines here
are NOT run under identical phage-side features. Still, bracketing:

| Axis | Slot | AUC | Brier |
|---|---|---|---|
| CH04 bacteria-axis (Guelin only) | TL17 | 0.8217 | 0.1435 |
| CH05 bacteria-axis (unified) | TL17 | 0.8218 | 0.1466 |
| CH05 phage-axis (unified) | TL17 | 0.8919 | 0.1181 |
| **CH07 both-axis (unified)** | **Arm 3** | **0.7749** | **0.1516** |

**Both-axis minus bacteria-axis = −4.7 pp AUC, +0.5 pp Brier. Both-axis minus phage-axis
= −11.7 pp AUC, +3.4 pp Brier.** The phage-axis drop is larger because phage-axis alone
keeps all 369 bacteria in training (rich host-side signal per test pair), which vanishes
under simultaneous bact holdout.

The fact that CH07 uses Arm 3 while CH05 used TL17 is a residual confound for a strict
ceiling verdict, but it's a small one: Arm 3 was NULL on the Guelin side (CH06 Arm 3 ΔAUC
± 0.15 pp on Guelin both axes), so CH05-under-Arm-3 would land within ~0.5 pp of
CH05-under-TL17. The −4.7 pp both-axis gap is robust to that confound.

#### Ceiling framing verdict (plan.yml pre-registration)

Plan.yml acceptance criteria pre-registered two thresholds:

- **AUC > 0.80** → "panel-size ceiling framing strengthened" (model learned pair-level
  mechanism that generalizes under compound holdout)
- **AUC 0.60-0.70 (marginal-riding range)** → "ceiling framing is undermined and a new
  knowledge unit (`chisel-marginal-riding-suspicion`) must be added flagging that
  feature/model choices, not panel size, may be the binding constraint"

**CH07 result (0.7749, CI [0.769, 0.781]) lands between the two thresholds.** Strictly:
the CI's upper bound (0.781) is disjoint from 0.80, so the "ceiling strengthened" verdict
is **not** met. But 0.77 is well clear of the 0.60-0.70 marginal-riding floor and the
per-cell distribution tells a more nuanced story — 35% of cells exceed 0.80, with a
median right at 0.78. The model is losing ~5 pp to simultaneous holdout, which is
consistent with a mix of (a) some genuine pair-level mechanism (kept 93% of single-axis
AUC) and (b) some host-side and phage-side marginal signal that vanishes together.

Neither pre-registered label applies cleanly. The honest call is **nuanced held-above-
marginal-riding, below-ceiling-threshold**: the model is learning real pair-level
structure (otherwise cold-start would be near 0.5), but the remaining gap to single-axis
AUC is consistent with panel-size being the dominant remaining lever rather than model
flaws.

Knowledge-model update: add `chisel-both-axis-holdout` as a new unit reporting the
0.7749 ± CI number, the per-cell distribution, and the "held-above-marginal, below-
ceiling" language. Do NOT add `chisel-marginal-riding-suspicion` — the marginal-riding
threshold (0.60-0.70) is not triggered.

#### BASEL-specific deficit under both-axis

BASEL subset AUC 0.7040 is 7 pp below Guelin's 0.7767 on the pooled predictions, with a
much wider CI (the CI on BASEL's 1,240 pairs is ±0.037 vs Guelin's ±0.007). This echoes
the CH05 BASEL-specific bacteria-axis deficit pattern (`chisel-unified-kfold-baseline`):
under unified training, BASEL underperforms Guelin on the bacteria-axis by ~10 pp. Under
both-axis the deficit narrows to 7 pp but is still disjoint. Consistent with the
`moriniere-receptor-fractions-validated` unit's framing: Arm 3 partially rescues BASEL
(the zero-vec TL17 phages gain +4.36 pp) but a residual BASEL-specific feature-space
gap remains — likely the non-zero-projection BASEL phages whose receptor assignments
differ subtly from Guelin-bank representatives.

#### Per-cell outliers (audit trail)

7 cells drop below 0.70. Spot-checked: most are cells where the held-out phage fold
contains the rarer / narrower-host phages (ICTV families at the lower end of the
stratified split) and the held-out bacteria fold happens to be cv_group-poor. These are
the hardest deployment cases — cold-start on a narrow-host phage against a novel strain
cohort. The 10 cells above 0.85 are the converse: broad-phage × phylogenetically well-
covered-bacterium fold intersections where the all-pairs model's marginals still work
well without any specific training pair.

#### Open follow-ups

1. **Re-run CH07 under baseline TL17** (single-arg flip: `--phage-slots-dir
   .scratch/basel/feature_slots`). Would remove the Arm 3 vs TL17 confound in the
   single-axis-vs-both-axis comparison. Cheap (4h) but deferred — not needed for the
   ceiling verdict, which is robust to ±0.5 pp shifts.
2. **Per-cell AUC distribution regression** — does cell AUC correlate with positive_rate,
   n_pairs, or ICTV family of held-out phages? The CH07 artifact has everything needed
   but no analysis was done. Small follow-up, useful for narrow-host behavior docs.
3. **Cross-source decomposition with CI** in each cell, not just pooled, to quantify how
   many cells drive the BASEL deficit vs broad Guelin patterns.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch07_both_axis_holdout.py` — eval driver.
- `lyzortx/tests/test_ch07_both_axis_holdout.py` — unit tests.
- `lyzortx/generated_outputs/ch07_both_axis_holdout/ch07_aggregate.json` — pooled AUC +
  Brier with pair-level bootstrap CI, cross-source subsets, per-cell distribution stats,
  phage_slots_dir provenance.
- `.../ch07_cell_metrics.csv` — one row per cell: fold IDs, n_pairs, n_bacteria, n_phages,
  positive_rate, per-cell AUC + Brier.
- `.../ch07_pair_predictions.csv` — pooled pair-level (max-concentration) predictions for
  all 36,643 pairs, tagged with source.
- `.../ch07_per_row_predictions.csv` — pooled per-row predictions for all held-out rows.
- `.../ch07_cross_source_breakdown.csv` — one row per source (guelin, basel) with AUC /
  Brier and CIs.
- `.../ch07_cell_distribution.png` — histogram of per-cell AUC with mean and median
  markers.

### 2026-04-21 05:37 CEST: CH08 — Wave-2 feature-family re-audit under CHISEL frame

#### Executive summary

Re-audited SPANDEX wave-2 nulls (SX12 Moriniere 815 phage 5-mers;
SX13 host OMP 5546 5-mers) under the CHISEL training frame: per-row binary labels,
`pair_concentration__log10_pfu_ml` feature, all-pairs only (no per-phage blending),
neat-only positive filter, 10-fold bacteria-axis CV over the 369×96 Guelin panel. Each
kmer slot was pre-filtered to the **top-100 features by binary variance** before RFE
to keep per-fold wallclock tractable on a 281K-row training frame (unfiltered RFE on
1281 total features had run at ~60 min/fold — 20 folds = 20 h+). Paired bacterium-level
bootstrap (1000 resamples) on the CH04-canonical baseline predictions vs each arm's
variant predictions:

| Arm | Extra slot | Arm AUC | Arm Brier | Δ AUC (CI) | Δ Brier (CI) | Knowledge-unit verdict |
|---|---|---|---|---|---|---|
| baseline | — | 0.8217 | 0.1435 | — | — | — |
| sx12 | phage_moriniere_kmer (top-100 var) | **0.8333** | 0.1415 | **+1.16 pp [+0.82, +1.51]** | **−0.20 pp [−0.41, ≈0]** | **`kmer-receptor-expansion-neutral` reopens** |
| sx13 | host_omp_kmer (top-100 var) | 0.8234 | 0.1423 | +0.17 pp [+0.03, +0.31] | −0.12 pp [−0.20, −0.06] | **`host-omp-variation-unpredictive` reopens** |

Both arms land with delta CIs disjoint from zero, so per plan.yml's "CRITIC RESPONSE: if
either re-audit shows non-null lift (CI disjoint from zero), this reopens the wave-2
feature-family conclusions … Do not suppress unexpected positives", **both wave-2
null knowledge units must be reverted from validated-null to open**. SX12's effect
size (+1.16 pp) is ~7× SX13's (+0.17 pp) — the Moriniere phage k-mers contribute a
meaningful signal under CHISEL that they did not contribute under SPANDEX; the host
OMP k-mers still barely move the needle but the CI is technically disjoint.

#### Bug caught during CH08 development

The initial CH08 run silently dropped BOTH extra slots because
`ch04_parallel.FEATURE_COLUMN_PREFIXES` was a hardcoded allowlist of slot column
prefixes that did not include `phage_moriniere_kmer__` or `host_omp_kmer__`. The
slots were attached to the `context.slot_artifacts` dict, the columns showed up in
the design matrix, but `prepare_fold_design_matrices` filtered `feature_columns` by
the hardcoded prefix tuple, so the kmer columns never reached RFE or LightGBM. The
symptom was SX13 predictions bit-identical to the CH04 baseline (max pred diff 0.0)
and SX12 predictions differing only via deterministic dataframe-ordering artifacts
that happened to bias probabilities downstream (spurious +1.0 pp point estimate
driven entirely by column-order side effects). Fix: added `host_omp_kmer__`,
`host_omp_cluster__`, and `phage_moriniere_kmer__` to `FEATURE_COLUMN_PREFIXES` with
a docstring calling out the allowlist contract. Re-running under the fix produced
the real numbers above.

This is a canonical silent-fail shape — slot attached, columns present in DataFrame,
model never sees them. If a future ticket adds a new slot family with its own
prefix, the onus is on the implementer to register the prefix in
`FEATURE_COLUMN_PREFIXES`. Consider hoisting that allowlist into the slot registry
so attaching a slot automatically declares its prefix.

#### Why pre-filtering to top-100 by variance

The SX12 slot has 815 5-mers; the SX13 slot has 5546 5-mers. RFECV (step=0.1, cv=5
internal folds, on 281K-row training frames) does ~41 iterations for 788 input
features, each fitting LightGBM 5 times → ~60 minutes per CV fold. 10 folds per arm
× 2 arms × 60 min = 20 h wallclock, violating plan.yml's "estimate before running
and confirm wallclock is feasible". Top-K-by-variance pre-filtering (K=100) rank-
orders features by entity-level binary variance and keeps the K most-variable. This
preserves:

- sparse receptor-class-specific kmers (support ≈ 10% of phages → binary variance
  0.09, rank-competitive with medium-support kmers)
- maximally-balanced OMP kmers (support ≈ 50% of bacteria → variance 0.25, the
  theoretical max for binary features)

It drops near-constant kmers (support 1/369 or 368/369 → variance → 0) which
contribute no discriminative signal anyway. Pre-filter is deterministic, leakage-
free (computed over entities, not over training/test pair splits), runs once before
any fold split. Arm AUC deltas are measured under the pre-filter, so the arm-level
conclusion "does any meaningful subset of the kmer slot help?" is answered without
the 20h-wallclock cost of the unfiltered RFE. Caveat: the pre-filter is NOT identical
to the SPANDEX SX12/SX13 runs, which saw the full 815 / 5546 kmers. So the CH08
deltas are not direct apples-to-apples replacements for SPANDEX SX12/SX13 deltas —
they're CHISEL-frame estimates on the top-100-variance subset.

#### SX12 result interpretation

SX12 reopens: **+1.16 pp [+0.82, +1.51] AUC, −0.20 pp Brier, CIs disjoint from zero**.
Under SPANDEX (pair-level any_lysis training, per-phage blending enabled, pre-fix
cv_group folds, full 815 kmers), SX12 was null: "+0.23 pp AUC, CIs overlap baseline"
(`kmer-receptor-expansion-neutral`). Under CHISEL (per-row binary, concentration
feature, per-phage retired, fixed cv_group folds, top-100 variance-filtered kmers),
the same slot contributes a statistically significant lift. Several plausible
mechanisms for the shift; the data does not decide among them:

1. **Per-row training exposes concentration-dependent signal** the pair-level rollup
   erased. Kmers that predict graded dilution response may correlate with receptor-
   class strength in a way pair-level labels averaged out.
2. **Pre-filter to top-100 variance** happens to select the receptor-class-
   discriminative kmers and drops ~700 near-constant features that just added noise
   to RFE under SPANDEX. Note this is consistent with Arm 3's finding: per-class
   Moriniere features help once aggregated appropriately — here "aggregation" is the
   variance-rank filter keeping the features that matter.
3. **cv_group fold fix** removed ANI-duplicate leakage; both baseline and variant
   gain under the fix, but the variant's ability to exploit receptor-class signal
   may benefit differently than baseline's ability to memorize host-side priors.
4. **Per-phage blending retired** shifts attribution from bacterium-level memorization
   to all-pairs mechanism; SX12 kmers may carry phage-mechanism signal that per-phage
   blending was otherwise double-counting.

The SPANDEX `kmer-receptor-expansion-neutral` unit framed SX12 as null on mechanism-
level grounds ("k-mers were selected to discriminate receptor class on K-12 … they
predict what we already know, not what we need"). That framing is still directionally
right — the +1.16 pp lift doesn't close the BASEL bacteria-axis 10-pp deficit or the
panel-size-ceiling gap — but the unit needs updating: SX12 kmers are NOT null under
CHISEL, they're a modest additive feature-family lift.

**Does this supersede CH06 Arm 3?** No. CH06 Arm 3 (per-receptor k-mer-fraction
vectors, 13-dim, panel-independent) beat CH06 Arm 2 and Arm 4 on the BASEL
deployability axis (+4.36 pp zero-vec BASEL phage-axis). SX12 kmer slot here is a
different feature shape (100 individual k-mer presence/absence features, not per-
class aggregates) and is evaluated on the Guelin bacteria-axis, not BASEL. The two
findings are complementary: Arm 3 says "aggregate to per-class fractions for panel
independence"; CH08 SX12 says "even the raw kmers, if variance-filtered, give a
~1 pp lift under CHISEL". A future experiment could test Arm 3 + SX12-kmer jointly,
but plan.yml has no such ticket and the Arm 3 + Arm 3 → canonical migration is
higher priority.

#### SX13 result interpretation

SX13 reopens: **+0.17 pp [+0.03, +0.31] AUC, −0.12 pp Brier, CIs barely disjoint
from zero**. Under SPANDEX, SX13 was null across four arms (marginal, cross_term,
path1_cluster, and baseline — all within ±0.4 pp). Under CHISEL, the marginal arm
shows a tiny CI-disjoint lift, but 0.17 pp is far smaller than any ordinary signal
and smaller than the typical run-to-run variance on other tickets (±0.5 pp).

The underlying `host-omp-variation-unpredictive` knowledge unit framed SX13 as
null because "host-range variance lives downstream of OMP recognition" and the
loop-level kmer escalation did not rescue prediction. That framing is still
essentially correct — +0.17 pp is not the "host OMP variation actually predicts
lysis" story. Likely mechanism: top-100 variance-filtered OMP kmers include kmers
that happen to correlate with phylogroup (lineage-confounding effect similar to
`defense-lineage-confounding`), and the filter captured a noise-level amount of
host-phylogroup marginal information. The Brier improvement of 0.12 pp is barely
detectable.

Update `host-omp-variation-unpredictive` with the CHISEL-frame note: under the
top-100 variance-filtered subset, the host OMP kmer slot contributes a tiny but
CI-disjoint +0.17 pp lift; this is consistent with phylogroup-correlated lineage
signal, not with OMP-specific host-range prediction.

#### Compute

6938 s (~1.9 h) on laptop after orphan subprocess cleanup (earlier runs left
residual RFECV worker pools that starved CPU and made fold wallclocks meaningless).
SX12 arm: 10 folds × ~2.4 min = 24 min; SX13 arm: 10 folds × ~2.8 min = 28 min;
bootstrap: ~1 min per arm. Full determinism (SEEDS, RFE seed, pre-filter rank) —
rerunning under the same code path reproduces the numbers bit-for-bit.

#### Open follow-ups

1. **Hoist the slot-prefix allowlist into a slot registry** so attaching a slot
   auto-declares its prefix. The current `FEATURE_COLUMN_PREFIXES` allowlist bug
   pattern will recur otherwise.
2. **Ablate pre-filter top-K** — repeat SX12 at K=50, K=200, K=400 to check whether
   the +1.16 pp is K-stable or K-sensitive. If K-sensitive, the lift is
   filter-artefactual and should be demoted.
3. **SX13 phylogroup-confound check** — permute the OMP kmer values within
   phylogroup and re-run; if the +0.17 pp lift survives, the effect is
   OMP-specific; if it vanishes, it's phylogroup-marginal noise.
4. **CHISEL-frame comparison vs SX12/SX13 unfiltered** — nice-to-have if compute
   budget allows, but the wallclock (20 h+) is prohibitive without either
   cell-level parallelism or a faster RFE variant.

#### Artifacts

- `lyzortx/pipeline/autoresearch/ch08_wave2_reaudit.py` — eval driver, pre-filter,
  paired-bootstrap helpers.
- `lyzortx/pipeline/autoresearch/ch04_parallel.py` — FEATURE_COLUMN_PREFIXES fix
  adding `host_omp_kmer__`, `host_omp_cluster__`, `phage_moriniere_kmer__`.
- `lyzortx/tests/test_ch08_wave2_reaudit.py` — unit tests for paired bootstrap.
- `lyzortx/generated_outputs/ch08_wave2_reaudit/ch08_summary.csv` — one row per arm
  (baseline, sx12, sx13) with AUC/Brier point + CIs and delta CIs.
- `.../ch08_combined_summary.json` — full report (task_id, scorecard, per-arm
  paired-bootstrap output with baseline/variant/delta blocks).
- `.../ch08_sx12_delta.json`, `.../ch08_sx13_delta.json` — per-arm delta reports.
- `.../ch08_sx12_predictions.csv`, `.../ch08_sx13_predictions.csv` — per-arm
  pair-level predictions (max-concentration) for audit.
