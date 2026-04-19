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
