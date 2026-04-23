### 2026-02-15: ST0.1 diagnostics, CI regression gate, and ST0.1b decision

#### What we implemented

1. Implemented ST0.1 label-policy step: `lyzortx/pipeline/steel_thread_v0/steps/st01_label_policy.py`.
2. Added ST0.1 regression checker and baseline: `lyzortx/pipeline/steel_thread_v0/checks/check_st01_regression.py`,
   `lyzortx/pipeline/steel_thread_v0/baselines/st01_expected_metrics.json`.
3. Added CI workflow to run the regression gate on push and pull request:
   `.github/workflows/steel-thread-regression.yml`.

#### ST0.1 findings from current internal data

- Raw rows: `318,816`.
- Observed bacteria-phage pairs: `35,424` (full 369 x 96 grid).
- Hard labels:
  - Positive: `9,720` (`27.44%`).
  - Negative: `25,546` (`72.12%`).
  - Unresolved: `158` (`0.45%`).
- Hard-label coverage: `99.554%`.
- Uncertainty flags:
  - Conflicting interpretable observations: `8,917` pairs (`25.17%`).
  - Has uninterpretable (`score='n'`): `4,537` pairs (`12.81%`).
  - High uninterpretable fraction: `1,444` pairs (`4.08%`).

#### Conflict interpretation summary

- Not all conflicts imply pure replicate noise.
- `2,359` conflicting pairs (`6.66%` of all pairs) are clean cross-dilution shifts.
- `6,558` conflicting pairs (`18.51%` of all pairs) include within-dilution replicate disagreement.
- Current positive rule (`any 1`) captures weak positives; `23.21%` of positives are single-hit (`1/9`).

#### Decision and plan impact

We will add **ST0.1b** as a parallel stricter label view with confidence tiers: `high_conf_pos`, `high_conf_neg`, and
`ambiguous`.

Why:

1. Preserve high recall behavior from ST0.1 for candidate generation.
2. Add a higher-trust slice for model debugging and honest early benchmarking.
3. Report metrics on both full-label and high-confidence slices to avoid noise-driven false gains.

### 2026-02-15: ST0.1b implemented (strict confidence tiers)

#### What we implemented in ST0.1b

1. Added ST0.1b step: `lyzortx/pipeline/steel_thread_v0/steps/st01b_confidence_tiers.py`.
2. Added ST0.1b regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st01b_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st01b_expected_metrics.json`.
4. Extended CI workflow to run both ST0.1 and ST0.1b regression gates.

#### ST0.1b policy (v1)

- `high_conf_pos`: hard label is positive, `score_1_count >= 2`, positive fraction among interpretable observations
  `>= 0.4`, and `score_n_count <= 1`.
- `high_conf_neg`: hard label is negative, `score_0_count >= 7`, and `score_n_count <= 1`.
- `ambiguous`: all remaining pairs.

#### ST0.1b output summary on current internal data

- Total pairs: `35,424`.
- `high_conf_pos`: `4,135`.
- `high_conf_neg`: `24,203`.
- `ambiguous`: `7,086`.
- Strict-slice coverage: `0.799966`.
- Strict positive fraction within strict slice: `0.145917`.

#### Interpretation

1. ST0.1b provides a narrower, more conservative training/evaluation slice while preserving broad coverage (~80%).
2. The strict slice remains class-imbalanced, so downstream modeling will need imbalance-aware training.
3. Positive-side conflict burden remains non-trivial in this first strict policy and should be stress-tested in ST0.2.

### 2026-02-15: ST0.2 implemented (canonical pair table)

#### What we implemented in ST0.2

1. Added ST0.2 step: `lyzortx/pipeline/steel_thread_v0/steps/st02_build_pair_table.py`.
2. Added ST0.2 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st02_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st02_expected_metrics.json`.
4. Extended CI workflow to run ST0.2 regression in addition to ST0.1 and ST0.1b.

#### ST0.2 output summary on current internal data

- Output rows: `35,424` (369 bacteria x 96 phages).
- Output schema: `64` columns in `st02_pair_table.csv`.
- Strict-slice rows: `28,338` (`0.799966` of all rows), inherited from ST0.1b.
- Join coverage:
  - Host metadata missing: `0`
  - Phage metadata missing: `0`
  - CV group missing: `0`
  - Interaction-matrix missing: `156` (auxiliary only; non-blocking).

#### ST0.2 Notes

1. ST0.2 treats `interaction_matrix.csv` as an auxiliary reference and explicitly marks it as non-feature to avoid
   leakage confusion.
2. A notable host metadata gap remains in `host_abc_serotype` (`22,848` missing row-values in the pair table).
3. ST0.2 is now stable and regression-gated; ST0.3 can consume `st02_pair_table.csv` as canonical input.

### 2026-02-15: ST0.3 implemented (leakage-safe split protocol)

#### What we implemented in ST0.3

1. Added ST0.3 step: `lyzortx/pipeline/steel_thread_v0/steps/st03_build_splits.py`.
2. Added ST0.3 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st03_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st03_expected_metrics.json`.
4. Extended CI workflow to run ST0.3 regression in addition to ST0.1 through ST0.2.

#### ST0.3 split summary on current internal data

- Rows assigned: `35,424`.
- Group key: `cv_group` from ST0.2.
- Fixed holdout groups: `57 / 283` (`0.201413`).
- Holdout rows: `6,240`.
- Non-holdout rows: `29,184`.
- CV fold rows on non-holdout:
  - fold 0: `5,760`
  - fold 1: `5,280`
  - fold 2: `6,624`
  - fold 3: `5,184`
  - fold 4: `6,336`

#### Leakage check results

- Holdout bacteria overlap count: `0`.
- Holdout cv_group overlap count: `0`.
- Cross-fold cv_group overlap count: `0`.

#### ST0.3 Notes

1. ST0.3 defines the locked v0 split contract for downstream model training and evaluation.
2. Assignment is deterministic via hash with salt `steel_thread_v0_st03_split_v1`.
3. ST0.4 should consume `st03_split_assignments.csv` directly and avoid custom split logic.

#### ST0.3 Interpretation

1. The leakage checks are clean (`0` overlap across holdout/train and across CV folds by `cv_group`), so this split
   protocol is suitable as the v0 benchmark contract.
2. Fold sizes are not perfectly balanced, but they are close enough for v0 model comparison; metrics should still be
   reported per fold and macro-averaged to reduce sensitivity to fold-size differences.
3. The strict trainable subset remains dominated by negatives in every fold, so ST0.4 should use class-imbalance-aware
   training and report both ranking metrics and calibration metrics.
4. Holdout size (~17.6% of rows) is large enough to be meaningful for a fixed benchmark while preserving enough
   non-holdout data for model development.

### 2026-02-15: ST0.4 implemented (baseline training)

#### What we implemented in ST0.4

1. Added ST0.4 step: `lyzortx/pipeline/steel_thread_v0/steps/st04_train_baselines.py`.
2. Added ST0.4 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st04_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st04_expected_metrics.json`.
4. Extended CI workflow to run ST0.4 regression and install `scikit-learn` explicitly.

#### ST0.4 output summary on current internal data

- Train rows (non-holdout hard-labeled): `29,031`.
- Holdout eval rows (hard-labeled): `6,235`.
- Vectorized feature count: `425`.
- Comparator model: `DummyClassifier(strategy='prior')`.
- Strong baseline: `LogisticRegression(class_weight='balanced', solver='liblinear')`.

#### Holdout metrics

- Dummy baseline:
  - Brier: `0.189304`
  - Log loss: `0.566574`
  - ROC-AUC: `0.500000`
  - Top-3 hit rate (all strains): `0.015385`
- Logistic baseline:
  - Brier: `0.171223`
  - Log loss: `0.521944`
  - ROC-AUC: `0.826948`
  - Top-3 hit rate (all strains): `0.846154`

#### ST0.4 Interpretation

1. The strong baseline materially outperforms the comparator on every tracked holdout metric, so ST0.4 clears the
   minimum "better than naive" bar for steel-thread viability.
2. Top-3 hit rate is substantially below the Tier 1 benchmark target; this is expected at this stage and motivates ST0.5
   calibration plus ST0.6 recommendation logic.
3. The feature-space and model artifact files are now stable inputs for downstream calibration/ranking.

### 2026-02-15: ST0.5 implemented (calibration and ranking)

#### What we implemented in ST0.5

1. Added ST0.5 step: `lyzortx/pipeline/steel_thread_v0/steps/st05_calibrate_rank.py`.
2. Added ST0.5 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st05_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st05_expected_metrics.json`.
4. Extended CI workflow to run ST0.5 regression in addition to ST0.1 through ST0.4.

#### ST0.5 output summary on current internal data

- Calibration rows (fold 0, non-holdout hard-labeled): `5,755`.
- Holdout eval rows (hard-labeled): `6,235`.
- Methods implemented per model: raw, isotonic calibration, and Platt scaling.
- Output artifacts:
  - `st05_calibration_summary.csv`
  - `st05_pair_predictions_calibrated.csv`
  - `st05_ranked_predictions.csv`

#### Holdout calibration metrics (logreg model)

- Raw: Brier `0.171223`, Log loss `0.521944`, ECE `0.176341`.
- Isotonic: Brier `0.140218`, Log loss `0.500302`, ECE `0.031802`.
- Platt: Brier `0.137795`, Log loss `0.430845`, ECE `0.029253`.

#### ST0.5 Interpretation

1. Calibration materially improved probabilistic quality, especially ECE, which dropped by an order of magnitude
   relative to raw logreg outputs.
2. Platt scaling slightly outperformed isotonic on holdout in this split configuration.
3. ST0.6 should use calibrated ranking scores (not raw model scores) for top-3 recommendation generation.

### 2026-02-15: ST0.6 implemented (top-3 recommendation generation)

#### What we implemented in ST0.6

1. Added ST0.6 step: `lyzortx/pipeline/steel_thread_v0/steps/st06_recommend_top3.py`.
2. Added ST0.6 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st06_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st06_expected_metrics.json`.
4. Extended CI workflow to run ST0.6 regression in addition to ST0.1 through ST0.5.

#### ST0.6 output summary on current internal data

- Recommended strains: `369`.
- Recommendation rows: `1,107` (`3` per strain).
- Diversity relaxation needed: `0` strains under current `max_per_family=2`.
- Holdout top-3 hit rate (all strains): `0.784615`.
- Holdout top-3 hit rate (susceptible strains only): `0.809524`.

#### ST0.6 Interpretation

1. The current simple recommendation layer underperforms the ST0.4 raw top-3 benchmark, indicating that ranking and
   recommendation objectives are not yet aligned end-to-end.
2. Diversity constraint did not bind in this dataset configuration (`0` relaxed strains), so current performance is not
   being driven by diversity tradeoffs.
3. ST0.7 should expose this gap explicitly in final reporting so the next iteration can focus on recommendation-quality
   optimization rather than only calibration quality.

### 2026-02-15: ST0.6b implemented (ranking-policy comparison)

#### What we implemented in ST0.6b

1. Added ST0.6b step: `lyzortx/pipeline/steel_thread_v0/steps/st06b_compare_ranking_policies.py`.
2. Compared six recommendation policies on the same holdout set: raw, platt, and isotonic scores; each with and without
   family-cap diversity.
3. Wrote outputs:
   - `lyzortx/generated_outputs/steel_thread_v0/intermediate/st06b_policy_comparison.csv`
   - `lyzortx/generated_outputs/steel_thread_v0/intermediate/st06b_recommendations_all_policies.csv`
   - `lyzortx/generated_outputs/steel_thread_v0/intermediate/st06b_top3_recommendations_best.csv`
   - `lyzortx/generated_outputs/steel_thread_v0/intermediate/st06b_summary.json`

#### ST0.6b holdout results

- `logreg_platt__none`: top-3 all `0.846154`, susceptible-only `0.873016` (best; tied with `logreg_raw__none`)
- `logreg_raw__none`: top-3 all `0.846154`, susceptible-only `0.873016`
- `logreg_platt__max_family_2`: top-3 all `0.815385`, susceptible-only `0.841270`
- `logreg_raw__max_family_2`: top-3 all `0.815385`, susceptible-only `0.841270`
- `logreg_isotonic__none`: top-3 all `0.800000`, susceptible-only `0.825397`
- `logreg_isotonic__max_family_2` (former ST0.6 policy): top-3 all `0.784615`, susceptible-only `0.809524`

#### ST0.6b Interpretation

1. The ST0.6 drop versus ST0.4 is primarily policy choice, not an implementation bug.
2. In this dataset/split, isotonic ranking is weaker than raw or platt ranking for top-3 hit-rate.
3. Family-cap diversity (`max_per_family=2`) reduces hit-rate for all three score variants in current holdout.
4. Next change should be to switch operational ranking from isotonic to platt (or raw), while keeping calibration
   outputs for probability-quality reporting.
5. Raw and Platt top-3 tie in this setup because Platt is a monotonic remapping of the same raw score; top-k lift should
   therefore be expected from better model signal, better labels, or new data, not monotonic recalibration alone.

### 2026-02-15: ST0.6 policy switch to `logreg_platt__none`

#### What we changed

1. Updated ST0.6 default ranking policy to `score_column=pred_logreg_platt`.
2. Disabled family-cap diversity by default (`max_per_family=0`, diversity mode `none`).
3. Updated ST0.6 and ST0.7 baselines to match the new default policy outputs.

#### Updated ST0.6 holdout metrics

- Top-3 hit rate (all strains): `0.846154` (`55/65`).
- Top-3 hit rate (susceptible-only): `0.873016` (`55/63`).
- Diversity relaxation count: `0` (diversity mode `none`).

#### Updated ST0.7 consequence

- `error_analysis.csv` holdout miss rows decreased from `14` to `10` under the new ST0.6 default policy.

### 2026-02-15: ST0.7 implemented (final report artifacts)

#### What we implemented in ST0.7

1. Added ST0.7 step: `lyzortx/pipeline/steel_thread_v0/steps/st07_build_report.py`.
2. Added ST0.7 regression check: `lyzortx/pipeline/steel_thread_v0/checks/check_st07_regression.py`.
3. Added baseline snapshot: `lyzortx/pipeline/steel_thread_v0/baselines/st07_expected_metrics.json`.
4. Extended CI workflow to run ST0.7 regression in addition to ST0.1 through ST0.6.

#### ST0.7 output summary on current internal data

- Metrics summary rows: `48`.
- Top-3 recommendation rows: `1,107`.
- Calibration summary rows: `12`.
- Error analysis rows (holdout misses): `14`.
- Output artifacts:
  - `lyzortx/generated_outputs/steel_thread_v0/metrics_summary.csv`
  - `lyzortx/generated_outputs/steel_thread_v0/top3_recommendations.csv`
  - `lyzortx/generated_outputs/steel_thread_v0/calibration_summary.csv`
  - `lyzortx/generated_outputs/steel_thread_v0/error_analysis.csv`
  - `lyzortx/generated_outputs/steel_thread_v0/run_manifest.json`

#### ST0.7 Interpretation

1. The steel thread is now complete end-to-end, with deterministic outputs and regression-gated final artifacts.
2. The recommendation gap observed in ST0.6 is now explicit in `error_analysis.csv` for targeted iteration.
3. The immediate next optimization target should be recommendation quality (Top-3 hit rate), not additional plumbing.

### 2026-03-15: ST08 Dual-Slice Reporting Verification (ST0.7)

#### What was implemented

- Added a dedicated regression test for ST0.7 report generation to verify that `metrics_summary.csv` contains separate
  holdout metric rows for both `full_label` and `strict_confidence` slices.
- The test enforces presence of the required acceptance metrics by slice:
  - `topk_hit_rate_all_strains`
  - `brier_score`
  - `ece`

#### Findings

- ST0.7 now has explicit test-level guarantees that dual-slice reporting is present in report outputs and keyed by
  `__full_label` / `__strict_confidence` suffixes for downstream parsing.

#### Interpretation

- This closes a reporting-audit gap between ST0.5/ST0.6 and ST0.7 by ensuring the final exported report preserves slice
  separation for recommendation quality and calibration quality metrics.

### 2026-03-16: ST09 holdout miss failure hypotheses

#### What we analyzed

1. Re-ran ST0.1 through ST0.7 and confirmed `lyzortx/generated_outputs/steel_thread_v0/error_analysis.csv` contains `10`
   current holdout misses under the active `logreg_platt__none` recommendation policy.
2. Ran ST0.6b and confirmed none of the `10` missed strains become hits under the six compared policy variants
   (`raw/platt/isotonic`, with and without family cap), so these are not simple monotonic-calibration or family-cap
   policy issues.
3. Reviewed miss context from `error_analysis.csv`, `st05_pair_predictions_calibrated.csv`,
   `st06_top3_recommendations.csv`, `st06b_recommendations_all_policies.csv`, and `st02_pair_table.csv`.

#### Error buckets

1. **No-susceptibility / abstention failure:** `FN-B4`, `NILS24`. These strains have zero labeled positive holdout
   phages, so the current always-return-top-3 policy is guaranteed to emit false positives.
2. **Single-weak-positive strains buried under broad-host-range priors:** `ECOR-06`, `H1-002-0060-C-T`. Each strain has
   exactly one positive holdout phage, and that positive sits below a dense block of higher-scored `Straboviridae`
   negatives.
3. **Within-family ordering misses among crowded `Straboviridae` candidates:** `ECOR-14`, `NILS41`, `NILS70`, `ROAR205`.
   At least one true positive is already near the top, but top-3 truncation inside a dense same-family score band still
   excludes it.
4. **Cross-family blind spot / `Straboviridae` collapse:** `ECOR-69`, `NILS53`. The recommender top-3 remains entirely
   `Straboviridae`, while most true positives for these strains live in `Other` or `Autographiviridae`.

#### Quantitative notes

- Miss strains are much narrower than hit strains: mean holdout positive count `5.8` vs `27.5` for holdout hits.
- `4/10` miss strains have `legacy_label_breadth_count = 0`; none of the `55` holdout-hit strains do.
- Current ST0.6 recommendations are `Straboviridae` for every holdout strain, not just the missed ones.
- Among true-positive pairs from the `10` missed strains, family counts are `Other = 40`, `Autographiviridae = 9`,
  `Straboviridae = 9`.

#### Strain-level hypotheses and next steps

<!-- pyml disable-num-lines 12 md013-->
| Strain            | Bucket                         | Failure hypothesis                                                                                                                                                                                                                                            | Actionable next step                                                                                                                                                                            |
| ----------------- | ------------------------------ | ------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- | ----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------- |
| `ECOR-06`         | Single-weak-positive           | Only one holdout positive (`AL505_Ev3`) exists, and it is buried under a broad `Straboviridae` block of slightly higher-scored negatives. This looks like narrow-susceptibility recall failure rather than a recommendation-policy bug.                       | Add a narrow-susceptibility audit slice (`n_true_positive_phages_holdout <= 3`) and test whether reweighting these strains or adding a recall-oriented reranker lifts single-positive recovery. |
| `ECOR-14`         | Within-family ordering         | Multiple true positives exist, but the best one lands just outside top-3 behind same-family negatives. The model appears to know the right family but not the right phage within that family.                                                                 | Add a second-stage within-family reranker using per-phage precision / empirical host-range features and check whether `LF73_P1` and related positives move into top-3.                          |
| `ECOR-69`         | Cross-family blind spot        | This strain has many positives, but the top ranks are dominated by `Straboviridae` negatives while most true positives are `Other`. That pattern suggests family-level prior collapse is overpowering host-specific compatibility.                            | Run a targeted ablation that weakens phage-family identity features, then measure whether non-`Straboviridae` positives recover on this strain and on the full miss bucket.                     |
| `FN-B4`           | No-susceptibility / abstention | Holdout has zero positives, so the model should have abstained instead of forcing three phages. This is a product-contract gap, not just a ranking gap.                                                                                                       | Prototype a no-recommend / low-confidence threshold using max score and score-margin features, then evaluate false-abstain vs false-positive tradeoffs on holdout.                              |
| `H1-002-0060-C-T` | Single-weak-positive           | The only positive (`DIJ07_P1`) is a low-support singleton far below a `Straboviridae` block. With `legacy_label_breadth_count = 0` and `host_n_defense_systems = 10`, the current feature set likely misses a rare host-compatibility signal.                          | Prioritize receptor/defense feature work for zero-infection, high-defense strains and test whether `DIJ07_P1` retrieval improves after those host features are added.                           |
| `NILS24`          | No-susceptibility / abstention | Like `FN-B4`, this strain has zero holdout positives, so any fixed top-3 output is necessarily wrong. Its zero prior infections and high defense burden make it a plausible genuinely hard-negative strain.                                                   | Include `NILS24` in the abstention-threshold benchmark set and require any future recommendation layer to permit an explicit no-match outcome.                                                  |
| `NILS41`          | Within-family ordering         | This strain is broadly susceptible, but positives sit inside a large, nearly saturated score band where many `Straboviridae` phages receive similarly high scores. The failure is phage selection inside the right region, not gross family misspecification. | Audit score ties / near-ties around ranks 1-10 and test a tie-aware reranker that uses within-family host-range evidence instead of pure score ordering.                                        |
| `NILS53`          | Cross-family blind spot        | Despite many positives, especially outside `Straboviridae`, the model still pushes a `Straboviridae`-only top-3. Its very high defense burden (`15`) suggests missing host-surface / defense context may be hiding the compatible non-`Straboviridae` phages. | Prioritize host receptor + defense features, then re-evaluate whether `Other` / `Autographiviridae` positives move upward for this strain and the rest of the family-collapse bucket.           |
| `NILS70`          | Within-family ordering         | This is a narrow near-miss: at least one true positive is very close to the cutoff, but top-3 truncation still favors nearby negatives. This looks like ranking-resolution failure around the decision boundary.                                              | Measure score margins for ranks 3-6 on holdout and test whether a boundary-aware reranker or expanded candidate set before reranking improves top-3 hit rate.                                   |
| `ROAR205`         | Within-family ordering         | The top region already contains the correct family and a true positive sits just beyond the cutoff, but sparse prior infection history suggests weak host evidence for choosing the correct member.                                                           | Add a host-neighbor retrieval diagnostic for zero-infection strains and test whether nearest-neighbor host evidence can rerank `BCH953_P4` into the recommended set.                            |

#### ST09 interpretation

1. The remaining `10` misses are primarily model/data failures, not recommendation-policy failures; ST0.6b does not
   rescue any of them.
2. The most important structural issue is family collapse toward `Straboviridae`, which suppresses many true positives
   from `Other` and `Autographiviridae`.
3. The second issue is lack of abstention: two misses are guaranteed errors because the current interface always emits
   three phages even when no labeled susceptible phage exists.
4. The next technically highest-value iteration is not more calibration tuning; it is adding host-compatibility signal
   plus an abstention/no-match mechanism, then re-checking this exact 10-strain bucket.
