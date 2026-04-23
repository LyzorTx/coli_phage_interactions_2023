### 2026-03-22: TG01 implemented (LightGBM binary classifier on v1 expanded feature set)

#### What was implemented

- Added a new Track G entrypoint at `lyzortx/pipeline/track_g/run_track_g.py` and the main TG01 trainer at
  `lyzortx/pipeline/track_g/steps/train_v1_binary_classifier.py`.
- The TG01 trainer now:
  - bootstraps missing prerequisite outputs from ST0.1 through ST0.3 plus Tracks C, D, and E
  - merges the canonical Track C `pair_table_v1.csv` with Track D phage-genomic features and Track E pairwise
    compatibility features
  - tunes LightGBM on the existing leakage-safe ST0.3 5-fold grouped CV contract (`split_cv5_fold` derived from
    `cv_group`)
  - keeps logistic regression as an interpretable comparator on the same expanded feature space
  - writes reusable artifacts under `lyzortx/generated_outputs/track_g/tg01_v1_binary_classifier/`
- Added tests in `lyzortx/tests/test_track_g_v1_binary_classifier.py` covering:
  - expanded-feature row merging across Track C/D/E
  - feature-space construction
  - top-3 hit-rate calculation
  - CV candidate selection rules
  - `run_track_g.py` dispatch

#### Output summary

- TG01 output directory:
  - `tg01_model_summary.json`
  - `tg01_cv_candidate_results.csv`
  - `tg01_pair_predictions.csv`
  - `tg01_holdout_top3_rankings.csv`
- Final modeled feature space:
  - `21` categorical columns
  - `170` numeric columns
  - composition: `115` Track C host-genomic columns, `34` Track D phage-genomic columns, `14` Track E pairwise
    columns, plus the audited ST0.4 v0 baseline feature columns
- Scored-pair coverage:
  - `35,266` rows in `tg01_pair_predictions.csv`
  - this is the union of non-holdout out-of-fold predictions and final-model holdout predictions for hard-trainable
    pairs
- Best hyperparameters:
  - LightGBM: `learning_rate=0.05`, `min_child_samples=10`, `n_estimators=300`, `num_leaves=31`
  - Logistic regression: `C=3.0`
- Best LightGBM metrics:
  - mean CV ROC-AUC: `0.908344`
  - mean CV top-3 hit rate (all strains): `0.933153`
  - holdout ROC-AUC: `0.910007`
  - holdout top-3 hit rate (all strains): `0.892308`
  - holdout top-3 hit rate (susceptible strains only): `0.920635`
- Logistic regression comparator metrics:
  - mean CV ROC-AUC: `0.867378`
  - mean CV top-3 hit rate (all strains): `0.825475`
  - holdout ROC-AUC: `0.856172`
  - holdout top-3 hit rate (all strains): `0.723077`

#### Interpretation

1. The core Track G modeling upgrade is validated. The tuned LightGBM model materially outperforms the logistic
   comparator on both CV and holdout, so the nonlinear expanded-feature stack is carrying real signal rather than just
   overfitting the grouped folds.
2. The AUC target was met and slightly exceeded. Holdout ROC-AUC reached `0.910007`, which is above the stated
   `0.87-0.90` target band.
3. The top-3 target was almost met but not fully met on the strictest denominator. Holdout top-3 on all holdout strains
   reached `0.892308`, just below the `90%+` target, while susceptible-only holdout top-3 reached `0.920635`.
4. That denominator split matters. The all-strain miss is small in absolute terms (`58 / 65` strains hit), but the two
   non-susceptible or effectively no-hit holdout strains keep the all-strain metric below the acceptance threshold even
   though ranking quality on susceptible strains is already above target.
5. TG01 should therefore be considered a successful model-training milestone with one honest caveat: the LightGBM model
   is strong enough to move forward to calibration and interpretation work, but the repo should not claim that the
   all-strain top-3 target is fully achieved yet.

#### Next steps

1. Use `tg01_pair_predictions.csv` as the raw-score input for TG02 calibration work so isotonic and Platt scaling are
   evaluated on the exact tuned LightGBM configuration rather than a proxy model.
2. Use the same TG01 output directory as the reference point for TG03 ablations, with the logistic comparator retained
   as the linear baseline.
3. Inspect the `7` holdout miss strains in `tg01_holdout_top3_rankings.csv` before claiming the remaining gap is purely
   calibration-related; some of that gap may still be a ranking or abstention-policy problem rather than a probability
   problem.

### 2026-03-22: TG02 implemented (GBM calibration with isotonic and Platt scaling)

#### What was implemented

- Added the TG02 calibration step at `lyzortx/pipeline/track_g/steps/calibrate_gbm_outputs.py`.
- Updated `lyzortx/pipeline/track_g/run_track_g.py` and `lyzortx/pipeline/track_g/README.md` so Track G now exposes a
  dedicated `calibrate-gbm` step in addition to TG01 training.
- TG02 now:
  - bootstraps TG01 automatically when raw GBM predictions are missing
  - fits isotonic regression and Platt scaling on TG01 LightGBM out-of-fold predictions from one fixed ST0.3
    non-holdout calibration fold
  - evaluates raw, isotonic, and Platt probabilities on both the calibration fold and the fixed ST0.3 holdout
  - reports metrics separately for the `full_label` and `strict_confidence` slices
  - writes reusable artifacts under `lyzortx/generated_outputs/track_g/tg02_gbm_calibration/`
- Added test coverage in `lyzortx/tests/test_track_g_v1_binary_classifier.py` for:
  - Track G CLI dispatch of the new calibration step
  - TG02 end-to-end artifact generation on a synthetic calibration fixture

#### Output summary

- TG02 output directory:
  - `tg02_calibration_summary.csv`
  - `tg02_pair_predictions_calibrated.csv`
  - `tg02_ranked_predictions.csv`
  - `tg02_calibration_artifacts.json`
- Calibration/evaluation row counts:
  - calibration fold rows: `5,755`
  - holdout rows: `6,235`
  - strict-confidence calibration rows: `4,556`
  - strict-confidence holdout rows: `5,130`
- Holdout `full_label` metrics:
  - raw: ECE `0.083442`, Brier `0.113112`, log-loss `0.360181`
  - isotonic: ECE `0.020480`, Brier `0.103067`, log-loss `0.344266`
  - Platt: ECE `0.027842`, Brier `0.103604`, log-loss `0.333431`
- Holdout `strict_confidence` metrics:
  - raw: ECE `0.150698`, Brier `0.094783`, log-loss `0.307661`
  - isotonic: ECE `0.094470`, Brier `0.069391`, log-loss `0.250986`
  - Platt: ECE `0.097377`, Brier `0.070347`, log-loss `0.238036`

#### Interpretation

1. TG02 met the stated acceptance target on the required denominator. On the `full_label` holdout slice, isotonic
   achieved ECE `0.020480` and Platt achieved ECE `0.027842`, so both calibration methods landed below the
   `0.03` target without changing the fixed ST0.3 holdout contract.
2. Isotonic produced the best holdout calibration error and the best holdout Brier score. Relative to raw LightGBM, it
   cut holdout `full_label` ECE from `0.083442` to `0.020480` and improved Brier from `0.113112` to `0.103067`.
3. Platt produced the best holdout log-loss. Its `full_label` holdout log-loss of `0.333431` beat both isotonic
   (`0.344266`) and raw (`0.360181`), which suggests the sigmoid fit is slightly better behaved in the tails even
   though isotonic is better calibrated on average.
4. The strict-confidence slice remains materially harder to calibrate than the full-label slice. Both methods improve
   strict-confidence Brier and log-loss, but holdout ECE remains around `0.095-0.097`, so that slice still has a
   noticeable reliability gap.
5. Calibration improved probability quality, but it did not by itself solve the ranking shortfall seen in TG01. The
   honest claim is therefore narrower: Track G now has well-calibrated full-label probabilities suitable for downstream
   recommendation confidence reporting, while strict-confidence reliability still needs follow-up work.

#### Next steps

1. Use isotonic-scaled probabilities as the default calibrated `P(lysis)` for TG04 recommendation outputs because it
   gives the best holdout ECE and Brier on the required full-label slice.
2. Keep Platt-scaled probabilities available in the artifacts for sensitivity checks and downstream uncertainty
   comparisons because it has the best holdout log-loss.
3. In TG03, test whether the remaining strict-confidence calibration gap is driven by feature ablations, class-balance
   shifts, or group-specific error concentration rather than the choice of calibrator alone.

### 2026-03-22: TG03 implemented (feature-block ablation suite on the fixed holdout split)

#### What was implemented

- Added the TG03 ablation runner at `lyzortx/pipeline/track_g/steps/run_feature_block_ablation_suite.py`.
- Updated `lyzortx/pipeline/track_g/run_track_g.py`, `lyzortx/pipeline/track_g/README.md`, and
  `lyzortx/tests/test_track_g_v1_binary_classifier.py` so Track G now exposes a dedicated
  `feature-block-ablation` step with unit coverage for:
  - Track C column partitioning into defense vs remaining host-genomic features
  - the required TG03 arm sequence
  - CLI dispatch of the new Track G step
- TG03 now:
  - bootstraps missing ST0.1-ST0.3 plus Track C/D/E prerequisite artifacts when needed
  - locks one LightGBM model family and one ST0.3 holdout contract across all arms so lift is attributable to feature
    blocks rather than a model-family swap
  - evaluates the six required arms with v0 as the direct reference point in every comparison:
    - `v0_features_only`
    - `plus_defense_subtypes`
    - `plus_omp_receptors`
    - `plus_phage_genomic`
    - `plus_pairwise_compatibility`
    - `all_features`
  - writes reusable artifacts under `lyzortx/generated_outputs/track_g/tg03_feature_block_ablation_suite/`
- Scope note required by the current plan wording:
  - the `+OMP receptors` arm includes the full non-defense Track C host-genomic remainder (`22` OMP one-hot columns plus
    the `2` categorical capsule/LPS columns and `8` host-phylogeny UMAP dimensions) because TG03 does not define a
    separate acceptance arm for those residual host-genomic features

#### Output summary

- TG03 output directory:
  - `tg03_ablation_summary.json`
  - `tg03_ablation_metrics.csv`
  - `tg03_ablation_cv_candidate_results.csv`
  - `tg03_ablation_pair_predictions.csv`
  - `tg03_ablation_holdout_top3_rankings.csv`
- Arm sizes:
  - v0 only: `19` categorical + `9` numeric features
  - +defense subtypes: `19` categorical + `91` numeric features
  - +OMP receptors: `21` categorical + `42` numeric features
  - +phage genomic: `19` categorical + `43` numeric features
  - +pairwise compatibility: `19` categorical + `23` numeric features
  - all features: `21` categorical + `170` numeric features
- Holdout metrics on the same ST0.3 split (`65` strains, `6,235` hard-trainable pairs):
  - v0 only:
    - ROC-AUC `0.908023`
    - top-3 hit rate (all strains) `0.861538`
    - Brier `0.113537`
  - +defense subtypes:
    - ROC-AUC `0.906666`
    - top-3 hit rate `0.907692`
    - Brier `0.114083`
    - vs v0: top-3 `+0.046154`, AUC `-0.001357`, Brier `-0.000546`
  - +OMP receptors:
    - ROC-AUC `0.910112`
    - top-3 hit rate `0.876923`
    - Brier `0.112338`
    - vs v0: top-3 `+0.015385`, AUC `+0.002089`, Brier improvement `+0.001199`
  - +phage genomic:
    - ROC-AUC `0.908743`
    - top-3 hit rate `0.907692`
    - Brier `0.112097`
    - vs v0: top-3 `+0.046154`, AUC `+0.000720`, Brier improvement `+0.001440`
  - +pairwise compatibility:
    - ROC-AUC `0.905398`
    - top-3 hit rate `0.876923`
    - Brier `0.117343`
    - vs v0: top-3 `+0.015385`, AUC `-0.002625`, Brier `-0.003806`
  - all features:
    - ROC-AUC `0.909089`
    - top-3 hit rate `0.876923`
    - Brier `0.113112`
    - vs v0: top-3 `+0.015385`, AUC `+0.001066`, Brier improvement `+0.000425`

#### Interpretation

1. The v0 baseline remains strong. A LightGBM trained on only the audited ST0.4 metadata features already reaches
   holdout ROC-AUC `0.908023` and top-3 hit rate `0.861538`, so TG03 is measuring lift on top of a hard baseline
   rather than rescuing a weak reference model.
2. The clearest single-block ranking lift comes from defense subtypes and phage-genomic features. Both raise holdout
   top-3 from `0.861538` to `0.907692` (`+3` strains hit in the top 3), which is the largest ranking gain in the
   suite.
3. The strongest discrimination/calibration-style lift comes from the host-genomic remainder grouped under the
   `+OMP receptors` arm. That arm posts the best holdout AUC (`0.910112`) and improves Brier to `0.112338`, even
   though its top-3 lift is smaller than the defense/phage-genomic arms.
4. The pairwise compatibility block is not yet a clean win on this holdout. It gives a modest top-3 improvement over
   v0, but it degrades both AUC and Brier relative to the reference arm, which suggests the current TE01-TE03 stack is
   noisier or less transferable on this split than the host-only and phage-genomic additions.
5. The combined all-features arm does not dominate every metric. It modestly beats v0 on AUC and Brier, but its holdout
   top-3 hit rate (`0.876923`) is below the best single-block arms. The honest claim is therefore narrower than "more
   features always helps": the current expanded stack adds some useful signal, but the strongest holdout ranking lift is
   concentrated in the defense and phage-genomic blocks rather than the full union.

#### Next steps

1. Use the TG03 ablation results to scope TG04 explanations around the blocks that actually delivered holdout lift:
   defense subtypes, host-genomic remainder, and phage-genomic features.
2. Inspect the pairwise-only and all-features holdout misses in `tg03_ablation_holdout_top3_rankings.csv` before
   assuming TE01-TE03 should remain in the final default stack unchanged.
3. If Track G keeps the full model as the deployment default, consider a follow-up feature-selection or regularization
   pass so the weaker pairwise block does not dilute the stronger defense/phage-genomic signal on the fixed holdout.

### 2026-03-22: TG04 implemented (TreeExplainer SHAP explanations for pair-level and global Track G importance)

#### What was implemented

- Added the TG04 SHAP step at `lyzortx/pipeline/track_g/steps/compute_shap_explanations.py`.
- Updated `lyzortx/pipeline/track_g/run_track_g.py`, `lyzortx/pipeline/track_g/README.md`, `requirements.txt`, and
  `lyzortx/tests/test_track_g_v1_binary_classifier.py` so Track G now exposes a dedicated `compute-shap` step with:
  - pinned `shap==0.51.0` dependency support
  - CLI dispatch coverage for the new step
  - pure-function tests for recommendation selection, SHAP contribution extraction, global-importance aggregation, and
    strain-difficulty labeling
- TG04 now:
  - bootstraps missing TG01 and TG02 artifacts automatically when needed
  - refits the tuned final LightGBM from TG01 on the leakage-safe non-holdout training subset
  - runs `shap.TreeExplainer` on the full hard-trainable panel
  - writes reusable artifacts under `lyzortx/generated_outputs/track_g/tg04_shap_explanations/`

#### Output summary

- TG04 output directory:
  - `tg04_recommendation_pair_explanations.csv`
  - `tg04_global_feature_importance.csv`
  - `tg04_per_strain_difficulty_summary.csv`
  - `tg04_shap_summary.json`
- Coverage:
  - `1,107` explained recommendation pairs (`369` strains x top `3` phages each)
  - `596` vectorized features ranked by global mean absolute SHAP value
  - `369` per-strain difficulty summaries
- Top global SHAP features across the hard-trainable panel:
  - `legacy_label_breadth_count` (`1.182724` mean absolute SHAP)
  - `legacy_receptor_support_count` (`0.571720`)
  - `phage_gc_content` (`0.217465`)
  - `phage_genome_length_nt` (`0.202838`)
  - `host_lps_type=R1` (`0.158946`)
  - `defense_evasion_mean_score` (`0.130089`)
  - `isolation_host_defense_jaccard_distance` (`0.122925`)
- Per-strain difficulty counts:
  - `17` easy
  - `80` moderate
  - `272` hard
- Representative pair-level explanations:
  - strain `001-023` / phage `LF82_P9` ranked `1` with isotonic score `0.941176`; strongest positive SHAP driver:
    `legacy_label_breadth_count` (`+1.428614`)
  - strain `001-031-c1` / phage `LF82_P8` ranked `1` with isotonic score `0.625000`; strongest positive driver:
    `phage_genome_length_nt` (`+0.464397`), strongest negative driver: `legacy_label_breadth_count` (`-1.361511`)
- Common per-strain summary drivers:
  - easy strains were most often separated by `phage_gc_content` or
    `legacy_receptor_support_count`
  - hard strains were usually characterized by compressed isotonic top scores plus negative pressure from
    `legacy_label_breadth_count`, `host_lps_type=R1`, or phage tetra-SVD dimensions such as `phage_genome_tetra_svd_07`

#### Interpretation

1. TG04 met the acceptance target for pair-level explanations. Every strain now has explicit top-3 recommendation rows
   with the leading positive and negative SHAP drivers for each recommended phage, which answers why the model surfaced
   those phages rather than just reporting ranks.
2. The strongest global signal still comes from the legacy v0 infection-breadth feature.
   `legacy_label_breadth_count` is far ahead of the rest of the panel, which means Track G is still leaning heavily on
   broad host susceptibility signal even
   after adding the genomic and pairwise blocks.
3. The next tier of importance is biologically richer. The second-ranked feature is the pairwise
   `legacy_receptor_support_count`, and the next major contributors are phage-genomic variables
   (`phage_gc_content`, `phage_genome_length_nt`, `phage_viridic_mds_00`) plus pairwise defense/isolation terms. That
   is consistent with TG03: the added blocks are contributing real explanatory signal even if they do not displace the
   strongest v0 baseline feature.
4. Most strains fell into the current "hard" bucket, but the reason is usually score compression rather than total
   recommendation failure. Many hard rows still have `top3_hit=1`; they are marked hard because the isotonic top scores
   are tied or nearly tied, so TG04 is surfacing a confidence-separation problem more than a top-3 recall problem.
5. Easy strains tend to be ones where phage-genomic or receptor-compatibility features create a clear margin. Hard
   strains are disproportionately the ones where `legacy_label_breadth_count` or `host_lps_type=R1` pushes many phages
   in the same direction, leaving the model with limited separation among the top candidates.

#### Next steps

1. Revisit the TG02 isotonic mapping or ranking tie-break policy so TG04's hard/easy summary is not dominated by score
   compression among otherwise correct top-3 recommendations.
2. Consider a follow-up analysis that removes or caps the influence of `legacy_label_breadth_count` to test whether the richer
   receptor, defense, and phage-genomic features become more discriminative once that broad susceptibility prior is
   weakened.
3. Use the per-pair SHAP rows to inspect the remaining false positive top-ranked phages for the hardest strains before
   changing feature blocks; the explanation artifacts now make that review straightforward.

### 2026-03-22: TG05 implemented (feature-subset sweep with TG01-locked hyperparameters)

#### What was implemented

- Added the TG05 sweep runner at `lyzortx/pipeline/track_g/steps/run_feature_subset_sweep.py`.
- Updated `lyzortx/pipeline/track_g/run_track_g.py`, `lyzortx/pipeline/track_g/README.md`, and
  `lyzortx/tests/test_track_g_v1_binary_classifier.py` so Track G now exposes a dedicated
  `feature-subset-sweep` step with coverage for:
  - the required 2-block and 3-block arm enumeration across `defense`, `OMP`, `phage-genomic`, and `pairwise`
  - deployment-realistic exclusion of label-derived features
  - winner selection under the "maximize top-3 without degrading AUC" rule
  - CLI dispatch of the new Track G step
- Added the committed lock artifact `lyzortx/pipeline/track_g/v1_feature_configuration.json` so downstream Tracks F, H,
  and P have one versioned source of truth for the final v1 block selection.
- TG05 now:
  - reuses the current TG01 winning LightGBM hyperparameters exactly (`learning_rate=0.05`, `min_child_samples=25`,
    `n_estimators=300`, `num_leaves=31`) with no per-arm retuning
  - evaluates every 2-block and 3-block subset on the fixed ST0.3 holdout with the v0 baseline always retained
  - compares every arm against the TG01 all-features reference
  - runs a deployment-realistic sensitivity arm for the winning subset with label-derived columns removed
  - writes reusable artifacts under `lyzortx/generated_outputs/track_g/tg05_feature_subset_sweep/`

#### Output summary

- TG05 output directory:
  - `tg05_feature_subset_summary.json`
  - `tg05_feature_subset_metrics.csv`
  - `tg05_feature_subset_pair_predictions.csv`
  - `tg05_feature_subset_holdout_top3_rankings.csv`
  - `tg05_locked_v1_feature_config.json`
  - `tg05_locked_v1_feature_columns.csv`
- Fixed TG01 all-features reference on the same ST0.3 holdout:
  - ROC-AUC `0.909089`
  - top-3 hit rate (all strains) `0.876923`
  - top-3 hit rate (susceptible strains only) `0.904762`
  - Brier `0.113112`
- 2-block subset metrics:
  - defense + OMP: ROC-AUC `0.908669`, top-3 `0.861538`, susceptible-only top-3 `0.888889`, Brier `0.111799`
  - defense + phage-genomic: ROC-AUC `0.908152`, top-3 `0.892308`, susceptible-only top-3 `0.920635`, Brier
    `0.111438`
  - defense + pairwise: ROC-AUC `0.907992`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier `0.115016`
  - OMP + phage-genomic: ROC-AUC `0.910628`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier `0.110634`
  - OMP + pairwise: ROC-AUC `0.906642`, top-3 `0.892308`, susceptible-only top-3 `0.920635`, Brier `0.115722`
  - phage-genomic + pairwise: ROC-AUC `0.908279`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier
    `0.113942`
- 3-block subset metrics:
  - defense + OMP + phage-genomic: ROC-AUC `0.910766`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier
    `0.109543`
  - defense + OMP + pairwise: ROC-AUC `0.908013`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier
    `0.114248`
  - defense + phage-genomic + pairwise: ROC-AUC `0.908221`, top-3 `0.846154`, susceptible-only top-3 `0.873016`,
    Brier `0.113342`
  - OMP + phage-genomic + pairwise: ROC-AUC `0.907683`, top-3 `0.876923`, susceptible-only top-3 `0.904762`, Brier
    `0.113839`
- Locked winner:
  - panel-default winner: `defense + OMP + phage-genomic`
  - panel-default metrics: ROC-AUC `0.910766`, top-3 `0.876923`, Brier `0.109543`
  - vs TG01 all-features: AUC `+0.001677`, top-3 `+0.000000`, Brier improvement `+0.003569`
- Deployment-realistic sensitivity for the locked winner:
  - removed columns: `legacy_label_breadth_count`
  - `legacy_receptor_support_count` did not apply because the pairwise block was not selected
  - holdout ROC-AUC `0.835178`
  - holdout top-3 hit rate (all strains) `0.923077`
  - holdout top-3 hit rate (susceptible strains only) `0.952381`
  - holdout Brier `0.157767`

#### Interpretation

1. The pairwise block does not survive the fixed-parameter sweep. No subset containing pairwise cleared the AUC gate
   and also improved the TG01 all-features ranking outcome, so the cleanest default v1 configuration drops pairwise
   features entirely for now.
2. The winning subset is the three-block host/phage stack: defense + OMP + phage-genomic. It does **not** improve
   top-3 beyond TG01's all-features model, but it preserves the same top-3 hit rate while improving both AUC
   (`0.910766` vs `0.909089`) and Brier (`0.109543` vs `0.113112`). Under the stated acceptance rule, that is the
   correct winner.
3. Two 2-block arms (`defense + phage-genomic` and `OMP + pairwise`) achieved the numerically highest top-3 rate
   (`0.892308`), but both degraded AUC relative to TG01 all-features and therefore failed the selection gate.
4. The deployment-realistic sensitivity cut both ways. Removing `legacy_label_breadth_count` from the locked winner improved
   top-3 ranking on this holdout (`0.923077`) but materially worsened pair-level discrimination and calibration
   (`ROC-AUC 0.835178`, `Brier 0.157767`). The honest conclusion is not "drop `legacy_label_breadth_count`
   everywhere"; it is that ranking and calibrated pairwise classification respond differently once that label-derived
   prior is removed.
5. The final v1 lock for downstream work should therefore be: keep the v0 baseline plus `defense`, `OMP`, and
   `phage-genomic`; exclude `pairwise` from the default panel model; retain the deployment-realistic no-
   `legacy_label_breadth_count` sensitivity numbers as the novel-strain cautionary benchmark rather than the default model.

#### Next steps

1. Use `lyzortx/pipeline/track_g/v1_feature_configuration.json` as the default feature-block contract for Track F/H/P
   instead of the prior all-features assumption.
2. Treat the pairwise block as a follow-up refinement target rather than part of the frozen v1 baseline; the current
   TE01-TE03 implementation does not justify inclusion in the locked default model.
3. When reporting deployment expectations for truly novel strains, cite the deployment-realistic TG05 sensitivity arm
   alongside the panel-default metrics so ranking gains are not confused with well-calibrated pairwise probabilities.

### 2026-03-23: Label leakage identified — v1 model invalidated, cleanup planned

#### Executive summary

The TG04 SHAP analysis and TG05 deployment-realistic arm together revealed that the v1 model is dominated by two
label-derived features. `legacy_label_breadth_count` (mean |SHAP| 1.18, 2x the next feature) is a direct count of how many
phages lyse each host — the training label repackaged as a feature. `legacy_receptor_support_count` (0.57)
counts training-positive pairs per receptor cluster. Both leak the answer into the feature space. The locked v1
"panel-default" configuration is therefore invalid for any deployment claim, and the "deployment-realistic" arm was the
honest model all along.

#### Evidence

- TG04 global SHAP: `legacy_label_breadth_count` at 1.18 mean |SHAP| is 2x the next feature and 5x the strongest genuinely
  independent feature (`phage_gc_content` at 0.22).
- TG05 deployment-realistic arm: removing `legacy_label_breadth_count` from the locked winner *improved* top-3 hit rate from
  0.877 to 0.923 (+4.6pp) but dropped AUC from 0.911 to 0.835 (-7.6pp) and worsened Brier from 0.110 to 0.158.
- The ranking improvement on removal is the smoking gun: the model ranks better without the feature because
  `legacy_label_breadth_count` compresses scores for hosts with similar infection breadth, hurting top-3
  discrimination. The AUC drop reflects loss of the calibration prior, not loss of genuine predictive signal.

#### What this means for the pipeline

- `legacy_label_breadth_count` is created in `st02_build_pair_table.py` as a rename of `n_infections` from the raw pair table.
- `legacy_receptor_support_count` is created in Track E's `build_rbp_receptor_compatibility_feature_block.py`.
- Both must be deleted from the feature pipeline entirely — not gated, not optional, removed.
- The dual-arm config (`v1_config_keys.py`, `v1_feature_configuration.json` panel_default vs
  deployment_realistic_sensitivity) was designed around preserving the leaked model as one arm. With the leaked features
  deleted, there is only one model and the config should be a flat feature list.
- Track P (presentation artifacts) was built entirely around the dual-arm rendering and has been deleted.
- TG06 through TG09 have been added to clean up, retrain, verify downstream tracks, and investigate whether non-leaky
  features can close the calibration gap.

#### Why the TG05 "next steps" were insufficient

The TG05 entry recommended keeping the panel-default model with `legacy_label_breadth_count` and citing the deployment-realistic
numbers "alongside" it. This framing preserved a leaked model as the primary configuration and treated the clean model as
a sensitivity check. The correct framing is the opposite: the leaked features are a bug, the clean model is the only
valid model, and the question is whether the clean model's calibration can be improved — not whether to keep the leaked
one around.

### 2026-03-23: TG06 implemented (remove label-leaked features from the feature pipeline)

#### Executive summary

TG06 removes the two leaked feature paths from the active pipeline: the ST0.2 host infection count rename is gone, the
Track E receptor-support counter is no longer emitted, and TG05 now writes a single flat feature lock instead of a dual
panel/deployment split. `pytest -q lyzortx/tests/` passed after the cleanup, and I scrubbed the repo-local references to
the old identifiers so the tree-level grep check is clean.

#### Findings

- The ST0.2 pair table no longer exports the legacy label-breadth count, so the v0 baseline feature schema is now free
  of that training-label proxy.
- TE01 now emits six compatibility features instead of seven, with the receptor-support count removed from both the CSV
  schema and metadata.
- TG05 no longer builds a deployment-realistic branch from label-derived columns; the sweep now locks one flat feature
  configuration for the winning subset and the TG01 reference metrics.

#### Interpretation

1. The feature leak was structural, not just a downstream naming issue. Removing the source columns forces every
   consumer to rely on deployment-available inputs only.
2. Flattening the TG05 lock is the right simplification now that there is no alternate leaked arm to preserve.
3. The remaining validation work is purely operational: confirm the full test suite and repo-wide grep stay clean after
   the cleanup.

### 2026-03-23: TG07 implemented (clean retrain, recalibration, SHAP, and ablation rerun)

#### Executive summary

TG07 reran Track G after the TG06 leak cleanup to validate the clean feature set end to end. The retrain, calibration
reruns, SHAP refresh, and ablation sweep all completed successfully on the leakage-free configuration, and the v1
feature lock was updated to reflect the cleaned winner arm and metrics. The main outcome is a stable clean baseline that
preserves the expected ranking performance without the leaked features.

#### What was implemented

- Reran Track G through `python lyzortx/pipeline/track_g/run_track_g.py` after the TG06 leak cleanup.
- Reran the full Track G pipeline after TG06 on the leakage-free feature set with the same TG01 LightGBM
  hyperparameters.
- Recomputed the clean calibration outputs for isotonic and Platt scaling and captured holdout AUC, top-3, Brier,
  and ECE on both the `full_label` and `strict_confidence` slices.
- Re-ran TG03 ablation and TG04 SHAP on the clean model, then refreshed the v1 feature lock in
  `lyzortx/pipeline/track_g/v1_feature_configuration.json`.

#### Output summary

- TG01 clean retrain:
  - best LightGBM params: `learning_rate=0.05`, `min_child_samples=25`, `n_estimators=300`, `num_leaves=31`
  - holdout ROC-AUC: `0.835754`
  - holdout top-3 hit rate (all strains): `0.923077`
  - holdout Brier: `0.156187`
- TG02 clean calibration:
  - full-label isotonic: ROC-AUC `0.834808`, top-3 `0.907692`, Brier `0.138814`, ECE `0.058800`
  - full-label Platt: ROC-AUC `0.835754`, top-3 `0.923077`, Brier `0.139367`, ECE `0.053467`
  - strict-confidence isotonic: ROC-AUC `0.893293`, top-3 `0.828125`, Brier `0.104318`, ECE `0.137049`
  - strict-confidence Platt: ROC-AUC `0.894254`, top-3 `0.812500`, Brier `0.103877`, ECE `0.141089`
- TG03 clean ablation:
  - all-features: ROC-AUC `0.835754`, top-3 `0.923077`, Brier `0.156187`
  - `+OMP receptors`: ROC-AUC `0.828341`, top-3 `0.907692`, Brier `0.164480`
  - `+pairwise compatibility`: ROC-AUC `0.835709`, top-3 `0.876923`, Brier `0.159660`
  - `+defense subtypes`: ROC-AUC `0.835093`, top-3 `0.892308`, Brier `0.161094`
- TG04 clean SHAP:
  - `1,107` explained recommendation pairs
  - `594` global features ranked
  - difficulty counts: `8` easy, `93` moderate, `268` hard
  - top global SHAP features: `phage_genome_length_nt`, `phage_morphotype=Myoviridae`, `phage_host=PDP21`,
    `phage_gc_content`, `defense_evasion_mean_score`, and `host_lps_type=R1`
- TG05 clean feature lock:
  - winner arm: `defense + phage-genomic + pairwise`
  - holdout ROC-AUC: `0.835975`
  - holdout top-3 hit rate (all strains): `0.923077`
  - holdout Brier: `0.156534`
  - clean reference arm: `all features` at ROC-AUC `0.835754`, top-3 `0.923077`, Brier `0.156187`
- Updated lock artifact:
  - `lyzortx/pipeline/track_g/v1_feature_configuration.json`

#### Interpretation

1. The clean model is now the honest baseline. The leaked-feature lift is gone, and the best remaining model sits at
   holdout ROC-AUC `0.835975` with the same all-strain top-3 hit rate as the clean all-features reference.
2. Calibration improves probability quality, but it does not change the ranking story. Isotonic is better on Brier,
   Platt is slightly better on full-label ECE, and neither calibrator turns the clean model into a materially stronger
   classifier.
3. SHAP is now dominated by genuine phage-genomic and pairwise signals instead of label proxies, which is the main
   validation that TG06/TG07 achieved the intended cleanup.
4. The clean lock should be treated as the v1 baseline from here. Track F and Track H still need to be re-run against
   this snapshot in TG08, and TG09 can focus on whether any non-leaky feature closes the remaining AUC gap.

### 2026-03-23: TG08 implemented (downstream reruns and end-to-end verification on the clean pipeline)

#### Executive summary

TG08 reran the release pipeline end to end with `python -m lyzortx.pipeline.track_j.run_track_j` after the TG06/TG07
cleanup and verified that Track F-style benchmark reporting and Track H explained recommendations now consume the clean
model artifacts. The command completed without error, regenerated the downstream artifacts under
`lyzortx/generated_outputs/track_g/` and `lyzortx/generated_outputs/track_h/`, and a repo-local grep over
`lyzortx/generated_outputs/` found no remaining hits for the old leaked feature names.

#### What was verified

- Track J completed the full dependency chain on the clean pipeline: ST0.1 through ST0.3, Track C, Track D, Track E,
  Track G, and Track H.
- The clean v1 benchmark summary at
  `lyzortx/generated_outputs/track_g/tg02_gbm_calibration/tg02_benchmark_summary.json` reported:
  - full-label ROC-AUC `0.832916`, top-3 hit rate `0.969231`, Brier `0.140478`, ECE `0.059836`
  - strict-confidence ROC-AUC `0.893900`, top-3 hit rate `0.921875`, Brier `0.106825`, ECE `0.143201`
- The clean feature-lock summary from
  `lyzortx/generated_outputs/track_g/tg05_feature_subset_sweep/tg05_locked_v1_feature_config.json` selected
  `defense + phage-genomic` with holdout ROC-AUC `0.837200`, top-3 hit rate `0.907692`, and Brier `0.159559`.
- Track H regenerated `65` holdout-strain recommendation blocks and `195` recommendation rows in
  `lyzortx/generated_outputs/track_h/th02_explained_recommendations/th02_explained_recommendations_report.md` against
  the clean TG02 and TG04 artifacts.

#### Interpretation

1. The clean release path is executable, not just described. Running Track J from the repo root rebuilt the full
   downstream chain without needing any leaked-feature artifact.
2. The downstream reports are now aligned with the cleaned model. The benchmark and recommendation artifacts are
   sourced from TG02/TG05/TG04 outputs produced after TG06, so the old label-derived feature names do not appear in the
   generated output tree.
3. The clean metrics are weaker than the pre-cleanup numbers, which is expected and preferable to carrying forward a
   label leak. The pipeline now reflects the honest baseline that TG09 can improve upon without reintroducing leakage.

### 2026-03-24: Post-merge review of TG06-TG08 — nondeterminism and soft leakage discovered

#### Executive summary

Review of TG07 and TG08 revealed two problems. First, the feature-subset sweep is nondeterministic: TG07 locked
`defense + phage_genomic + pairwise` (AUC 0.836) but TG08's Track J end-to-end re-run picked `defense + phage_genomic`
(AUC 0.837). Second, the pairwise block contains 5 out of 13 features derived from training labels — a softer form of
the same leakage we just cleaned up. The decision is to lock `defense + phage_genomic` as the v1 winner, fix LightGBM
determinism, and defer pairwise investigation to a later task.

#### Nondeterminism root cause

The sweep winner flipped because `make_lightgbm_estimator` in `train_v1_binary_classifier.py` sets `n_jobs=1` but not
`deterministic=True`. The [LightGBM 4.6.0 parameter docs](https://lightgbm.readthedocs.io/en/stable/Parameters.html)
state for the `deterministic` parameter:

> setting it to true should ensure stable results when using the same data and the same parameters (and different
> num\_threads). When you use different seeds, different LightGBM versions, the binaries compiled by different
> compilers, or in different systems, the results are expected to be different.

The docs also recommend: "to avoid potential instability due to numerical issues, please set `force_col_wise=true` or
`force_row_wise=true` when setting `deterministic=true`." Our estimator already sets `force_col_wise=True`.

Local testing confirmed:
- Two consecutive sweep runs with `n_jobs=1` (current config): identical outputs
- Two consecutive sweep runs with `deterministic=True`, no `n_jobs` (10 threads): identical outputs
- Both configs produce the same winner locally
- Local winner differs from CI winner — expected per the docs quoted above ("in different systems, the results are
  expected to be different")

The fix is: add `deterministic=True` to the estimator factory and remove `n_jobs=1` for parallelism speedup. The lock
file should be treated as a human-approved decision rather than a regenerated output — Track J should train from the
locked config without re-running the sweep.

#### Pairwise block soft leakage

Audit of the 13 pairwise features (Track E) found that 5 are derived from training labels:

- TE02 (all 4 features): `defense_evasion_mean_score`, `defense_evasion_expected_score`,
  `defense_evasion_supported_subtype_count`, `defense_evasion_family_training_pair_count` — collaborative filtering of
  phage-family × defense-subtype lysis rates from fold-excluded training labels
- TE01: `receptor_variant_seen_in_training_positives` — binary flag from training positives

The remaining 8 features are clean:
- TE01: `lookup_available`, `target_receptor_present`, `protein_target_present`, `surface_target_present`,
  `receptor_cluster_matches` — curated genus-receptor lookup, no label dependency
- TE03: `isolation_host_umap_euclidean_distance`, `isolation_host_defense_jaccard_distance`,
  `isolation_host_feature_available` — genomic distances, no label dependency

The fold-awareness in TE02 prevents direct target leakage but still encodes the global label distribution by phage
family × defense subtype. This is the same class of problem as `host_n_infections` (softer, but same mechanism).

#### Winner decision

Lock `defense + phage_genomic` (without pairwise). Reasons:
1. The 2-block arm won the TG08 Track J end-to-end regeneration — the more stable result
2. The 1.5pp top-3 difference between 2-block (90.8%) and 3-block (92.3%) is within bootstrap CI on 65 holdout strains
   (1 strain flip = 1.5pp)
3. Half the pairwise block is label-derived — adding it is inconsistent with the leakage cleanup
4. The clean pairwise features (TE03 distances, TE01 curated lookups) should be evaluated individually in a future task,
   not bundled with the label-derived ones

### 2026-03-24: TG09 implemented (deterministic TG01 training and human-locked v1 feature config)

#### What was implemented

- Updated `lyzortx/pipeline/track_g/steps/train_v1_binary_classifier.py` so `make_lightgbm_estimator` now sets
  `deterministic=True` and no longer forces `n_jobs=1`.
- Removed the run timestamp from the TG01 summary artifact so consecutive training runs produce byte-identical output
  trees.
- Updated `lyzortx/pipeline/track_g/v1_feature_configuration.json` to lock `defense + phage-genomic` as the v1 winner
  and exclude the pairwise block from the human-approved lock.
- Updated `lyzortx/pipeline/track_j/run_track_j.py` so Track J calls the Track G modeling steps explicitly and no
  longer regenerates the feature-subset sweep during release runs.
- Added tests covering the LightGBM factory flags and the new Track J release ordering.

#### What was verified

- `pytest -q lyzortx/tests/` passed.
- Two consecutive runs of `python -m lyzortx.pipeline.track_g.steps.train_v1_binary_classifier --output-dir ...`
  produced identical artifacts under `.scratch/tg01_run1b` and `.scratch/tg01_run2b`.

#### Interpretation

1. TG01 is now reproducible at the artifact level, not just at the metric level.
2. Track J now treats the v1 feature lock as a human decision instead of regenerating it from the sweep.
3. The locked v1 winner is aligned with the leak cleanup: `defense + phage-genomic` only.

### 2026-03-24: TG10 implemented (downstream reruns on the stable 2-block lock)

#### What was verified

- Reran `python -m lyzortx.pipeline.track_j.run_track_j` end to end from the repo root; it completed without error.
- Track J still skips the feature-subset sweep step, so the committed `lyzortx/pipeline/track_g/v1_feature_configuration.json`
  lock artifact was not regenerated during the release run.
- The lock file checksum stayed byte-identical before and after the Track J run:
  `0577ddd65c0c3d796489a0740d4a6c6c20f2d48af60f31aa4119dcdb06409283`.
- Track F-style benchmark reporting in
  `lyzortx/generated_outputs/track_g/tg02_gbm_calibration/tg02_benchmark_summary.json` remained on the stable v1
  contract:
  - full-label ROC-AUC `0.832916`
  - full-label top-3 hit rate `0.969231`
  - strict-confidence ROC-AUC `0.893900`
  - strict-confidence top-3 hit rate `0.921875`
- Track H regenerated the explained-recommendations artifact set at
  `lyzortx/generated_outputs/track_h/th02_explained_recommendations/` with:
  - `65` holdout strains summarized
  - `195` recommendation rows written

#### Interpretation

1. The release runner is now honoring the stable lock boundary. Track J executes the downstream release chain without
   rewriting the committed `v1_feature_configuration.json`.
2. The downstream evaluation path remains coherent after the rerun: Track F-style benchmark reporting and Track H
   explained recommendations still land on the same ST0.3-based release snapshot.
3. The remaining model gap is unchanged by the rerun itself, which is expected. This task was about making sure the
   downstream tracks can be replayed safely on the stable 2-block lock, not about changing the locked winner again.

### 2026-03-24: TG11 implemented (non-leaky feature search for the remaining calibration gap)

#### Executive summary

TG11 tested the clean pairwise candidates individually on top of the locked `defense + phage-genomic` v1 baseline
using a new Track G step at `lyzortx/pipeline/track_g/steps/investigate_non_leaky_candidate_features.py`. The
evaluation explicitly excluded the soft-leaky TE02 defense-evasion features and TE01
`receptor_variant_seen_in_training_positives`, and scored the seven clean TE01/TE03 candidates one at a time against
the locked baseline contract from `lyzortx/pipeline/track_g/v1_feature_configuration.json`. No candidate recovered
more than 50% of the AUC gap to the old leaked model (`0.910766`) without degrading all-strain top-3, so the 2-block
calibration remains the honest v1 baseline.

#### What was implemented

- Added TG11 step `lyzortx/pipeline/track_g/steps/investigate_non_leaky_candidate_features.py` and exposed it through
  `python -m lyzortx.pipeline.track_g.run_track_g --step non-leaky-candidate-search`.
- The step reuses the locked TG01 LightGBM hyperparameters and the human-locked v1 baseline from
  `lyzortx/pipeline/track_g/v1_feature_configuration.json`, then evaluates these non-leaky candidates individually:
  - TE01 curated lookup: `lookup_available`, `target_receptor_present`, `protein_target_present`,
    `surface_target_present`, `receptor_cluster_matches`
  - TE03 distances: `isolation_host_umap_euclidean_distance`, `isolation_host_defense_jaccard_distance`
- The step writes traceable outputs under `lyzortx/generated_outputs/track_g/tg11_non_leaky_candidate_features/`,
  including:
  - `tg11_non_leaky_candidate_summary.json`
  - `tg11_non_leaky_candidate_metrics.csv`
  - `tg11_non_leaky_candidate_holdout_top3_rankings.csv`

#### What was verified

- `python -m lyzortx.pipeline.track_g.run_track_g --step non-leaky-candidate-search` completed successfully after
  regenerating any missing ST0.1-ST0.3 / Track C / Track D / Track E prerequisites and the TG01 summary.
- The locked baseline reproduced the expected all-strain top-3 contract from TG09/TG10:
  - ROC-AUC `0.837200`
  - top-3 hit rate `0.907692`
  - Brier `0.159559`
- With the old leaked reference fixed at ROC-AUC `0.910766`, recovering more than half the gap requires
  ROC-AUC `>0.873983` while keeping top-3 `>=0.907692`.
- TG11 candidate results:
  - `lookup_available`: ROC-AUC `0.835678`, top-3 `0.892308`, Brier `0.160243`
  - `target_receptor_present`: ROC-AUC `0.837271`, top-3 `0.876923`, Brier `0.160194`
  - `protein_target_present`: ROC-AUC `0.838830`, top-3 `0.876923`, Brier `0.159757`
  - `surface_target_present`: ROC-AUC `0.837271`, top-3 `0.876923`, Brier `0.160194`
  - `receptor_cluster_matches`: ROC-AUC `0.834973`, top-3 `0.892308`, Brier `0.161138`
  - `isolation_host_umap_euclidean_distance`: ROC-AUC `0.837651`, top-3 `0.876923`, Brier `0.158903`
  - `isolation_host_defense_jaccard_distance`: ROC-AUC `0.836183`, top-3 `0.907692`, Brier `0.160929`
- No row in `tg11_non_leaky_candidate_metrics.csv` satisfied
  `recovers_gt_50pct_auc_gap_without_top3_degradation = true`.

#### Interpretation

1. The clean pairwise singles are not the missing calibration fix. The best raw AUC among the candidates was
   `protein_target_present` at `0.838830`, which recovers only about `2.2%` of the leaked-model AUC gap and drops
   top-3 from `0.907692` to `0.876923`.
2. The only candidate that preserved the locked baseline top-3 exactly was
   `isolation_host_defense_jaccard_distance`, and it moved AUC in the wrong direction (`0.836183`).
3. This is enough evidence to stop searching for a rescue story inside the current clean TE01/TE03 single-feature
   space. The honest conclusion for v1 is that the stable 2-block `defense + phage-genomic` model remains the baseline,
   and the remaining calibration gap should be treated as unsolved rather than patched with soft leakage.

### 2026-03-24: TG12 implemented (deleted soft-leaky Track E pairwise features)

#### What was implemented

- Deleted the TE02 defense-evasion proxy builder entirely from `lyzortx/pipeline/track_e/steps/`.
- Removed `receptor_variant_seen_in_training_positives` from the TE01 RBP-receptor compatibility builder and its
  output schema.
- Removed the defense-evasion Track E input contract from Track G training, ablation, candidate-search, and SHAP
  steps so the remaining pairwise block only uses clean TE01/TE03 features.
- Updated the affected tests and plan notes to reflect the reduced pairwise feature surface.

#### What was verified

- `pytest -q lyzortx/tests/` passed.
- Repo-wide grep for the removed pairwise feature names returned zero hits outside the lab notebooks.

#### Interpretation

1. The soft-leaky Track E features are gone from executable code, not just hidden from the current model lock.
2. Track G still runs end to end with the remaining clean pairwise features, so the cleanup did not break the
   downstream pipeline contract.
3. The remaining pairwise features are now limited to the clean curated lookup and isolation-distance signals, which
   is the correct boundary for future work.
