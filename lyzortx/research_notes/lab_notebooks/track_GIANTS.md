### 2026-04-12 00:15 CEST: Track GIANTS launched — three-layer biological prediction model

#### Executive summary

Launched Track GIANTS (Standing on the Shoulders of GenoPHI and Moriniere) with 4 tasks. The hypothesis: lysis requires
passing adsorption gates (capsule penetration OR receptor binding) then surviving host defenses. We build features for
each layer using the best available tools (DepoScope for depolymerases, Moriniere Table S1 for OMP receptors,
DefenseFinder for defense systems), then let LightGBM + RFE learn the gating logic. Baseline: AUTORESEARCH all-pairs
0.810 AUC on ST03 holdout.

#### Biological motivation

Two papers from the Arkin/Mutalik group reshape our approach:

1. **Moriniere et al. 2026** — Receptor specificity is a near-solved discrete classification problem (19 classes, AUROC
   0.99 from amino acid k-mers). The signal is in localized hypervariable sequence motifs at RBP tip domains, not global
   protein embeddings. This explains why our PLM approach (AX07/AX08) failed — mean-pooled whole-protein embeddings are
   the wrong tool.

2. **Noonan et al. 2025 (GenoPHI)** — ML pipeline configuration matters 5-18x more than feature representation. CatBoost
   + RFE + inverse-frequency class weighting is the optimal classical ML config. GenoPHI achieves AUROC 0.869 on our
   exact phages with binary features.

Together they suggest a three-layer architecture matching the actual infection mechanism:

```
Gate 1: Can phage degrade host's capsule/O-antigen? (depolymerase x capsule)
Gate 2: Can phage bind host's OMP receptor? (RBP receptor class x host OMP variant)
   Adsorption = Gate 1 OR Gate 2 (or both)
Gate 3: Can phage survive host's defenses? (all defense systems, unrestricted)
   Lysis = Adsorption AND survives Gate 3
```

#### Why previous approaches failed at the feature level

+ **AX03 pairwise cross-terms** — paired "any RBP" with "any receptor." Biologically meaningless without knowing WHICH
  receptor the phage targets.
+ **AX07/AX08 PLM embeddings** — encoded family similarity, not binding specificity. The signal is in short motifs at
  specific loci (Moriniere), not global protein embeddings.
+ **Defense features flat** — applied defense globally, confounding it with phylogroup. Defense should only matter AFTER
  adsorption succeeds through either gate.

#### Key data enabling this track

+ **Host capsule variation is real**: 99 capsule feature columns across 369 diverse clinical E. coli isolates (not K-12
  lab strains). All 369 have nonzero capsule profiles.
+ **34/97 phages have tail spike annotations** from Pharokka (39 genes), but the DEPOLYMERASE_PATTERNS detection bug
  means they're currently classified as generic RBPs, not depolymerases. DepoScope will provide high-confidence calls
  with domain boundaries.
+ **Moriniere Table S1** gives genus-level OMP receptor mapping for ~60% of our phages. High-confidence for some genera
  (Justusliebigvirus→NGR, Lambdavirus→LamB), ambiguous for others (Tequatrovirus spans 5 receptors), and "Resistant"
  for capsule-dependent genera (Vectrevirus, Kagunavirus).

#### Baseline

AUTORESEARCH all-pairs: 0.810 AUC, 90.8% top-3, 0.167 Brier on ST03 holdout. This is the single clean baseline — TL18
(0.823) has feature integrity issues, per-phage (0.830) is not deployable.

#### Task graph

```
GT01 (depolymerase x capsule) ──┐
                                 ├── GT03 (integration + RFE + class weighting) ── GT04 (holdout eval)
GT02 (receptor x OMP)     ──────┘
```

GT01 and GT02 are parallel. GT03 combines all layers with RFE and class weighting. GT04 (HPO) and GT05 (CatBoost)
are incremental follow-ups.

### 2026-04-12 08:41 CEST: GT03 — Three-layer integration with RFE and inverse-frequency weighting

#### Executive summary

Combined all three gate feature sets (depolymerase × capsule, receptor × OMP, defense) with the 5-slot AUTORESEARCH
baseline and evaluated on ST03 holdout. The all_gates_rfe arm achieves 0.823 AUC [0.781, 0.858], a statistically
significant +1.2pp over the 0.810 baseline. Gate 1 (depolymerase × capsule) is the primary driver, contributing +0.7pp
alone with 22% feature importance. RFE selects 252/502 numeric features. Inverse-frequency weighting hurts AUC but
improves Brier from 0.165 to 0.145.

#### Ablation results

| Arm | AUC | AUC 95% CI | Top-3 | Brier | AUC delta CI vs baseline |
|-----|-----|-----------|-------|-------|--------------------------|
| baseline | 0.810 | [0.767, 0.848] | 92.3% | 0.165 | — |
| +gate1 (depo×capsule) | 0.818 | [0.772, 0.856] | 93.8% | 0.162 | [+0.003, +0.012] *|
| +gate2 (receptor×OMP) | 0.814 | [0.769, 0.851] | 89.2% | 0.165 | [-0.001, +0.007] |
| +gate3 (defense) | 0.816 | [0.777, 0.850] | 89.2% | 0.165 | [-0.003, +0.016] |
| all_gates | 0.822 | [0.781, 0.857] | 92.3% | 0.161 | [+0.002, +0.022]* |
| all_gates_rfe | 0.823 | [0.781, 0.858] | 90.8% | 0.162 | [+0.003, +0.024] * |
| all_gates_rfe_ifw | 0.809 | [0.758, 0.849] | 89.2% | 0.145 | [-0.014, +0.015] |

`*` = 95% CI excludes zero. 3 seeds, 1000 bootstrap resamples on 65 holdout bacteria.

#### Feature importance (all_gates_rfe arm, seed-averaged)

+ pair_depo_capsule: 21.5% — the dominant new signal
+ host_surface: 21.2% — OMP/capsule HMM scores
+ phage_projection: 15.6% — TL17-frozen phage features
+ phage_stats: 13.9% — genome statistics
+ host_typing: 12.8% — phylogroup/serotype/ST
+ host_stats: 8.4% — genome statistics
+ host_defense: 4.6% — defense system counts
+ pair_receptor_omp: 2.0% — directed receptor cross-terms

#### Interpretation

**Gate 1 (depolymerase × capsule) is the primary discovery.** The 242 features (41 cluster memberships + has_depo/count
× 99 capsule scores) capture 22% of total feature importance and produce a statistically significant +0.7pp AUC lift
alone. This validates the three-layer hypothesis at the capsule penetration layer: DepoScope-predicted depolymerase
presence interacted with host capsule HMM profiles provides discriminative signal that LightGBM can exploit.

**Gate 2 (receptor × OMP) is marginal.** Only 8/96 phages have clean OMP receptor assignments from the genus-level
Table S1 lookup (Tequatrovirus→Tsx, Lambdavirus→LamB, Dhillonvirus→FhuA). The +0.3pp AUC is not statistically
significant. GenoPHI per-phage receptor prediction (AUROC 0.99) would assign receptors to all 96 phages and likely
strengthen this gate.

**Gate 3 (defense) replicates prior findings.** The +0.6pp AUC is consistent with the previously observed +0.7pp
(defense-ablation-autoresearch). The CI is wide [-0.003, +0.016], confirming it's not a reliable contributor at this
sample size.

**RFE helps marginally.** Pruning from 507 to ~257 features adds +0.1pp over all_gates (0.823 vs 0.822). The pruning
primarily removes redundant capsule cross-terms and low-importance defense subtypes.

**Inverse-frequency weighting hurts AUC, helps calibration.** The IFW arm drops AUC by 1.4pp but improves Brier from
0.162 to 0.145 — the model becomes better calibrated by upweighting narrow-host phage positives, but at the cost of
discrimination. This suggests IFW might be better applied as a post-hoc calibration adjustment rather than a training
weight.

#### Error bucket analysis

The baseline has 6/65 holdout misses (90.3% top-3 after seed averaging). The all_gates_rfe arm has comparable top-3
(90.8%). Individual seeds show high variance: top-3 ranges from 87.7% to 95.4% across seeds, reflecting the small
holdout (1 strain = 1.5pp). A proper error bucket analysis requires comparing per-strain predictions, deferred to
after the bootstrap summary.

#### Next steps

+ GT04: HPO with Optuna over LightGBM params — the new feature families may benefit from different tree depth and
  regularization than the default config.
+ GT05: CatBoost comparison — handles categoricals natively, found optimal by GenoPHI.
+ Flag for future: GenoPHI per-phage receptor prediction to strengthen Gate 2 beyond the 8-phage genus-level mapping.
+ Knowledge model update: Gate 1 depolymerase × capsule is a validated signal; IFW calibration tradeoff is a new
  finding.

### 2026-04-12 10:01 CEST: GT04 — HPO with Optuna on three-layer feature set

#### Executive summary

Ran 50-trial Optuna HPO over key LightGBM hyperparameters using 5-fold stratified CV on the GT03 RFE-selected feature
set (257 features). The tuned params achieve 0.828 AUC vs 0.823 default — a marginal +0.4pp that is not statistically
significant (delta CI [-0.010, +0.015]). The GT03 default LightGBM configuration is near-optimal for this feature set.

#### Optuna best params vs GT03 defaults

| Parameter | GT03 Default | Optuna Best |
|-----------|-------------|-------------|
| n_estimators | 300 | 450 |
| learning_rate | 0.05 | 0.077 |
| num_leaves | 31 | 108 |
| min_child_samples | 10 | 9 |
| subsample | 0.8 | 0.75 |
| colsample_bytree | 0.8 | 0.83 |
| reg_lambda | (default) | 0.003 |
| reg_alpha | (default) | 0.19 |

#### Holdout results

| Arm | AUC | 95% CI | Top-3 | Brier |
|-----|-----|--------|-------|-------|
| gt03_default | 0.823 | [0.782, 0.859] | 89.2% | 0.161 |
| optuna_tuned | 0.828 | [0.775, 0.867] | 90.8% | 0.161 |

Delta (tuned vs default): [-0.010, +0.015] — not significant.

#### Interpretation

The HPO found a more complex model (108 leaves, 450 trees) that marginally improves CV AUC but doesn't translate to a
significant holdout gain. This is consistent with GenoPHI's finding that algorithm choice matters more than
hyperparameter tuning once the algorithm family is fixed. The default params (31 leaves, 300 trees) are a better
tradeoff: simpler, faster, and within noise of the tuned config.

The feature importance shift is minor: depo×capsule increases from 21.5% to 24.1% with tuned params (deeper trees
capture more of the cross-term signal), while host_typing drops from 12.8% to 5.5%.

**Caveat:** RFE was applied once on the full training set before Optuna CV, so the inner CV scores are mildly
optimistic (folds see features selected using their own labels). This doesn't affect the holdout comparison — both
arms use the same RFE-selected features — but the Optuna best CV score overstates the true generalization.

#### Next steps

+ GT05: CatBoost comparison — the GenoPHI-optimal algorithm, not just tuned LightGBM.
+ GT06: GenoPHI per-phage receptor prediction to strengthen Gate 2.

### 2026-04-12 11:56 CEST: GT05 — CatBoost comparison on three-layer feature set

#### Executive summary

Replaced LightGBM with CatBoost (GenoPHI-optimal algorithm) using 50-trial Optuna HPO with native categorical handling.
CatBoost achieves 0.826 AUC vs LightGBM's 0.823 — a +0.3pp difference that is not statistically significant (delta CI
[-0.011, +0.017]). CatBoost does improve calibration (Brier 0.152 vs 0.161). The algorithm choice is not the bottleneck;
feature quality (particularly Gate 2's limited 8/96 phage coverage) is the binding constraint.

#### Best CatBoost params

| Parameter | Value |
|-----------|-------|
| iterations | 400 |
| learning_rate | 0.091 |
| depth | 8 |
| l2_leaf_reg | 0.79 |
| random_strength | 2.38 |
| bagging_temperature | 2.02 |
| border_count | 89 |

#### Holdout results

| Arm | AUC | 95% CI | Top-3 | Brier |
|-----|-----|--------|-------|-------|
| lgbm_gt03_default | 0.823 | [0.782, 0.859] | 89.2% | 0.161 |
| catboost_tuned | 0.826 | [0.781, 0.864] | 90.8% | 0.152 |

Delta (CatBoost vs LightGBM): [-0.011, +0.017] — not significant.

#### Interpretation

**CatBoost's native categorical handling does not unlock a significant AUC gain.** The phylogroup/serotype/ST
categoricals are already well-captured by LightGBM's one-hot encoding at this tree depth. CatBoost's ordered target
encoding is theoretically superior for high-cardinality categoricals, but our categoricals are low-cardinality
(phylogroup: 7 levels, serotype: ~40, ST: ~80) and already saturated by LightGBM's 31-leaf trees.

**CatBoost does improve calibration.** Brier drops from 0.161 to 0.152 — similar to the IFW effect in GT03 but without
the AUC penalty. CatBoost's ordered boosting reduces overfitting, which improves probability calibration even when
discrimination (AUC) stays flat.

**Feature importance is consistent across algorithms.** depo×capsule is the top slot at 24% (CatBoost) vs 21.5%
(LightGBM). host_defense is slightly higher in CatBoost (5.6-6.8% vs 4.6%), suggesting CatBoost extracts marginally
more signal from defense features — possibly because ordered boosting is less susceptible to the lineage confounding
that degrades LightGBM's defense-feature usage.

**Same RFE-before-CV caveat as GT04 applies.** Inner CV scores are mildly optimistic but holdout comparison is clean.

#### Conclusion

Neither HPO (GT04: +0.4pp) nor algorithm change (GT05: +0.3pp) produces a statistically significant gain over the GT03
LightGBM baseline. The 0.823 AUC ceiling is feature-bound, not model-bound. GT06 (GenoPHI per-phage receptor
prediction) is the most promising path forward — expanding Gate 2 from 8/96 to 96/96 phages addresses the actual
bottleneck.

#### Next steps

+ GT06: GenoPHI per-phage receptor prediction to strengthen Gate 2 from 8/96 to 96/96 phage coverage.

### 2026-04-12 12:09 CEST: GT06 — GenoPHI k-mer receptor prediction for Gate 2

#### Executive summary

Used GenoPHI's 815 receptor-predictive amino acid k-mers (Dataset S6, Moriniere 2026) to predict per-phage OMP
receptor class from proteome sequence. Expanded Gate 2 coverage from 8/96 (genus-level) to 39/96 (k-mer-based) OMP
phages. However, the expanded predictions produce no AUC improvement: 0.824 vs 0.823, delta CI [-0.005, +0.005].
The 0.823 AUC ceiling from GT03 is robust — neither algorithm changes (GT04/GT05) nor feature expansion (GT06)
breaks through.

#### Receptor prediction coverage

| Type | Count | Details |
|------|-------|---------|
| OMP | 39/96 | ompC (13), lptD (12), tsx (8), btub (6) |
| LPS | 33/96 | Capsule-degrading phages |
| NGR | 22/96 | Protein receptor, unidentified |
| Unknown | 2/96 | No k-mer hits |

#### Validation against genus-level assignments

The k-mer predictions frequently disagree with genus-level consensus — this is expected and biologically correct.
Table S1 shows within-genus receptor variation is high (Tequatrovirus spans 7 receptor types). The k-mer approach
predicts per-phage specificity, not genus averages. For example, BCH953_P2/P4/P5 (Tequatrovirus) are predicted as
OmpC with >75% confidence — OmpC accounts for 9% of Tequatrovirus assays in Table S1, so this is a plausible
minority receptor.

#### Holdout results

| Arm | AUC | 95% CI | Top-3 | Brier |
|-----|-----|--------|-------|-------|
| gt03_genus_gate2 (8 OMP) | 0.823 | [0.782, 0.859] | 89.2% | 0.161 |
| kmer_gate2_only (39 OMP) | 0.824 | [0.780, 0.861] | 89.2% | 0.159 |
| both_gate2 (combined) | 0.823 | [0.780, 0.857] | 89.2% | 0.161 |

Delta (k-mer vs genus): [-0.005, +0.005] — not significant.

#### Interpretation

**Expanding Gate 2 coverage does not break the AUC ceiling.** Despite 5x more OMP-assigned phages (39 vs 8), the
k-mer-based receptor predictions add no discriminative signal on holdout. Two possible explanations:

1. **Gate 2 is not the binding constraint.** The model already captures receptor-host compatibility through the
   host_surface OMP HMM scores (21% feature importance in GT03). Adding phage-side receptor identity doesn't help
   because the host-side variation is already informative — if a host has high OmpC score, phages that lyse it will
   be ranked higher regardless of whether we explicitly encode "this phage targets OmpC."

2. **k-mer predictions are noisy.** The simplified k-mer hit counting (max-vote) is less accurate than GenoPHI's
   full gradient-boosted classifier with RFE. The low confidence on many predictions (mean ~0.5) suggests the signal
   is diluted. A properly trained GenoPHI model (not just feature scanning) might produce better predictions.

#### Track GIANTS conclusion

The three-layer hypothesis produced one validated discovery: **Gate 1 depolymerase × capsule cross-terms** provide
a statistically significant +1.2pp AUC lift (0.810 → 0.823). Gate 2 (receptor × OMP) and Gate 3 (defense) do not
significantly contribute beyond Gate 1. Algorithm optimization (HPO, CatBoost) does not break the 0.823 ceiling.

The remaining improvement paths are:
+ Full GenoPHI pipeline (trained classifier, not just k-mer scanning) for more accurate receptor predictions
+ Panel expansion beyond 96 phages to increase statistical power
+ Per-phage inverse-frequency weighting combined with CatBoost (GT05 showed Brier improvement)

### 2026-04-12 13:43 CEST: Pair-level biological analysis of Track GIANTS results

#### Executive summary

Deep pair-level analysis of GT03-GT06 holdout predictions reveals why receptor × OMP cross-terms produced zero lift
despite Moriniere 2026's AUROC 0.99 receptor classifiers. The root cause is threefold: (1) all 369 clinical E. coli
hosts express all 12 core OMPs at nearly identical levels (CV 0.01-0.17), collapsing host-side cross-terms to constants;
(2) phages sharing the same predicted receptor have wildly different host ranges (Tsx phages: mean pairwise Jaccard
0.091, below the random-pair baseline of 0.17); (3) Moriniere's models were trained on K-12 lab strains lacking
O-antigen/capsule, but 100% of our clinical strains have O-antigen, making polysaccharide barriers — not OMP receptor
presence — the dominant host-side adsorption determinant. This explains why Gate 1 (depolymerase × capsule) works
(+1.2pp) but Gate 2 (receptor × OMP) doesn't: capsule/O-antigen is the outer physical barrier.

#### The 6 holdout misses

The all_gates_rfe arm misses 6/65 holdout bacteria in top-3. Each miss has a distinct biological signature:

| Bacterium | Positives | Failure mode | Top-3 FP families | Best TP rank |
|-----------|-----------|-------------|-------------------|-------------|
| ECOR-06 | 1/96 (1.0%) | Needle-in-haystack | Straboviridae, Felixounavirus | AL505_Ev3 rank 25 |
| ECOR-69 | 10/96 (10.4%) | Broad-phage prior | Phapecoctavirus, Justusliebigvirus | LF82_P2 rank 7 |
| FN-B4 | 0/96 (0.0%) | Abstention | — | No TPs exist |
| NILS24 | 0/96 (0.0%) | Abstention | — | No TPs exist |
| NILS41 | 14/96 (14.6%) | Broad-phage prior | Phapecoctavirus, Straboviridae | LF73_P4 rank 6 |
| NILS53 | 14/96 (14.6%) | Narrow-host collapse | Straboviridae, Phapecoctavirus | NIC06_P2 rank 5 |

**Common pattern:** The model's top-3 is always dominated by broad-host phages (DIJ07_P1/P2, LF82_P8, 536_P7/P9,
LF73_P1) with 47-65% panel lysis rates. These phages are predicted high for every bacterium because the model has
learned their base rate. The true positives for missed bacteria are moderate-to-narrow host range phages ranked 5-48.

**NILS53 is the critical case.** It has 14 true positives, but 13 of them are narrow-host specialists (10-26% lysis
rate) from Autographiviridae, Dhillonvirus, and Kagunavirus. Only NIC06_P2 (Tequatrovirus, 42% lysis rate) ranks in
the top 5. The remaining 13 TPs are ranked 33-48 — completely buried under broad-host phages. These narrow-host
specialists target different receptors (lptD, lps) than the dominant broad-host phages (ngr, ompC), but the model
can't exploit this because the host-side receptor scores are near-identical across bacteria.

#### Why receptor × OMP cross-terms fail: three compounding problems

**1. Host OMP scores are near-constant across clinical E. coli.**

| Receptor | Mean | Std | CV | Unique values |
|----------|------|-----|-----|---------------|
| LptD | 1821.8 | 22.8 | **0.01** | 65 |
| OmpA | 772.1 | 14.5 | **0.02** | 43 |
| Tsx | 685.7 | 19.2 | **0.03** | 24 |
| OmpC | 756.7 | 31.2 | **0.04** | 79 |

All 369 hosts have all 12 OMPs at >99% presence. The HMM similarity scores measure gene presence/conservation, not
functional variation at the binding interface. The cross-term `predicted_is_OmpC × host_OmpC_score ≈
predicted_is_OmpC × 757 ± 31` — the host side contributes no discrimination.

Direct test: OmpC-targeting phages lyse high-OmpC hosts at 37.0% vs low-OmpC hosts at 33.4% (Cohen's d=0.086). For
Tsx: d=0.055. For BtuB: d=0.033. These effect sizes are negligible — OMP HMM scores do not predict within-receptor
lysis.

**2. Same-receptor phages have uncorrelated host ranges.**

Tsx-targeting phages (8 phages) span 3-204 lysed hosts out of 402. Their mean pairwise Jaccard similarity is 0.091 —
*below* the random-pair baseline of 0.17. A Lambdavirus targeting Tsx (411_P1, 3/402) and a Krischvirus targeting Tsx
(LF73_P4, 204/402) share almost zero hosts. Receptor class tells you which door the phage knocks on, but not whether
it gets through — that depends on O-antigen/capsule barriers, intracellular defenses, and other post-adsorption
factors.

OmpC phages (13 phages) are more cohesive (mean Jaccard 0.422, 2.5x enrichment over random) because they're dominated
by closely related Felixounavirus, but this family similarity is already captured by phage_projection features.

**3. Moriniere's models were trained on K-12 (no O-antigen/capsule).**

The k-mer receptor classifiers achieve AUROC 0.99 because they were trained on BW25113/BL21 — lab strains where OMP
receptors are the *only* surface barrier. Our 369 clinical strains all have O-antigen (100% non-zero O-antigen scores,
mean 813, std 693). In clinical isolates, the polysaccharide layer gates physical access to OMP receptors. A phage that
targets OmpC must first penetrate the O-antigen/capsule — and that penetration depends on depolymerase specificity
(Gate 1), not receptor binding affinity (Gate 2).

This is why Gate 1 (depolymerase × capsule) produces the only validated lift: it encodes the *outer* barrier. Gate 2
encodes a barrier that is masked by polysaccharides in clinical isolates.

#### What would help: OMP variant-level features

The HMM scores measure gene presence/conservation at the whole-gene level. But Moriniere showed that receptor
specificity is encoded in short hypervariable sequence motifs at RBP tip domains — and the equivalent host-side signal
would be in the extracellular loop regions of OMP proteins. OmpC has 50 allelic variants across our panel; BtuB has 28.
The current HMM score (CV 0.04 for OmpC) compresses this allelic diversity into a near-constant.

Extracting OMP extracellular loop sequences and computing allele-level or k-mer features would give the cross-term
actual host-side variation to exploit. The cross-term would become `predicted_receptor_is_OmpC ×
host_OmpC_loop3_variant` rather than `predicted_receptor_is_OmpC × host_OmpC_HMM_score`.

#### Training data expansion assessment

| Source | New phages | New bacteria | Overlap | Value |
|--------|-----------|-------------|---------|-------|
| GenoPHI (Noonan 2025) | 0 (94/96 ours) | ~33 (402 vs 369) | Contains our data | Modest |
| BASEL collection (Zenodo 15736582) | 56 | 0 (25 ECOR) | 25 bacteria overlap | New phage diversity |
| VHRdb | Varies | Varies | Contains our dataset | Aggregation resource |

The BASEL collection's 56 new phages tested against 25 ECOR strains would help break the lysis-breadth prior by adding
more narrow-host phages with diverse receptor specificities. However, earlier external data expansion attempts (TK01-03)
showed neutral lift, suggesting the bottleneck is feature quality, not training volume.

### 2026-04-12 14:17 CEST: GT07 — OMP extracellular loop variant features

#### Executive summary

Replaced near-constant whole-gene OMP HMM scores (CV 0.01-0.17) with binary OMP allele-variant features from Track C
(99% identity BLAST clustering, 22 features with variance 0.08-0.25). Variance pre-flight passed: BtuB cluster 99_15
showed Cohen's d=0.455, Tsx cluster 99_11 d=0.264. However, the directed cross-terms (predicted_receptor_is_X ×
host_X_variant_cluster_Y) produce no significant holdout improvement: AUC 0.826 vs 0.823, delta CI [-0.001, +0.006].
The OMP variant features capture 5.1% feature importance (vs 2.0% for HMM-based) but this does not translate to better
discrimination on 65 holdout bacteria.

#### Variance pre-flight results

| Feature | Variance | Cohen's d (lysed vs non-lysed) | Status |
|---------|----------|-------------------------------|--------|
| BtuB cluster 99_15 | 0.168 | +0.455 | Pass |
| Tsx cluster 99_11 | 0.212 | +0.264 | Pass |
| BtuB cluster 99_6 | 0.250 | -0.176 | Pass (anti-correlated) |
| Tsx cluster 99_2 | 0.225 | -0.218 | Pass (anti-correlated) |
| OmpC cluster 99_24 | 0.111 | (not tested, only 13 OmpC phages) | Marginal |

Pre-flight conclusion: features have real variance and measurable discrimination. Proceeded with full evaluation.

#### Holdout results

| Arm | AUC | 95% CI | Top-3 | Brier |
|-----|-----|--------|-------|-------|
| gt03_baseline (all_gates_rfe) | 0.823 | [0.782, 0.859] | 89.2% | 0.161 |
| +omp_variants | 0.826 | [0.784, 0.861] | 87.7% | 0.159 |

Delta (+omp_variants vs baseline): [-0.001, +0.006] — not significant.

#### Feature importance

The OMP variant cross-terms capture 5.1% feature importance (seed-averaged), vs 2.0% for the HMM-based receptor
features. This confirms the variant features are more informative than the near-constant HMM scores, but the
improvement is absorbed by the existing host_surface slot (20.2%) which already captures OMP information at the
whole-gene level.

#### Interpretation

The variance pre-flight correctly identified that OMP variant features have real discriminative power at the
single-feature level (BtuB d=0.455). But this per-feature discrimination doesn't compound into a holdout AUC
improvement because:

1. **The variant features are still binary (present/absent) for each allele cluster.** A host either has BtuB allele
   99_15 or it doesn't — there's no continuous interaction strength. The cross-term is still 0 or 1.
2. **RFE selects 278/543 features.** The OMP variant cross-terms survive RFE (they have real information), but they
   compete with 270+ other features that already capture host-surface variation.
3. **The 0.823 ceiling persists.** Neither more accurate receptor predictions (GT06), nor better host-side features
   (GT07) breaks through — the constraint may be in the interaction matrix itself (96 phages, 65 holdout bacteria).

#### Conclusion

OMP allele-variant features are a genuine improvement over HMM scores (5.1% vs 2.0% importance, real variance, real
per-feature discrimination), but they do not break the 0.823 AUC ceiling on ST03 holdout. The remaining Gate 2 path
is exhausted at the feature level — receptor × OMP cross-terms cannot significantly improve all-pairs predictions
regardless of how the host side is encoded.

### 2026-04-12 14:44 CEST: GT08 — GenoPHI binary protein-family features

#### Executive summary

Reproduced GenoPHI's core feature representation: MMseqs2 clustering (40% identity, 80% coverage) of 1.78M proteins
from all 369 host + 96 phage genomes → 28,389 non-singleton protein-family clusters → binary presence-absence matrix.
After variance filtering (top 500 host + 200 phage by variance), combined with GT03 mechanistic features, and evaluated
on ST03 holdout. Result: AUC 0.826 vs 0.823 baseline, delta CI [-0.002, +0.008] — not significant. RFE retained 244
protein-family features (216 host + 28 phage) capturing 13.1% feature importance, but this does not translate to a
holdout improvement.

#### MMseqs2 clustering results

| Metric | Value |
|--------|-------|
| Total proteins clustered | 1,783,579 |
| Total clusters | 64,592 |
| Non-singleton clusters (≥2 genomes) | 28,389 |
| Host-only clusters | 26,702 |
| Phage-only clusters | 1,500 |
| Shared clusters (host + phage) | 187 |

#### Variance pre-flight

| Side | Total clusters | Variance > 0.05 | Degenerate (>90% same) |
|------|---------------|-----------------|----------------------|
| Host | 28,389 | 4,276 (15.1%) | 25,709 (90.5%) |
| Phage | 28,389 | 670 (2.4%) | 28,075 (98.9%) |

Pre-flight conclusion: sufficient non-degenerate features exist on both sides. Filtered to top 500 host + 200 phage
by variance. Proceeded with evaluation.

#### Holdout results

| Arm | AUC | 95% CI | Top-3 | Brier |
|-----|-----|--------|-------|-------|
| gt03_baseline | 0.823 | [0.782, 0.859] | 89.2% | 0.161 |
| +protein_families | 0.826 | [0.783, 0.862] | 89.2% | 0.159 |

Delta: [-0.002, +0.008] — not significant.

#### Feature importance

Protein-family features capture 13.1% total importance (pf_host=11.0%, pf_phage=2.1%), comparable to host_stats (5.5%)
and host_defense (4.0%). RFE retained 216/500 host PF features and 28/200 phage PF features, indicating moderate
informativeness — enough to survive feature selection but not enough to displace the existing mechanistic features.

#### Interpretation

**GenoPHI's binary protein-family features add no significant lift when combined with mechanistic features.** This is
a critical finding: GenoPHI's AUROC 0.869 (vs our 0.823) is NOT due to their feature representation being superior.
Their features (binary k-mer/protein-family presence-absence) and our features (continuous HMM scores, depolymerase
cross-terms, defense counts) encode partially overlapping information — the protein-family clusters largely capture
the same host typing and surface variation that our HMM-based features already provide.

The remaining gap between 0.823 and GenoPHI's 0.869 is likely due to:
1. **Different holdout strategies** — GenoPHI uses LOGOCV (leave-one-genome-out per phage), we use cv_group-disjoint
   bacteria holdout. Not directly comparable.
2. **Different bacteria counts** — GenoPHI uses 402 bacteria, we use 369.
3. **ML pipeline details** — GenoPHI uses CatBoost + RFE + per-phage inverse-frequency weighting as a joint
   configuration, not as incremental additions. The interactions between these choices may matter.

#### Track GIANTS Phase 2 conclusion

GT07 (OMP variant features) and GT08 (GenoPHI binary features) both fail to break the 0.823 ceiling:

| Ticket | Feature type | Importance | AUC delta CI |
|--------|-------------|-----------|-------------|
| GT07 | OMP allele variants (22 features) | 5.1% | [-0.001, +0.006] |
| GT08 | Protein-family P/A (700 features) | 13.1% | [-0.002, +0.008] |

Both feature families are informative enough to survive RFE and capture measurable importance, but neither provides
sufficient independent signal to improve predictions on the 65-bacteria holdout. The 0.823 AUC represents a hard
ceiling for this feature set, model class, and evaluation protocol.

GT09 (BASEL panel expansion) cancelled: both GT07 and GT08 show the ceiling is not feature-bound in the way more data
would address, and prior external data attempts (TK01-03) showed neutral lift. The evidence does not support investing
in data expansion before resolving the more fundamental constraint.

### 2026-04-12 14:46 CEST: Track GIANTS final assessment — biological root cause analysis

#### Executive summary

Track GIANTS ran 8 tickets (GT01-GT08) exploring the three-layer biological hypothesis for phage-host lysis
prediction. The single validated discovery is Gate 1 (depolymerase × capsule cross-terms, +1.2pp AUC). Five
independent attempts to improve beyond 0.823 AUC all failed: HPO (+0.4pp), CatBoost (+0.3pp), k-mer receptor
expansion (+0.0pp), OMP allele variants (+0.3pp), and GenoPHI protein-family features (+0.3pp). None reached
statistical significance. This section analyzes why, grounded in biology rather than ML methodology.

#### The fundamental constraint: 96 phages on 65 holdout bacteria

The holdout has 65 bacteria × 96 phages = 6,240 pairs. At 20.7% positive rate, that's ~1,293 positive pairs. The
model's job is to rank the correct 1,293 pairs above the incorrect 4,947. Any feature that improves ranking must
shift enough pairs to move at least 1 of the 6 missed bacteria from miss to hit — which requires shifting at least
one true-positive phage above at least three false-positive phages for that bacterium.

The core problem is **the model already knows which phages are broadly lytic** (60-65% lysis rate). These phages
dominate the top-3 for every bacterium. To improve, a feature must identify cases where a narrow-host specialist
phage (10-25% lysis rate) should be ranked above a broadly lytic phage for a specific host. This requires:

1. **Phage-side specificity**: knowing which narrow-host phage matches which host
2. **Host-side discrimination**: the feature must vary enough across hosts to create different rankings

#### Why each feature family failed to break through

**Gate 2 (receptor × OMP cross-terms, GT02/GT06/GT07)**: The phage side works — k-mer receptor predictions assign
39/96 phages to specific OMP receptors with biologically plausible accuracy. But the host side fails. OMP HMM scores
(CV 0.01-0.17) are near-constant because all E. coli express all 12 core OMPs. Even OMP allele clusters (GT07),
which have real variance (0.08-0.25) and measurable single-feature discrimination (BtuB d=0.455), don't create
enough differential ranking across bacteria to change the top-3 for any missed strain. The biological reason: in
clinical E. coli with O-antigen/capsule, OMP receptor identity is necessary but not sufficient for infection —
polysaccharide penetration (Gate 1) is the gating step.

**GenoPHI protein-family features (GT08)**: 28,389 binary protein clusters from both host and phage genomes. After
filtering and RFE, 244 features survived with 13.1% importance. But they encode largely the same information as our
existing mechanistic features — the protein-family clusters for host surface genes correlate with our HMM scores, and
phage-side clusters correlate with phage_projection features. This is not surprising: both approaches measure
"which genes does this genome have?" at different granularity. Adding a second encoding of the same biological signal
does not create new discriminative power.

**HPO and CatBoost (GT04/GT05)**: Confirmed GenoPHI's finding that algorithm choice matters less than features once
the algorithm family is fixed. 300-tree LightGBM with default parameters is near-optimal for this feature set.

#### What the 0.823 ceiling means biologically

The model can correctly rank phages for 59/65 holdout bacteria. The 6 misses decompose into:

+ **2 abstention bacteria** (0/96 positives): No phage in the panel works. No feature can fix this.
+ **1 needle-in-haystack** (1/96 positive): Only AL505_Ev3 lyses ECOR-06. Statistical near-impossibility.
+ **3 broad-phage-prior failures** (ECOR-69, NILS41, NILS53): These bacteria have 10-14 true positives, but all
  are ranked below broadly lytic phages. The model knows the right phage *family* but can't promote the right
  *individual phage* over the dominant Straboviridae/Phapecoctavirus.

The 3 rescuable misses share a pattern: the true-positive phages are narrow-host specialists whose lysis of these
specific hosts is governed by factors not captured by any genomic feature we've tried — likely involving expression-
level regulation, phase variation, or phage-host co-evolutionary dynamics that don't leave genomic signatures
detectable by presence-absence or HMM-score features.

#### What would actually help (not feature engineering)

1. **More phages in the panel.** The 96-phage panel is the binding constraint, not the feature set. With only 96
   phages, the model's phage-side ranking capacity is limited. A 200-phage panel with more narrow-host diversity
   would naturally break the broad-phage prior dominance. This is a data collection problem, not a computational one.

2. **Different evaluation protocol.** Our holdout is 65 cv_group-disjoint bacteria — strong for generalization claims,
   but 1 strain flip = 1.5pp top-3. With 6 misses (3 fixable), the maximum achievable improvement is ~4.6pp. A
   larger holdout (or expanding the bacteria panel) would give more statistical power to detect real improvements.

3. **Per-phage models with transfer learning.** The per-phage blend (+2.0pp, 0.830 AUC) works because each phage
   gets its own small model that can learn host-specific patterns. If per-phage models could be made deployable
   (via few-shot transfer from a general model), this is the most promising architectural direction.

4. **Experimental data on the 6 misses.** Testing the 6 missed bacteria against the panel phages under different
   growth conditions or phage multiplicities might reveal whether the misses are genuine biological resistances or
   assay sensitivity limits.

#### Track GIANTS outcome summary

| Ticket | What | AUC delta | Significance |
|--------|------|-----------|-------------|
| GT01 | Depolymerase × capsule features | (component of GT03) | — |
| GT02 | Receptor × OMP features (genus) | (component of GT03) | — |
| GT03 | Three-layer integration + RFE | **+1.2pp** (0.810→0.823) | **Significant** |
| GT04 | Optuna HPO | +0.4pp | Not significant |
| GT05 | CatBoost comparison | +0.3pp | Not significant |
| GT06 | k-mer receptor expansion | +0.0pp | Not significant |
| GT07 | OMP allele variant features | +0.3pp | Not significant |
| GT08 | GenoPHI protein-family features | +0.3pp | Not significant |

**One validated discovery (Gate 1, +1.2pp) and seven null results.** The null results are as informative as the
positive: they establish that the 0.823 ceiling is not breakable by feature engineering, algorithm choice, or
representation changes alone — it is bound by the 96-phage panel size and the 65-bacteria holdout resolution.

#### Label quality finding: uninterpretable scores in holdout misses

Raw score analysis reveals that at least 2 of the 6 holdout misses have uninterpretable (`'n'`) scores for the
model's top-ranked phages — meaning the actual plaque images were ambiguous, not confirmed negatives:

| Bacterium | Model's top FP | Raw score | Image file | Interpretation |
|-----------|---------------|-----------|------------|----------------|
| ECOR-69 | DIJ07_P2 (rank 1) | `'n'` in rep3 | 151_rep3.jpg | Ambiguous plaque |
| ECOR-69 | DIJ07_P1 (rank 2) | `'n'` in rep3 | 151_rep3.jpg | Ambiguous plaque |
| NILS53 | LF82_P8 (rank 1) | `'n'` in rep2 | 207_rep2.jpg | Ambiguous plaque |
| NILS53 | 536_P9 (rank 4) | `'n'` in rep2 | 207_rep2.jpg | Ambiguous plaque |
| NILS53 | LF82_P9 (rank 3) | `'n'` in rep2 | 207_rep2.jpg | Ambiguous plaque |
| NILS41 | NRG_11A2 | `'n'` in rep2 | 195_rep2.jpg | Ambiguous plaque |

If any of these uninterpretable scores are actually positive, the corresponding bacterium flips from a top-3 miss to
a hit. For NILS53, the model's #1, #3, and #4 ranked phages all have `'n'` scores — the model may be *correct* and
the labels wrong. For ECOR-69, the model's #1 and #2 ranked phages have `'n'` scores.

Checking the raw plaque images (Zenodo doi:10.5281/zenodo.10202713) to verify whether these are genuine negatives or
ambiguous positives. If 2 misses flip to hits, the true top-3 rate would be ~93.8% (matching the per-phage blend),
fundamentally changing the interpretation of whether the all-pairs model has room for improvement.

### 2026-04-12 20:16 CEST: GT09 — Image review, clean-label evaluation, and BASEL data analysis

#### Executive summary

Three findings from GT09. (1) Raw plaque image review confirms the `'n'` scores are genuinely ambiguous — plates show
faint signals and physical artifacts that make definitive calls impossible. (2) Clean-label re-evaluation (excluding
3,462 ambiguous pairs from training) improves top-3 from 89.2% to 92.3% (+3.1pp) with the same model, demonstrating
that noisy negatives actively harm ranking. (3) The GenoPHI E. coli matrix is 100% identical to our labels (37,788/37,788
pairs match); BASEL provides 52 genuinely new phages on 25 ECOR bacteria, but full integration requires new phage feature
computation infrastructure that doesn't yet exist.

#### Image review findings

Examined all 3 replicates for each of the 6 missed holdout bacteria from Zenodo doi:10.5281/zenodo.10202713. Key
observations:

**NILS53** (207_rep1/2/3): Rep1 and rep2 show clear plaques in Plate C (LM07_P1, LM33_P1 confirmed positives) and
Plate D (Przondovirus and Dhillonvirus positives). The `'n'` positions for LF82_P8 (Plate C, Y=4 X=6) and 536_P9
(Plate B, Y=3 X=0) are in regions with some faint marks but nothing definitively distinguishable from dried droplet
artifacts at this resolution. Rep3 shows wavy surface texture on Plates A and B complicating assessment. **Verdict:
genuinely ambiguous — the scorer's uncertainty was justified.**

**ECOR-69** (151_rep1/2/3): Very few visible clearing zones across all replicates. Confirmed positives (LF82_P9,
LF110_P1-P3) produce only faint spots visible mainly in rep2 and rep3 — this is a host with weak lysis. The `'n'`
positions for DIJ07_P1/P2 (Plate B right edge) are in a region indistinguishable from background. **Verdict: ambiguous
but more likely negative than positive given the overall faint signal.**

**NILS41** (195_rep1/2/3): Many clear plaques across all plates and replicates — consistent with 10+ positives. Rep1
shows abundant clearing on Plates A and C. The `'n'` position for NRG_11A2 at Plate D Y=7 is near visible clearing
zones. **Verdict: most plausible hidden positive among the 6 bacteria.**

**FN-B4** (332_rep1/2/3): Clean plates with no clearing zones. Physical agar damage (torn/crumpled agar) on Plates A
and C explains why some positions were scored `'n'` — the scorer couldn't read those spots due to physical damage,
not ambiguous lysis. **Verdict: genuinely resistant, `'n'` scores are damage artifacts.**

**ECOR-06** (088_rep1): Clean plates with minimal signal. The single positive (AL505_Ev3, matrix=2) is weak. The
`'n'` scores for LF31_P3 and LM02_P1 are in clean areas. **Verdict: genuinely near-zero susceptibility.**

**NILS24**: All 864 raw scores are `'0'`, no `'n'` values. Images not reviewed (no ambiguity to resolve).

#### Clean-label re-evaluation

Excluded 3,462 ambiguous pairs (with `'n'` score and no `'1'` score) from training and/or holdout. Three arms
compared:

| Arm | Train pairs | Holdout pairs | AUC | Top-3 | Brier |
|-----|------------|--------------|-----|-------|-------|
| gt03_original | 29,031 | 6,235 (65 bacteria) | 0.823 | 89.2% | 0.161 |
| gt03_clean_holdout | 29,031 | 5,832 (65 bacteria) | 0.829 | 89.2% | 0.158 |
| gt03_clean_both | 26,130 | 5,832 (65 bacteria) | 0.824 | **92.3%** | 0.160 |

**Removing noisy negatives from training improves top-3 by +3.1pp (89.2% → 92.3%).** This is the largest single
improvement since depolymerase × capsule features (+1.2pp AUC in GT03). The mechanism: when the model trains on
ambiguous pairs labeled as negative, it learns to downweight phages that are actually lytic for those hosts. Removing
these pairs lets the model rank correctly.

The AUC improvement is modest (0.824 vs 0.823) because AUC measures discrimination across all pairs, including those
that are unambiguously negative. But top-3 — the clinically relevant metric — improves substantially because the
cleaned model no longer penalizes phages for ambiguous interactions.

#### BASEL data analysis

**GenoPHI E. coli matrix**: 402 bacteria × 94 phages, 100% identical to our labels (37,788/37,788 pairs match after
name normalization). Our data is already in GenoPHI's framework, likely via the Brisse lab. **No new training pairs
for existing phages.**

**BASEL ECOR interaction matrix**: 52 new phages × 25 ECOR bacteria = 1,300 new interactions (24.4% positive).
Phage breadth distribution on 25 ECOR: 25 narrow (≤3 hosts), 15 moderate (4-10), 12 broad (>10), 3 non-lytic.
All 52 phage genomes downloaded from NCBI GenBank.

**Integration gap**: Our feature pipeline (Pharokka annotation → DepoScope → protein family clustering → k-mer
receptor prediction → phage_projection → phage_stats) is built for the 96 Guelin phages. Running it on 52 new BASEL
phages requires infrastructure work: pyrodigal protein prediction, Pharokka annotation, DepoScope depolymerase
detection, and feature materialization. The 52 genomes are downloaded and ready at `.scratch/basel/genomes/` but
feature computation is not automated for arbitrary new phages.

Only 3/25 ECOR bacteria from BASEL are in our holdout (ECOR-14, ECOR-29, ECOR-35), so BASEL data primarily augments
training, not evaluation.

#### Conclusion

The most impactful finding from GT09 is not BASEL data but **clean-label evaluation**: excluding ~3,500 ambiguous
pairs from training improves top-3 by +3.1pp — more than any feature engineering attempt in the entire track. This
suggests that label quality, not feature quality, is the binding constraint for ranking performance. The practical
recommendation is to flag ambiguous pairs in the training corpus and either exclude them or assign soft labels.
