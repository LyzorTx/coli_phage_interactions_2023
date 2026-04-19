# Project Knowledge Model

<!-- Last consolidated: 2026-04-19T16:00:00+02:00 -->
<!-- Source: lyzortx/research_notes/lab_notebooks -->

**65 knowledge units** across 7 themes (46 active, 19 dead ends)

## Data & Labels

Labeling policy, data quality, split contracts, and training corpus.

- **`label-policy-binary`**: CHISEL label frame: the atomic training unit is (bacterium, phage, concentration,
  replicate) → {0, 1} — each raw observation in raw_interactions.csv stands on its own. Concentration enters the model
  as a numeric feature (log_dilution ∈ {0, -1, -2, -4}); rows with score='n' are dropped from training, not silently
  treated as negative. This supersedes the SPANDEX-era any_lysis rollup, which collapsed all concentrations and
  replicates for a pair into a single binary label and absorbed ~7% single-hit noise (matrix_score=0 despite positive
  any_lysis) and ~10% 'n'-as-negative ambiguity into the training set. [validated; source: ST0.1, 2026-04-18 CHISEL
  design; see also: raw-interactions-authority, mlc-dilution-potency, ambiguous-label-noise]
  - *Rationale: the any_lysis rollup was a derived summary imposed on top of the raw data. Training on the raw
    observations naturally integrates BASEL (single concentration → one row) and any future dilution-series panels
    without a rollup convention. The SPANDEX label-policy-binary statement (any_lysis + matrix_score noise) describes
    how SX10 and earlier models were trained and stays valid for historical comparison; CH02 keeps any_lysis semantics
    as a regression safety-net before the CH04 behavioral flip.*
- **`split-contract`**: ST02 defines 294 training bacteria; ST03 reserves 65 cv_group-disjoint bacteria for sealed
  holdout. Deterministic cv_group hashing with salt ensures reproducibility. [validated; source: ST0.2, ST0.3, AR01; see
  also: raw-interactions-authority, cv-group-leakage-fixed]
- **`cv-group-leakage-fixed`**: CHISEL CH02 fixed a cross-validation fold-hashing bug in
  lyzortx/pipeline/autoresearch/sx01_eval.assign_bacteria_folds. Folds were hashed on bacterium name instead of
  cv_group, so 45 of 48 multi-bacterium cv_groups (≤1e-4 ANI clusters) split across folds — near-duplicate bacteria
  appeared on both sides of the train/test boundary. The fix hashes on cv_group: all bacteria in a cv_group land in the
  same fold. Revalidating SX10 canonical (GT03 all_gates_rfe + AX02 per-phage blending, any_lysis labels, SX10 feature
  bundle, 10-fold bacteria-axis CV, 369×96 panel) under the corrected scheme gave AUC 0.8521 [0.8381, 0.8649] and Brier
  0.1317 [0.1253, 0.1381] — a −1.78 pp AUC and +0.69 pp Brier shift vs the SPANDEX SX10 numbers (AUC 0.8699, Brier
  0.1248). All SPANDEX-era aggregates that went through assign_bacteria_folds (SX03/SX04/SX10/SX11/SX12/SX13/SX14/SX15)
  inherit a comparable small positive bias; their per-arm null conclusions remain valid because the leakage was uniform
  across arms, but headline numbers are subject to this shift. [validated; source: CH02, 2026-04-19 CHISEL CV fix; see
  also: split-contract, spandex-final-baseline, spandex-unified-kfold-baseline]
  - *The AUC drop is larger than the 17%-of-bacteria leakage footprint alone would predict because changing the hash
    input from bacterium name to cv_group id reshuffles every singleton cv_group's fold too (323 / 369 bacteria change
    fold, not just the ~70 in split cv_groups). The shift therefore captures both leakage correction and
    fold-composition variance on singletons. The Brier deterioration with disjoint CIs is the cleanest signal:
    calibration is genuinely worse once near-duplicates are prevented from leaking. The new numbers are not yet
    canonical — CH03 must first verify row-expansion preserves any_lysis semantics (safety-net) before CH04 flips to
    per-row binary training. Canonical artifacts: lyzortx/generated_outputs/ch02_cv_group_fix/sx10_revalidated/
    (predictions + bootstrap), ch02_sx10_revalidated_metrics.json (summary + fold diff), ch02_fixed_folds.csv (bacterium
    → fold mapping under CHISEL hashing).*
- **`raw-interactions-authority`**: raw_interactions.csv is the authoritative training corpus; all derived training
  cohorts and evaluation splits trace back to this file. [validated; source: ST0.2, AR01, DEPLOY01]
- **`mlc-dilution-potency`**: HISTORICAL (SPANDEX-era). MLC is a rule-based rollup over the raw (bacterium, phage,
  concentration, replicate) → {0, 1} observations, not a raw field in raw_interactions.csv. The paper's MLC ∈ {0..4}
  distinguishes plaque morphology at 5x10^6 pfu/ml (individual plaques vs lawn lysis); our binary spot data cannot
  reconstruct that split, so SX05 capped the SPANDEX label at {0..3}. CHISEL abandons MLC as a training/evaluation label
  (see label-policy-binary); this unit remains only to document the SPANDEX-era convention. [validated; source:
  Gaborieau 2024 Methods, Gaborieau 2024 Fig 2b legend, 2026-04-13 raw CSV verification; see also: label-policy-binary,
  raw-interactions-authority, ambiguous-label-noise]
- **`basel-binary-only`**: BASEL × ECOR interaction data (52 phages × 25 ECOR bacteria from GenoPHI) is confirmed binary
  only — no graded upstream data exists. BASEL used a single high-titer spot test (>10^9 pfu/ml), not a dilution series.
  BASEL EOP data on Zenodo covers K-12 defense experiments only, not the ECOR host range panel. [validated; source:
  2026-04-12 SPANDEX design, Maffei 2021, Maffei 2025, Noonan 2025; see also: mlc-dilution-potency,
  external-data-neutral, genophi-data-identical]
  - *BASEL spot test at >10^9 pfu/ml is ~2x higher than our maximum concentration (~5x10^8). BASEL positive maps to MLC
    >= 1 (lysis at high titer confirmed). BASEL negative maps to high-confidence MLC=0 (failed at even higher titer than
    our max). 1,240 observed pairs (302 positive, 24.4%), 160 missing. The matrix is already sparse — not all 52 phages
    tested on all 25 bacteria.*

## Features & Feature Engineering

What works, what doesn't, leakage risks, and encoding decisions.

- **`host-trait-field-effect`**: Phylogroup, ST, and serotype structure hard-to-lyse behavior at field level
  (Kruskal-Wallis p=1.24e-09) but no individual trait values survive FDR correction. [preliminary; source: TB03]
- **`narrow-susceptibility-rescue`**: Narrow-susceptibility rescue is concentrated in 19/96 phages (19.8%); Myoviridae
  dominate, and no single universal rescuer exists. [validated; source: TB04; see also: family-bias-straboviridae]
- **`receptor-variant-richness`**: Receptors are richer than binary presence: OmpC spans 50 variants, BtuB spans 28. The
  TL18 pipeline discarded this structure via binary thresholds; the AUTORESEARCH pipeline uses continuous HMM scores
  (e.g., OmpC: 79 unique values in [509, 818]). [validated; source: TC01, 2026-04-09 AUTORESEARCH host_surface audit;
  see also: defense-integer-counts, receptor-specificity-solved]
- **`rbp-receptor-unknowns`**: Only 77/96 phages have curated receptor data; 19 are explicitly uncovered rather than
  guessed, requiring explicit unknown handling in RBP-receptor compatibility features. [validated; source: TE01; see
  also: receptor-specificity-solved]
  - *Moriniere et al. 2026 can predict receptor class from phage proteome (AUROC 0.99) via GenoPHI, which could fill
    this gap for the 19 uncovered phages.*
- **`defense-integer-counts`**: Defense system gene counts reflect real biological redundancy (e.g., 2 vs 1 MazEF
  copies); HMM detection scores reflect only tool confidence. Integer counts are the correct encoding. [validated;
  source: DEPLOY06, DEPLOY track design]
- **`defense-lineage-confounding`**: Defense subtypes correlate with phylogroup, causing lineage confounding: defense
  features rerank borderline phages by host lineage rather than mechanistic evasion signal. [validated; source:
  2026-04-08 defense ablation, 2026-04-08 defense top-3 deep dive; see also: defense-ablation-autoresearch,
  adsorption-first-strategy, ml-pipeline-over-features]
  - *Specific mechanism: lytic phages dropped from top-3 when the host defense profile resembled phylogroups where those
    phages don't lyse. The top-3 regression only manifests in multi-seed aggregated predictions, not individual seeds.
    GenoPHI's RFE feature selection may automatically mitigate this by excluding confounded features.*
- **`adsorption-dominates-paper`**: Adsorption factors dominate lysis prediction across independent studies: Gaborieau
  2024 found 30 significant adsorption traits vs 2 defense; Noonan/GenoPHI 2025 SHAP analysis replicates this ratio;
  Boeckaerts 2024 found no defense correlation for Klebsiella. [validated; source: 2026-04-08 paper analysis, Nature
  Microbiology 2024, Noonan 2025 GenoPHI; see also: adsorption-first-strategy, defense-lineage-confounding]
- **`receptor-specificity-solved`**: Phage receptor specificity is a near-solved discrete classification problem: 19
  receptor classes for E. coli, predicted from amino acid 5-mer presence-absence on the phage proteome at median AUROC
  0.99 with zero false positives on independent validation. [validated; source: Moriniere 2026, 2026-04-11 literature
  review; see also: receptor-variant-richness, pairwise-cross-terms-dead-end, plm-rbp-redundant, rbp-receptor-unknowns]
  - *Receptor specificity is encoded in localized hypervariable sequence domains (HVSs) in RBP tip regions, not in
    diffuse structural features. Single amino acid changes can switch receptor class (Q206L: OmpF to OmpW). This makes
    general protein embeddings (PLM, mean-pooled) the wrong tool — the signal is in short motifs at specific loci, not
    global protein similarity. Limitation: training on BW25113/BL21 which lack O-antigen/capsule, so
    polysaccharide-mediated specificity is not covered.*
- **`three-layer-hypothesis`**: Lysis requires passing adsorption gates (capsule penetration OR receptor binding) then
  surviving host defenses; features should be directed per-layer, not flat. Gate 1 (depolymerase × capsule) is validated
  (+1.2pp AUC); Gate 2 (receptor × OMP) fails due to host-side OMP score homogeneity; Gate 3 (defense) is marginal.
  [validated; source: Moriniere 2026, Noonan 2025 GenoPHI, 2026-04-11 literature review, GT03; see also:
  receptor-specificity-solved, defense-lineage-confounding, adsorption-dominates-paper, depo-capsule-validated,
  omp-score-homogeneity]
  - *GT03 ablation on ST03 holdout: Gate 1 alone +0.7pp AUC (significant), Gate 2 alone +0.3pp (not significant), Gate 3
    alone +0.6pp (not significant), all gates with RFE +1.2pp AUC (significant, 0.810→0.823). Gate 1 works because
    capsule profiles vary richly across clinical strains (99 features, all nonzero). Gate 2 fails because all E. coli
    have all 12 core OMPs at near-identical HMM scores (CV 0.01-0.17). The polysaccharide layer gates physical access to
    OMP receptors in clinical isolates — Gate 1 is the outer barrier.*
- **`host-capsule-variation`**: Our 369 clinical E. coli hosts have rich capsule variation (99 capsule feature columns,
  all nonzero) — they are NOT K-12 lab strains and can support depolymerase-capsule learning. [validated; source:
  2026-04-11 host_surface audit; see also: three-layer-hypothesis, receptor-variant-richness, depo-capsule-validated]
- **`depo-capsule-validated`**: Depolymerase × capsule pairwise cross-terms (242 features: 41 DepoScope cluster
  memberships + has_depo/count × 99 capsule HMM scores) provide a statistically significant +1.2pp AUC lift
  (0.810→0.823) with 22% feature importance — the only validated pairwise feature family. [validated; source: GT03; see
  also: three-layer-hypothesis, host-capsule-variation, omp-score-homogeneity]
  - *Gate 1 alone gives +0.7pp AUC (significant). Combined with RFE pruning (502→252 features), the three-layer model
    reaches 0.823 AUC [0.781, 0.858]. The cross-terms capture a specific mechanism: this phage has a depolymerase that
    can degrade this host's capsule type. This works because capsule profiles vary richly across clinical E. coli (99
    features with distinct variation), unlike OMP receptor scores which are near-constant.*
- **`omp-score-homogeneity`**: Whole-gene OMP HMM scores are near-identical across 369 clinical E. coli (CV 0.01-0.17),
  making whole-gene receptor × OMP cross-terms collapse. Loop-level allelic variation is substantial (BTUB 28 MMseqs2
  99% clusters, OMPC 49, across 5546 retained 5-mers) but does not predict lysis either — see
  host-omp-variation-unpredictive. [validated; source: GT03, GT06, 2026-04-12 pair-level analysis, SX13; see also:
  receptor-variant-richness, receptor-specificity-solved, three-layer-hypothesis, host-omp-variation-unpredictive]
  - *LptD CV=0.01, OmpA CV=0.02, Tsx CV=0.03, OmpC CV=0.04 on HMM scores. Direct test: OmpC phages lyse high-OmpC hosts
    at 37.0% vs low-OmpC at 33.4% (Cohen d=0.086, negligible). SX13 tested whether loop-level k-mer or cluster
    representations would restore host-side signal — they did not (AUC deltas ≤ +0.3 pp across marginal, cross-term, and
    cluster arms). The homogeneity observation at the HMM level is real; the escalation to finer representations doesn't
    rescue prediction because host-range variance lives downstream of OMP recognition in clinical strains
    (polysaccharide access + post-adsorption factors).*
- **`same-receptor-uncorrelated-hosts`**: Phages sharing a predicted receptor class have weakly correlated or
  uncorrelated host ranges: Tsx phages (8 phages) have mean pairwise Jaccard 0.091 (below random-pair baseline of 0.17);
  only OmpC phages show above-random cohesion (Jaccard 0.42) due to Felixounavirus family dominance. [validated; source:
  GT06, 2026-04-12 pair-level analysis; see also: omp-score-homogeneity, receptor-specificity-solved]
  - *Receptor class tells which OMP a phage binds, but not whether it infects a given host. Post-adsorption factors
    (O-antigen/capsule blocking, intracellular defenses, injection efficiency) dominate host-range determination in
    clinical isolates. This means receptor identity is necessary but far from sufficient for predicting strain-level
    lysis.*

## Model Architecture & Performance

Architecture choices, calibration, and performance bounds.

- **`tl18-flawed-baseline`**: TL18 model (0.823 AUC) is not a valid baseline: DefenseFinder version drift inflated 17.3%
  of feature importance, and 5 soft-leaky pairwise features contributed ~5.5%. [validated; source: TL18 audit; see also:
  autoresearch-baseline]
- **`chisel-baseline`**: CHISEL canonical baseline (CH04, 2026-04-19): per-row binary training on every interpretable
  (bacterium, phage, log_dilution, replicate, X, Y) raw observation with `pair_concentration__log_dilution` as a numeric
  feature, SX10 feature bundle otherwise unchanged (host_surface + host_typing + host_stats + host_defense +
  phage_projection + phage_stats + pair_depo_capsule + pair_receptor_omp, RFE-selected), **all-pairs only** (AX02
  per-phage blending retired, see per-phage-retired-under-chisel), 10-fold bacteria-axis CV under CH02 cv_group hashing,
  369×96 panel. Evaluation: each held-out pair scored at its max observed log_dilution (all 35,266 pairs tested at
  log_dilution=0, so evaluation is at neat concentration) with bacterium-level bootstrap CIs (1000 resamples). Result:
  **AUC 0.8084 [0.7944, 0.8217], Brier 0.1750 [0.1677, 0.1824]**. This is the active reference point for all future
  CHISEL arms. [validated; source: CH04, 2026-04-19 CHISEL baseline; see also: spandex-final-baseline,
  cv-group-leakage-fixed, label-policy-binary, ranking-metrics-retired, per-phage-retired-under-chisel, deployment-goal]
  - *CH04 drops substantially from the CH02 revalidated baseline (AUC 0.8521, Brier 0.1317): ΔAUC = −4.37 pp, ΔBrier =
    +4.33 pp, with disjoint 95% CIs on both. Three changes compound between the two numbers: (a) per-row training
    replaces pair-level any_lysis — every (bacterium, phage, log_dilution, replicate, X, Y) raw observation with score ∈
    {0, 1} becomes its own training row (rows with score = "n" are dropped as missing, not negative); (b) log_dilution
    enters as a numeric feature so the model must predict at a specific concentration rather than "does this pair ever
    lyse?"; (c) AX02 per-phage blending is retired (see per-phage-retired-under-chisel), removing the bacterium-level
    memorization head that contributed to CH02's AUC under the bacteria-axis setup. Diagnostic decomposition: evaluation
    label distribution is nearly unchanged (27.4% positive at max-conc vs 27.6% pair-level any_lysis — ~47 pair labels
    flip), so the drop is training-side, not evaluation-side. The concentration feature is the #4 ranked feature by mean
    LightGBM importance (328.70, retained by RFE in all 30 fold × seed fits, above every receptor-OMP cross-term).
    CHISEL takes the 4.4 pp aggregate-AUC cost in exchange for a deployable per-observation predictor aligned with
    deployment-goal. Subsequent CHISEL tickets (CH05 phage-axis, CH06 both-axis, CH07 feature re-audit) are anchored to
    this baseline, not to SPANDEX numbers. Canonical artifacts:
    lyzortx/generated_outputs/ch04_chisel_baseline/ch04_aggregate_metrics.json, ch04_predictions.csv (pair-level),
    ch04_per_row_predictions.csv (per-row), and ch04_feature_importance.csv.*
- **`spandex-final-baseline`**: HISTORICAL (SPANDEX-era). Superseded as the active canonical by chisel-baseline (CH04).
  SPANDEX final configuration — GT03 all_gates_rfe + AX02 per-phage blending on SX05-corrected MLC 0-3 labels, 10-fold
  CV bacteria-axis on the 369×96 panel: **AUC 0.8699 [0.8570, 0.8819], Brier 0.1248 [0.1187, 0.1309]** within-panel,
  inflated by the cv_group fold-hashing bug (see cv-group-leakage-fixed). CH02 revalidation under fixed folds lands at
  **AUC 0.8521 [0.8381, 0.8649], Brier 0.1317 [0.1253, 0.1381]**. Cross-panel Arm C (train Guelin → predict BASEL ×
  ECOR): AUC 0.7607 [0.6886, 0.8307], Brier 0.1844 [0.1426, 0.2213]. Historical nDCG/mAP numbers (nDCG 0.7958, mAP
  0.7111 within-panel; nDCG 0.7619, mAP 0.5186 cross-panel) are retained for SPANDEX-era comparison but are no longer on
  the scorecard. SX07 and SX09 cancelled (plm-rbp-redundant); SX08 null. Wave-2 (SX11–SX13) also null under AUC+Brier —
  see ordinal-regression-not-better, kmer-receptor-expansion-neutral, and host-omp-variation-unpredictive. [validated;
  source: SX05, SX06, SX08, SX10; see also: chisel-baseline, autoresearch-baseline, mlc-dilution-potency,
  new-phage-generalization, plm-rbp-redundant, panel-size-ceiling, label-policy-binary, cv-group-leakage-fixed]
  - *The 10.9 pp AUC gap between within-panel (0.87) and cross-panel (0.76) is the load-bearing generalization problem
    inherited by CHISEL; closing it still requires panel expansion rather than richer phage-side features. MLC-derived
    label scoring and nDCG narratives are historical artifacts of the SPANDEX scorecard and do not constrain CHISEL.
    Per-family null conclusions remain valid — the leakage bug was a uniform inflation across arms, not an arm-selective
    effect. CHISEL's per-row training (CH04) lands at 4.4 pp lower AUC than CH02 revalidated; the SPANDEX within-panel
    number is not a ceiling CHISEL needs to re-hit, because the training unit and evaluation question are both
    different. Canonical SPANDEX artifacts: lyzortx/generated_outputs/sx05_sx01_eval/ (within-panel) and
    lyzortx/generated_outputs/sx06_sx03_eval/ (cross-panel Arm C).*
- **`stratified-eval-framework`**: CHISEL rescope: per-stratum reporting is a narrow-host diagnostic, not a primary
  scorecard. Under AUC+Brier, aggregate and stratified metrics are more aligned than under nDCG — the SPANDEX marquee
  case (SX11 LambdaRank +3.5 pp within-family nDCG, null in aggregate) does not survive the metric change. CHISEL
  retains the four-stratum decomposition (within-family / cross-family / narrow-host-phage / phylogroup-orphan) for
  error-bucket diagnosis, especially the narrow-host-phage stratum where broad-phage priors dominate the panel and
  rescue is the deployment question. It is no longer required reporting for every ticket. [validated; source: SX14,
  2026-04-18 CHISEL design; see also: panel-size-ceiling, spandex-unified-kfold-baseline, narrow-host-prior-collapse]
  - *The stratum definitions and bootstrap machinery are still useful diagnostics. Stratum routing at inference still
    requires computing training-positive family overlap and phylogroup sibling counts from the training fold (no test
    leakage). Future CHISEL tickets evaluating a single candidate arm can reuse
    `.agents/skills/case-by-case/compare_predictions.py` for per-bacterium audit and the sx14_eval.py pipeline for full
    four-stratum decomposition when narrow-host behaviour is specifically under investigation.*
- **`chisel-unified-kfold-baseline`**: CHISEL unified Guelin+BASEL k-fold baseline (CH05, 2026-04-19): per-row binary
  training on the unified 148-phage × 369-bacteria panel (36,643 pairs: 35,403 Guelin + 1,240 BASEL), SX10 feature
  bundle, all-pairs only (per-phage blending retired track-wide per `per-phage-retired-under-chisel`). Two axes:
  **bacteria-axis AUC 0.8061 [0.7917, 0.8199], Brier 0.1778 [0.1702, 0.1853]** (10-fold CH02 cv_group hash; all 148
  phages in training per fold); **phage-axis AUC 0.8850 [0.8617, 0.9062], Brier 0.1348 [0.1219, 0.1495]** (10-fold
  StratifiedKFold by ICTV family + "other" <10-phage bucket + "UNKNOWN" no-family bucket — 40% of folds are
  pseudo-family catch-alls; calling it "ICTV-stratified" without that qualifier misleads). This unit replaces the
  earlier "cross-source AUC parity" headline with three separate findings: **(1) phage-axis discrimination parity**
  (Guelin 0.8863 vs BASEL 0.8818, |ΔAUC| 0.0045 — a weak non-rejection on 52 BASEL phages with CI 3× wider than
  Guelin's, not positive evidence of transfer); **(2) phage-axis calibration divergence** (Guelin Brier 0.1329 vs BASEL
  0.1884, disjoint CIs, BASEL mid-P reliability gap 21-27 pp wider than Guelin's in the 0.5-0.9 predicted-probability
  bins); **(3) BASEL bacteria-axis deficit** (BASEL-only bacteria-axis AUC 0.7152 on the 1,240 BASEL pairs vs
  Guelin-only 0.8098 on the same axis — a 9.5 pp BASEL-specific deficit invisible in the 96.6% Guelin-weighted
  aggregate, a standalone deployability finding separate from the phage-axis parity story). Root-cause diagnostic
  (post-hoc, no model rerun): the mid-P miscalibration localises to the 39/52 BASEL phages whose `phage_projection`
  vectors are non-zero (Brier 0.31 bacteria-axis) — these phages map into Guelin-derived TL17 neighborhoods associated
  with broad-host lysis but carry narrower actual host ranges; the 13/52 BASEL phages with zero-vector projection
  calibrate correctly (Brier 0.12) because the model has no phage signal to misuse and falls back to the host-side
  prior. Straboviridae exclusion closes only 1.5 pp of the 9.5 pp bacteria-axis BASEL deficit — family bias is not the
  driver. This is the active CHISEL reference for two-axis generalization and cross-source behaviour. [validated;
  source: CH05, 2026-04-19 CHISEL unified k-fold; see also: chisel-baseline, spandex-unified-kfold-baseline,
  per-phage-retired-under-chisel, cv-group-leakage-fixed, new-phage-generalization, deployment-goal, plm-rbp-redundant,
  panel-size-ceiling]
  - *Phage-axis AUC exceeds bacteria-axis AUC by 7.9 pp (CIs disjoint). The gap is structural, not a deployment-value
    signal: phage-axis folds hold out entire phages but keep all 369 bacteria in training; bacteria-axis folds hold out
    bacteria and remove host-side signal for those test pairs. The SX15 "per-phage blending tax" framing (which
    interpreted this gap under per-phage blending enabled) is retired — under all-pairs-only evaluation there is no tax;
    the gap is purely structural training-data coverage. Bacteria-axis aggregate AUC 0.8061 matches CH04's 0.8084 within
    0.25 pp because adding 1,240 BASEL pairs to 35K Guelin pairs barely shifts the 96.6%-Guelin-weighted aggregate —
    which is exactly why the BASEL-specific 9.5 pp bacteria-axis deficit had to be reported separately. The earlier
    draft's "BASEL phages generalize essentially identically to Guelin phages" claim was a discrimination-only finding,
    not deployment readiness — retired. The TL17-bias root cause supersedes the earlier encoding-hypothesis framing:
    CH05's reliability analysis showed BASEL over-predicted in mid-P, the opposite direction an encoding mismatch
    predicts, so the pre-existing relative-log_dilution encoding of BASEL (at Guelin neat) was not the dominant driver;
    the encoding has been fixed track-wide to absolute log10_pfu_ml (Guelin {4.7, 6.7, 7.7, 8.7}; BASEL 9.0 per Maffei
    2021 Fig. 12 + Maffei 2025 Fig. 13 `>10⁹ pfu/ml if possible`) for semantic correctness; CH04 rerun under the new
    encoding is bit-identical at fold level (affine invariance confirmed). Full CH05 rerun under new encoding deferred
    to a follow-up ticket. Connects to `plm-rbp-redundant`: same Guelin-bank-dependent-phage-features failure mode, here
    on a genuine cross-panel split rather than cross-family within Guelin. Supports `panel-size-ceiling`: the fix is
    panel expansion or panel-independent phage features (dispatched as new CH06 with four candidate arms including
    OOD-aware inference, pairwise proteome similarity, Moriniere receptor-class probabilities, tail-protein-restricted
    TL17 projection), not richer engineered features. Canonical artifacts:
    lyzortx/generated_outputs/ch05_unified_kfold/ch05_combined_summary.json, ch05_{bacteria,phage}_axis_metrics.json,
    ch05_cross_source_breakdown.csv, ch05_{bacteria,phage}_axis_predictions.csv, ch05_per_family_breakdown.csv,
    ch05_straboviridae_exclusion.csv, ch05_reliability_tables.csv, ch05_basel_feature_variance.csv,
    ch05_basel_zero_vector_split.csv.*
- **`spandex-unified-kfold-baseline`**: HISTORICAL (SPANDEX-era). Superseded as the active unified-panel reference by
  chisel-unified-kfold-baseline (CH05). SX15 unified Guelin+BASEL k-fold baseline (2026-04-15, default BASEL+→MLC=2;
  2026-04-18 re-headline under CHISEL frame): bacteria-axis AUC 0.8685, phage-axis AUC 0.8988, cross-source near-parity
  AUC 0.896 (BASEL) vs 0.899 (Guelin). Phage-axis (all-pairs only; held-out phages unseen) was the first honest
  deployability estimate for unseen phages. The BASEL-vs-Guelin nDCG comparison (0.8332 vs 0.7193) was dropped — it was
  a metric artifact because BASEL binary labels have fewer nDCG rungs than Guelin's MLC grades, so direct nDCG
  comparison is not interpretable. AUC is comparable across sources and shows BASEL phages generalize essentially
  identically to Guelin phages. [validated; source: SX15; see also: chisel-unified-kfold-baseline,
  spandex-final-baseline, stratified-eval-framework, new-phage-generalization, external-data-neutral,
  per-phage-retired-under-chisel, cv-group-leakage-fixed]
  - *Panel: 369 bacteria (all 25 BASEL ECOR overlap with Guelin; no new bacteria) × 148 phages (96 Guelin + 52 BASEL) =
    33,202 observed pairs. The 3 pp bacteria-axis vs phage-axis AUC gap (0.8685 → 0.8988) looks like a "phage-axis is
    easier" effect but really reflects the deployment mix — phage-axis folds hold out entire phages so each test pair
    has more host-side signal available per training positive. Under CHISEL, this baseline is superseded by
    chisel-unified-kfold-baseline (established in CH05) after the cv_group leakage fix (CH02), concentration-feature
    flip (CH04), and per-phage retirement (per-phage-retired-under-chisel). SX15 numbers were computed with pair-level
    any_lysis training and per-phage blending enabled under the leaky cv_group folds; all three of those confounders are
    closed under CHISEL. Artifact paths:
    lyzortx/generated_outputs/sx15_eval/sx15_{bacteria,phage}_axis_stratified_metrics.csv and
    sx15_{bacteria,phage}_axis_predictions.csv.*
- **`autoresearch-baseline`**: AUTORESEARCH all-pairs model (0.810 AUC, 90.8% top-3 on ST03 holdout) is the canonical
  clean baseline: derived from raw FASTA, no leakage, no feature mismatch, no per-phage blending. Track GIANTS improved
  this to 0.823 AUC with depolymerase × capsule features and RFE. [validated; source: 2026-04-08 AUTORESEARCH eval,
  GT03; see also: tl18-flawed-baseline, defense-ablation-autoresearch, per-phage-blending-dominant,
  depo-capsule-validated]
  - *The 0.810→0.823 gain is entirely from Gate 1 (depolymerase × capsule) cross-terms + RFE feature selection. Neither
    HPO (+0.4pp, GT04), CatBoost (+0.3pp, GT05), nor expanded k-mer receptor predictions (+0.0pp, GT06) breaks through
    the 0.823 ceiling — it is feature-bound, not model-bound.*
- **`defense-ablation-autoresearch`**: Adding 79 DefenseFinder host defense features to AUTORESEARCH improves AUC from
  0.810 to 0.817 (+0.7pp) but regresses top-3 hit rate from 90.8% to 86.2% (-4.6pp). [validated; source: 2026-04-08
  defense ablation; see also: defense-lineage-confounding, autoresearch-baseline, ml-pipeline-over-features]
  - *Defense features help discrimination (AUC) but hurt ranking (top-3) due to lineage confounding. LightGBM's native
    feature selection is not aggressive enough to ignore noisy defense subtypes. GenoPHI found RFE is the optimal
    feature selection method.*
- **`per-phage-blending-dominant`**: HISTORICAL (SPANDEX-era, retired under CHISEL — see
  per-phage-retired-under-chisel). Per-phage LightGBM sub-models blended with all-pairs predictions were the dominant
  AUTORESEARCH architectural gain: +2.0pp AUC (0.810->0.830) on ST03 holdout, +3.1pp top-3 (90.8%->93.8%), and -2.3pp
  Brier (0.167->0.144). Surpasses TL18 on AUC (+0.7pp) and matches top-3 (93.8% vs 93.7%). [validated; source:
  2026-04-09 APEX ablation, 2026-04-09 APEX holdout; see also: per-phage-retired-under-chisel, per-phage-not-deployable,
  cv-group-leakage-fixed, chisel-baseline, deployment-goal]
  - *Each phage gets its own 32-tree LightGBM on host-only features (surface + typing + stats), blended 50/50 with
    all-pairs predictions. Bootstrap CIs overlap with TL18 — differences not statistically significant on 65-bacteria
    holdout. The +2 pp gain was inflated by two confounders that are both closed under CHISEL: (1) the pre-CH02
    fold-leakage bug (45 of 48 multi-bacterium cv_groups split across folds) let per-phage models memorize
    bacterium-level priors and recognize near-duplicate bacteria in held-out folds; (2) CH04's per-row diagnostic showed
    per-phage deflates positive predictions by ~5 pp under per-row training because the per-phage head sees 9 rows per
    bacterium with identical host features but mixed labels (no visibility to log_dilution). Not to be re-enabled in
    CHISEL or successor tracks.*
- **`adsorption-first-strategy`**: Adsorption-first modeling (host surface + typing features) is the correct critical
  path; defense features contribute but are not gate-critical for first baseline. [validated; source: 2026-04-05 replan,
  antiphage-landscape reading; see also: autoresearch-baseline]
- **`ml-pipeline-over-features`**: ML pipeline configuration (algorithm, training strategy, feature selection) matters
  5-18x more than genomic feature representation for strain-level phage-host prediction quality. RFE is the most
  impactful single pipeline change we've validated. [validated; source: Noonan 2025 GenoPHI, 2026-04-11 literature
  review, GT03, GT04, GT05; see also: defense-ablation-autoresearch, defense-lineage-confounding, genophi-benchmark,
  depo-capsule-validated]
  - *GenoPHI's 13.2M training run sweep: algorithm choice DMCC=0.365, training strategy DMCC=0.203, feature selection
    DMCC=0.184, but genomic representation only DMCC=0.020-0.072. Our GT03-GT05 confirm: RFE adds +1.2pp over baseline
    (validated), but HPO (+0.4pp, GT04) and CatBoost (+0.3pp, GT05) are not significant — the 0.823 ceiling is
    feature-bound, not model-bound. CatBoost does improve Brier (0.152 vs 0.161) without AUC penalty.*
- **`genophi-benchmark`**: GenoPHI achieves AUROC 0.869 on E. coli using 94/96 of our Guelin phages with binary
  protein-family and k-mer features, directionally outperforming our 0.830 with richer engineered features. [validated;
  source: Noonan 2025 GenoPHI, 2026-04-11 literature review; see also: ml-pipeline-over-features, autoresearch-baseline,
  tl18-flawed-baseline]
  - *Not directly comparable (different holdout strategies and bacteria counts). GenoPHI uses 402 strains; we use 369.
    Their features are annotation-free binary presence-absence from both phage and host proteomes. Our features are
    mechanistically grounded continuous scores. The approaches are complementary — combining ML pipeline optimizations
    with richer feature engineering may outperform either alone.*

## Evaluation & Benchmarking

Holdout protocol, benchmark methodology, and error analysis.

- **`st03-canonical-benchmark`**: ST03 grouped host split (65 cv_group-disjoint bacteria) is the canonical v1 benchmark;
  both TL18 and AUTORESEARCH are evaluated on it for honest comparison. [validated; source: TF01, 2026-04-08
  AUTORESEARCH eval; see also: split-contract]
- **`bootstrap-strain-level`**: Bootstrap CIs must be computed at holdout-strain level (not pair level) to align
  evaluation denominator with recommendation metric; 1000 resamples on 65 strains. [validated; source: TF01]
- **`error-buckets`**: The all-pairs model misses 6/65 holdout bacteria in top-3. Failure modes: Straboviridae prior
  collapse (2 strains rescued by per-phage, not deployable), abstention (2 strains with zero positives),
  needle-in-haystack (1 strain with 1/96 positive), and narrow-host prior collapse (1 strain with 14 narrow-host
  positives). [validated; source: ST09, TF02, 2026-04-09 APEX holdout; see also: family-bias-straboviridae,
  narrow-host-prior-collapse, per-phage-not-deployable]
  - *Per-phage blend reduced misses to 4/65, but per-phage is not deployable for unseen phages. The 6/65 all-pairs error
    rate is the operationally relevant baseline for the deployment goal.*
- **`family-bias-straboviridae`**: Straboviridae prior collapse suppresses cross-family true positives: broad-host
  Straboviridae (62-71% lysis rate) dominate model rankings, pushing narrow-host true positives (10-53% lysis rate)
  below the top-3 cutoff. [validated; source: ST09, 2026-04-09 APEX holdout NILS53 analysis]
- **`narrow-host-prior-collapse`**: Per-phage models help broad-to-moderate phages distinguish their hosts but cannot
  override the broad-phage prior for narrow-host phages (<30% lysis rate) with few training positives; NILS53 (14
  narrow-host TPs) remains a top-3 miss even with directed receptor features (GT06). [validated; source: 2026-04-09 APEX
  holdout NILS53 analysis, GT06, 2026-04-12 pair-level analysis; see also: family-bias-straboviridae,
  receptor-specificity-solved, per-phage-blending-dominant, omp-score-homogeneity]
  - *Receptor-class-directed features do not help because host OMP scores are near-constant (CV 0.01-0.17). NILS53's 13
    narrow-host TPs (Autographiviridae, Dhillonvirus, Kagunavirus; 10-26% lysis rate) are ranked 33-48, buried under
    broad-host phages. OMP allele-level features or depolymerase-capsule specificity for these narrow-host phages are
    the remaining paths.*
- **`inner-val-unreliable-top3`**: Inner-val top-3 predictions do not reliably transfer to holdout: per-phage blend
  predicted a top-3/AUC trade-off on inner-val (93.2% vs 94.6%) that did not replicate on holdout (93.8% vs 90.8% — both
  improved). [validated; source: 2026-04-09 APEX holdout; see also: st03-canonical-benchmark, bootstrap-strain-level]
  - *On 65-74 bacteria evaluation sets, 1 strain flip = 1.4-1.5pp top-3. Inner-val bacteria overlap with training
    distribution; holdout bacteria are cv_group-disjoint. Holdout is the only honest top-3 test.*
- **`ranking-metrics-retired`**: Top-3 hit rate, nDCG, and mAP are all retired from the scorecard starting with Track
  CHISEL (2026-04-18). Ranking is a product-layer concern — how many phages to administer, what selection policy to
  apply — not a property of the biological model. The biological model predicts P(lysis | bacterium, phage,
  concentration); downstream product code turns that calibrated probability into rankings or recommendations. The CHISEL
  scorecard is AUC + Brier only: AUC measures pairwise discrimination (does the model rank a true positive above a true
  negative?), Brier measures calibration (are the predicted probabilities honest?). Those are the two questions a
  pairwise lysis predictor can answer. Top-3 specifically was retired earlier in SPANDEX for separate reasons (collapses
  to binary signal, can't handle mixed-source panels) — CHISEL generalises that argument to all ranking metrics.
  [validated; source: 2026-04-12 SPANDEX design, 2026-04-18 CHISEL design; see also: st03-canonical-benchmark,
  bootstrap-strain-level, label-policy-binary]
- **`kfold-cv-replaces-fixed-holdout`**: 10-fold bacteria-stratified cross-validation replaces the single fixed ST03
  holdout (65 bacteria) as the primary evaluation protocol in Track SPANDEX. Every bacterium is evaluated exactly once;
  performance estimates are robust to holdout composition. [preliminary; source: 2026-04-12 SPANDEX design; see also:
  st03-canonical-benchmark, split-contract, ranking-metrics-retired]
  - *ST03 remains available as a single-fold comparison column for backwards compatibility with GIANTS results. The 10x
    model fitting cost is acceptable with LightGBM (minutes, not hours). k-fold also enables cross-source evaluation
    when BASEL bacteria overlap with different folds.*
- **`new-phage-generalization`**: The model generalizes moderately to unseen BASEL phages (AUC 0.72, nDCG 0.65) trained
  only on Guelin panel data — well above chance but below within-panel performance (AUC 0.87). Host-side features
  transfer; phage-specific features (zero-filled projection/PLM for BASEL) cap ranking quality. [validated; source:
  SX03; see also: deployment-goal, per-phage-not-deployable, basel-binary-only]
  - *BASEL phage features were zero-filled for phage_projection (no TL17 BLAST DB) and PLM PCA. Only phage_stats (GC,
    genome length) and depo cross-terms provided non-zero phage differentiation. Computing proper RBP family features
    for BASEL phages would likely close part of the AUC gap.*

## Deployment & Train/Inference Parity

Feature derivation parity, raw-input pipeline, and pre-computation.

- **`deployment-goal`**: The overarching deployment goal is a model that produces reliable lysis-likelihood inference on
  unseen E. coli strains, ranks or recommends cocktails from a set of potentially unseen phages, and generalizes along
  both the host and phage axes simultaneously. [validated; source: 2026-04-09 project direction; see also:
  per-phage-not-deployable, receptor-specificity-solved]
- **`per-phage-not-deployable`**: Per-phage sub-models (0.830 AUC) are not deployable and not a valid comparison
  baseline: they require training-time interaction data for each phage and cannot produce predictions for unseen phages.
  [validated; source: 2026-04-09 APEX holdout, 2026-04-09 project direction; see also: per-phage-blending-dominant,
  per-phage-retired-under-chisel, deployment-goal, autoresearch-baseline]
  - *The +2.0pp AUC gain over all-pairs is real on the fixed-panel holdout but non-transferable. All future track
    comparisons should use the all-pairs 0.810 baseline, not per-phage 0.830.*
- **`per-phage-retired-under-chisel`**: CHISEL retires per-phage blending (AX02) entirely. Do not re-enable it in CHISEL
  or successor tracks. Rationale: (1) `per-phage-not-deployable` — per-phage models memorize phage-specific
  bacterium-level priors and cannot transfer to unseen phages, which contradicts the CHISEL `deployment-goal` of
  generalization along both host and phage axes. (2) The SPANDEX +2 pp bacteria-axis gain was partly a leakage artifact
  — `cv-group-leakage-fixed` closed a fold-hashing bug that let per-phage models recognize near-duplicate bacteria
  across the train/test boundary; the post-fix per-phage gain was never cleanly measured. (3) CH04's per-row diagnostic
  found per-phage deflates positive predictions by ~5 pp under per-row training (the per-phage head sees 9 rows per
  bacterium with identical host features but mixed labels and no log_dilution visibility, collapsing to mean lysis rate
  across dilutions). (4) The phage-axis generalization question per-phage was silently answering is measured directly by
  CH05 (phage-axis k-fold) with held-out phages — a cleaner, deployment-aligned measurement. [validated; source: CH04,
  2026-04-19 CHISEL design; see also: per-phage-blending-dominant, per-phage-not-deployable, chisel-baseline,
  cv-group-leakage-fixed, deployment-goal]
- **`picard-assemblies`**: Picard Figshare assemblies (403 FASTAs, 1.9GB, CC BY 4.0) are the authoritative source for
  raw-input feature derivation; not all 403 have interaction labels. [validated; source: DEPLOY01]

## Dead Ends

Compressed lessons from approaches that didn't work.

- **`external-data-neutral`**: VHRdb, BASEL, KlebPhaCol, and GPB external interaction datasets showed neutral cumulative
  lift over the internal-only baseline; adding them did not improve predictions. [validated; source: TK01, TK02, TK03,
  TI09, SX03; see also: basel-binary-only, genophi-data-identical, new-phage-generalization]
  - *TK02 was invalidated (zero joined rows). SX03 is the first proper BASEL test with full feature computation
    (Pharokka + DepoScope on 52 genomes): Arm B (our data + BASEL training) gave nDCG +0.3pp with overlapping CIs vs
    baseline — confirmed neutral. 1,240 BASEL pairs (3.8% of training) is too small to move the needle.*
- **`ordinal-regression-not-better`**: Five loss formulations for ordinal potency prediction — vanilla regression
  (SX04), hurdle two-stage, LambdaRank, ordinal all-threshold (SX11), plus the SX11 binary control — are all null under
  AUC. SX11 ordinal all-threshold was -0.2 pp AUC vs binary baseline; LambdaRank was a -3.4 pp AUC regression. Every
  loss formulation tested so far fails to improve discrimination over a plain binary LightGBM on the panel. The
  within-positive reranking effect SX11 captures (Kendall tau 0.208→0.290 among positives) is a product-layer ranking
  concern that does not translate to discrimination or calibration gains. [validated; source: SX04, SX11; see also:
  mlc-dilution-potency, ranking-metrics-retired, label-policy-binary]
- **`label-derived-features-leaky`**: Label-derived features (legacy_label_breadth_count, defense_evasion_*,
  receptor_variant_ seen_in_training_positives) caused severe leakage and were removed entirely from TL18. [validated;
  source: TG04, TG05, TG06, TG08, TG12; see also: pairwise-block-leaky]
- **`pairwise-block-leaky`**: 5/13 pairwise compatibility features were soft-leaky (all 4 defense-evasion features + 1
  receptor training-positive flag); the 8 clean ones individually cannot recover the AUC gap without degrading top-3.
  [validated; source: TG08, TG11, TG12]
- **`mechanistic-features-no-lift`**: Mechanistic pairwise features (RBP-receptor, anti-defense pairs) showed no
  statistically significant lift on honest holdout rerun; kept in TL18 bundle but contribute only 3.5% and 2.0% to
  feature importance. [validated; source: TL12, TL18 audit; see also: adsorption-first-strategy]
- **`physicochemical-rbp-insufficient`**: Physicochemical RBP protein descriptors (AA composition +
  MW/GRAVY/pI/aromaticity/charge, 28 features from 80/96 phages) add no predictive signal; bulk protein properties
  cannot capture binding-interface specificity. [validated; source: 2026-04-09 APEX ablation, 2026-04-09 APEX holdout;
  see also: receptor-specificity-solved, adsorption-dominates-paper]
- **`phage-functional-noise`**: Phage functional gene repertoire features (PHROG category counts/fractions,
  anti-defense, depolymerase; 25 features) degrade top-3 from 94.6% to 87.8% on inner-val and are effectively noise for
  this prediction task. [validated; source: 2026-04-09 APEX ablation; see also: defense-lineage-confounding,
  defense-ablation-autoresearch]
- **`pairwise-cross-terms-dead-end`**: Both undirected (AX03: has_rbp × receptor_score) and directed (GT02/GT06:
  predicted_receptor_is_OmpC × host_OmpC_score) receptor × OMP cross-terms show no improvement on ST03 holdout — OMP HMM
  scores are too homogeneous to provide host-side discrimination. [validated; source: 2026-04-09 APEX ablation,
  2026-04-09 APEX holdout, GT03, GT06, 2026-04-12 pair-level analysis; see also: inner-val-unreliable-top3,
  receptor-specificity-solved, omp-score-homogeneity]
  - *AX03 failed because it paired "any RBP" with "any receptor" without binding specificity. GT02/GT06 directed
    cross-terms failed for a different reason: the host-side OMP scores have CV 0.01-0.17, so the cross-term collapses
    to a phage-only feature. The replacement approach must use OMP allele/variant-level features (extracellular loop
    regions) rather than whole-gene HMM scores.*
- **`hpo-catboost-not-bottleneck`**: Neither Optuna HPO (+0.4pp AUC, GT04) nor CatBoost algorithm switch (+0.3pp AUC,
  GT05) significantly improves over LightGBM defaults on the GT03 feature set — the 0.823 AUC ceiling is feature-bound,
  not model-bound. [validated; source: GT04, GT05; see also: ml-pipeline-over-features, depo-capsule-validated]
  - *GT04 Optuna found a more complex model (108 leaves, 450 trees) with delta CI [-0.010, +0.015]. GT05 CatBoost with
    native categoricals and 50-trial HPO gave delta CI [-0.011, +0.017]. CatBoost improves Brier (0.152 vs 0.161)
    without AUC penalty — useful for calibration but not discrimination. Confirms GenoPHI's finding that algorithm
    choice matters less than feature selection once the algorithm family is fixed.*
- **`genophi-features-redundant`**: GenoPHI's binary protein-family features (MMseqs2 clustering, 28,389 clusters) add
  no significant lift when combined with mechanistic features (AUC 0.826 vs 0.823, delta CI [-0.002, +0.008]). The
  protein-family clusters encode the same biological information as our HMM-based features at different granularity.
  [validated; source: GT08; see also: genophi-benchmark, autoresearch-baseline, panel-size-ceiling]
  - *RFE retained 244/700 protein-family features (216 host + 28 phage) capturing 13.1% importance. But host
    protein-family clusters correlate with HMM scores (both measure gene presence), and phage clusters correlate with
    phage_projection features. GenoPHI's AUROC advantage (0.869 vs 0.823) is not due to feature representation — it
    comes from different holdout strategy, bacteria counts, or ML pipeline interactions.*
- **`panel-size-ceiling`**: The 0.823 AUC ceiling is bound by the 96-phage panel and 65-bacteria holdout, not by feature
  engineering, algorithm choice, or feature representation. Seven independent attempts (GT04-GT08) all failed to break
  through, each adding 0.0-0.4pp (none significant). [validated; source: GT03, GT04, GT05, GT06, GT07, GT08, 2026-04-12
  final assessment; see also: autoresearch-baseline, error-buckets, narrow-host-prior-collapse, depo-capsule-validated,
  ambiguous-label-noise]
  - *The 6/65 holdout misses decompose into: 2 abstention (0 positives), 1 needle-in-haystack (1/96 positive), 3
    broad-phage-prior failures. The 3 rescuable misses require promoting narrow-host specialist phages (10-25% lysis
    rate) above broadly lytic phages (60-65% lysis rate) — a ranking challenge that genomic features alone cannot
    resolve because the host-specificity factors (expression regulation, phase variation, co-evolutionary dynamics)
    leave no detectable genomic signatures in presence-absence or HMM-score features. Expanding the phage panel or the
    holdout is the path forward, not more features.*
- **`ambiguous-label-noise`**: 3,462 pairs (10% of training, 6.5% of holdout) are labeled as negative but have
  uninterpretable ('n') raw scores — the actual experimental result was ambiguous. Excluding these from training
  improves top-3 from 89.2% to 92.3% (+3.1pp), the largest single ranking improvement in Track GIANTS. [validated;
  source: GT09, 2026-04-12 image review; see also: error-buckets, panel-size-ceiling, label-policy-binary]
  - *Raw plaque image review (Zenodo 10.5281/zenodo.10202713) confirms the 'n' scores are genuinely ambiguous — plates
    show faint signals and physical artifacts that make definitive calls impossible. For NILS53, the model's top-ranked
    phages (LF82_P8, 536_P9, LF82_P9) all have 'n' scores. For ECOR-69, DIJ07_P1 and DIJ07_P2 have 'n' scores. The model
    may be correct and the labels wrong. Label quality — not feature quality — is the binding constraint for ranking
    performance.*
- **`genophi-data-identical`**: GenoPHI's E. coli interaction matrix (402 × 94 phages) is 100% identical to ours
  (37,788/37,788 pairs match after name normalization). Our data is already in their framework. No new training pairs
  available from GenoPHI for existing phages. [validated; source: GT09; see also: genophi-benchmark,
  raw-interactions-authority]
- **`kmer-receptor-expansion-neutral`**: Moriniere 2026's 815 receptor-predictive 5-mers fail to lift performance
  regardless of encoding path: as intermediate-classifier features (GT06: AUC 0.824 vs 0.823, delta CI [-0.005, +0.005])
  or as direct phage-side features (SX12: AUC 0.8722 vs 0.8699, delta +0.23 pp, CIs overlap). The Moriniere 815-kmer
  approach is exhausted on this panel. [validated; source: GT06, SX12; see also: omp-score-homogeneity,
  pairwise-cross-terms-dead-end, receptor-specificity-solved, plm-rbp-redundant, narrow-host-prior-collapse,
  panel-size-ceiling]
  - *Two failure modes overlap. (1) The k-mers were selected to discriminate receptor class on K-12 derivatives
    (BW25113/BL21) which lack capsule/O-antigen — they predict what we already know (receptor identity), not what we
    need (strain-level capsule penetration). (2) The k-mers are information-redundant with phage_projection (TL17 BLAST)
    — both encode phage sequence similarity at different granularities. SX12 RFE keeps 95/815 k-mers (11.7%) but they
    contribute ~5% of total importance vs 22% for depo×capsule; zero k-mers appear in the top-20 features. Per-bacterium
    analysis: ~10 genuine wins (e.g., IAI78 +0.12 nDCG, ECOR-25 +0.17) exactly offset by ~10 genuine losses (e.g.,
    ECOR-19 -0.36, EDL933 -0.24). NILS53 (the canonical narrow-host case) gains +0.011 nDCG — k-mers do not break
    narrow-host prior collapse. Not a LightGBM shortcoming; a feature-redundancy + panel-size ceiling.*
- **`host-omp-variation-unpredictive`**: Host-side OMP allelic variation is substantial (369 clinical E. coli span BTUB
  28 MMseqs2 99%-identity clusters, OMPC 49, 5546 retained 5-mers across 12 core OMPs) but does not predict lysis. Four
  arms tested in SX13 (k-mers marginal, k-mers × phage k-mers cross-term, cluster IDs, plus baseline) all land within
  ±0.4 pp of SX10 on all metrics. [validated; source: SX13; see also: omp-score-homogeneity,
  pairwise-cross-terms-dead-end, kmer-receptor-expansion-neutral, narrow-host-prior-collapse,
  same-receptor-uncorrelated-hosts]
  - *Permutation test on cross_term aggregate delta: 73% of random prediction swaps are as extreme — signal
    indistinguishable from noise. Small directional signal in the marginal arm for moderate-narrow-host deciles (lysis
    5-11%, mean +1-2 pp nDCG across ~70 bacteria), biologically consistent with OMP-variation-matters-for-narrow-phages
    but not significant by sign test (75/137 positive, p=0.31). NILS53 specifically improved +2.59 pp under cross_term
    but sits at the 76th percentile of its lysis-rate peers (mean Δ = -0.001) — one outlier draw, not a reproducible
    rescue. Top-3 hit rate moves from 92.2% to 93.3% under cross_term (+4 strains). The finding refines
    omp-score-homogeneity: HMM-score homogeneity was real, but escalating to finer representations doesn't rescue
    prediction — host-range variance lives downstream of OMP recognition (polysaccharide access, intracellular defenses,
    co-evolutionary dynamics).*
- **`label-vision-reading-spot-checked-dead`**: Using a vision model to re-read the ambiguous 'n' plaque-image scores
  (the plate crops backing the ~10% of training rows labeled negative but with uninterpretable raw scores) was evaluated
  via manual spot checks before 2026-04 and did not look promising enough to justify a full re-read pipeline. Do not
  propose vision-based 'n' resolution as a feature / label fix; the label noise floor is accepted as-is for now.
  [validated; source: 2026-04 spot checks before SPANDEX wave-2 plan; see also: ambiguous-label-noise,
  panel-size-ceiling]
  - *The underlying ambiguity persists — ~3,462 pairs have score='n' with no interpretable positives and are labeled 0.
    GT09 showed excluding these pairs alone improves top-3 by +3.1 pp (ambiguous-label-noise). But a vision-driven label
    correction loop isn't a viable shortcut given the plate image quality. Future tracks should either expand the
    wet-lab panel or accept the label noise as a fixed ceiling.*
- **`plm-rbp-redundant`**: ProstT5+SaProt PLM embeddings of RBP sequences (1280-dim, PCA to 32) achieve 33.9% LightGBM
  feature importance but zero predictive lift on ST03 holdout; they cannibalize existing phage family features without
  adding new discriminative signal. [validated; source: 2026-04-10 AX08 holdout, 2026-04-16 APEX retro stratified audit;
  see also: physicochemical-rbp-insufficient, receptor-specificity-solved, kmer-receptor-expansion-neutral,
  phage-functional-noise]
  - *PLM embeddings encode protein-level similarity that correlates with genome-level family membership. Mean-pooling
    across RBPs further dilutes binding-specific signal. Moriniere 2026 confirms the signal is in localized
    hypervariable sequence motifs at RBP tip domains, not global protein embeddings — k=5 amino acid k-mers predict
    receptor class at AUROC 0.99. Retroactive stratified audit (2026-04-16) confirmed null across all SX14 strata: PLM
    features give +0.89 pp within-family nDCG (CIs overlap baseline) but -1.22 pp cross-family nDCG — phage features
    that encode family structure actively hurt cold-start predictions because they inject within-family priors the model
    misapplies to unseen families. Same cross-family damage pattern observed for phage_functional (AX04, -0.40 pp
    cross-family) and the combined AX06 config (-0.86 pp cross-family). Unlike SX11's loss- function changes (+3.5 pp
    within-family with disjoint CI), feature-level APEX additions are null at every stratum level. AUC flat at 0.849
    across all arms — no discrimination signal.*

## Open Questions

Unresolved items that still matter for the project direction.

- **`narrow-susceptibility-features`**: Can host-compatibility features (receptor/defense) lift narrow-susceptibility
  recovery for the 12/36 resolved narrow strains not rescued by any panel phage? [preliminary; source: ST09, TB04; see
  also: narrow-susceptibility-rescue, error-buckets]
- **`within-family-reranking`**: Can within-family reranking improve phage selection inside saturated score bands?
  Receptor-class features cannot help (GT06: same-receptor phages have uncorrelated host ranges, Jaccard 0.091 for Tsx).
  OMP allele-level features or GenoPHI binary protein-family features are the remaining paths. [preliminary; source:
  ST09, GT06, 2026-04-12 pair-level analysis; see also: error-buckets, family-bias-straboviridae,
  receptor-specificity-solved, omp-score-homogeneity, same-receptor-uncorrelated-hosts]
  - *Receptor-class-directed features failed because: (1) host OMP scores are near-constant, and (2) same-receptor
    phages have wildly different host ranges. GenoPHI's binary protein-family features from both host and phage
    proteomes (MMseqs2 clustering) are the most validated alternative — they achieve AUROC 0.869 with annotation-free
    features.*
- **`defense-top3-mitigation`**: RFE mitigates the defense top-3 regression: GT03 all_gates_rfe achieves 90.8% top-3 (vs
  baseline 92.3% and defense-only 89.2%), while capturing the AUC gain (0.823 vs 0.810). The regression is reduced but
  not eliminated. [validated; source: 2026-04-08 defense ablation, GT03; see also: defense-ablation-autoresearch,
  defense-lineage-confounding, ml-pipeline-over-features, depo-capsule-validated]
  - *RFE prunes ~250 of 502 features, removing redundant defense subtypes that cause lineage confounding. The
    all_gates_rfe arm shows defense at only 4.6% feature importance (vs 21% for depo×capsule), suggesting RFE correctly
    deprioritizes confounded features.*
- **`receptor-directed-features`**: Directed receptor × OMP cross-terms (predicted_receptor=OmpC × host_OmpC_score) do
  not break the feature plateau — tested with both genus-level (GT02, 8/96 phages) and k-mer (GT06, 39/96 phages)
  receptor predictions, both show zero lift. The host-side OMP HMM scores are too homogeneous. [validated; source: GT02,
  GT03, GT06, 2026-04-12 pair-level analysis; see also: receptor-specificity-solved, pairwise-cross-terms-dead-end,
  narrow-host-prior-collapse, omp-score-homogeneity]
  - *Answered: No. The cross-term collapses because host OMP scores have CV 0.01-0.17. The replacement path is OMP
    allele-level features (extracellular loop variants) that provide actual host-side variation. OmpC has 50 allelic
    variants, BtuB has 28 — the information exists but is compressed away by whole-gene HMM scoring.*
