# Project Knowledge Model

<!-- Last consolidated: 2026-04-13T01:40:00+02:00 -->
<!-- Source: lyzortx/research_notes/lab_notebooks -->

**58 knowledge units** across 7 themes (44 active, 14 dead ends)

## Data & Labels

Labeling policy, data quality, split contracts, and training corpus.

- **`label-policy-binary`**: Binary lysis labels use any_lysis rule (matrix_score > 0); ~7% of positives are single-hit
  noise with matrix_score=0 despite positive labels. [validated; source: ST0.1]
- **`split-contract`**: ST02 defines 294 training bacteria; ST03 reserves 65 cv_group-disjoint bacteria for sealed
  holdout. Deterministic cv_group hashing with salt ensures reproducibility. [validated; source: ST0.2, ST0.3, AR01; see
  also: raw-interactions-authority]
- **`raw-interactions-authority`**: raw_interactions.csv is the authoritative training corpus; all derived training
  cohorts and evaluation splits trace back to this file. [validated; source: ST0.2, AR01, DEPLOY01]
- **`mlc-dilution-potency`**: MLC (minimum lytic concentration) is a derived score, not a raw feature —
  raw_interactions.csv contains only binary scores (0/1/n) at four tested log_dilutions (0, -1, -2, -4). The paper
  defines MLC over three replicated concentrations (5x10^8, 5x10^7, 5x10^6 pfu/ml) and explicitly excludes the
  unreplicated 5x10^4 (log_dilution=-4) data. Paper's MLC=3 vs MLC=4 is a morphological distinction at 5x10^6
  (individual plaques vs entire lawn lysis) that our binary spot data cannot reproduce. SX05 applied the fix:
  DILUTION_WEIGHT_MAP = {0: 1, -1: 2, -2: 3} with EXCLUDED_LOG_DILUTIONS = {-4}; interaction_matrix.csv was regenerated
  by capping MLC=4 cells to MLC=3 (1288 pairs collapsed, MLC=0/1/2 counts preserved). Our MLC range is now {0, 1, 2, 3}.
  MLC=1 is a standard low-potency interaction — not a problem despite the paper's LFW/Abi caveat on lawn clearing at
  high phage concentration. [validated; source: Gaborieau 2024 Methods "Evaluating phage-bacteria interaction outcomes
  by plaque assay experiments", Gaborieau 2024 Fig 2b legend, 2026-04-13 raw CSV verification; see also:
  label-policy-binary, raw-interactions-authority, ambiguous-label-noise, top3-metric-retired]
  - *EXACT PAPER QUOTES (Methods, "Evaluating phage-bacteria interaction outcomes by plaque assay experiments" section):
    - MLC definition: "We encoded each phage-bacteria interaction using the MLC   score, which corresponds to the lowest
    concentration of the phage at which   a lytic interaction is observed." - What counts as lytic: "We considered an
    interaction to be lytic when we   observed individual lysis plaques or full clearing at any phage   concentration.
    Individual lysis plaques ascertain the lysis of the   bacterial population with production of new virions. Clearing
    of the   bacterial lawn at high phage concentration could result from productive   lysis of the bacterial population
    by the phage, or from another mechanism   such as lysis from without, or abortive infection." - Guelin MLC range:
    "The MLC score is null in the case of non-lytic   interaction and ranges from 1 (lytic interaction at the highest
    phage   titre) to 4 (uncountable number of lysis plaques at 5 x 10^6 pfu/ml)." - Replicate structure: "5 x 10^8
    pfu/ml (replicates R1, R2 and R3),   5 x 10^7 pfu/ml (replicates R2 and R3), 5 x 10^6 pfu/ml (replicates R1, R2
    and R3) and 5 x 10^4 pfu/ml (replicate R1)." - Exclusion of 5x10^4: "The outcome of interaction at 5 x 10^4 pfu/ml
    was   not taken into account in the calculation of the MLC score because it was   not verified by a replicate."
    EXACT PAPER QUOTE (Fig 2b legend): "0, no lytic interaction was observed; 1, lytic interaction at the highest phage
    titre (5 x 10^8 pfu/ml); 2, lytic interaction at the middle phage concentration (5 x 10^7 pfu/ml); interactions
    observed at the lowest phage concentration (5 x 10^6 pfu/ml) are distinguished between 3 (individualized lysis
    plaque) and 4 (entire lysis of the bacterial lawn)." RAW DATA SHAPE (raw_interactions.csv, verified 2026-04-13):
    columns are bacteria, phage, image, replicate, plate, log_dilution, X, Y, score. Only binary scores {0, 1, n} and
    only four log_dilutions {0, -1, -2, -4} — no -3 (5x10^5 was never tested). No plaque-morphology column, so the
    paper's MLC=3 vs MLC=4 morphological split cannot be reconstructed. SX05 IMPLEMENTATION (2026-04-13):
    DILUTION_WEIGHT_MAP reduced to {0: 1, -1: 2, -2: 3}; EXCLUDED_LOG_DILUTIONS = {-4} added to
    build_track_a_foundation.py and applied at raw-row ingestion in both track_a and build_contract (autoresearch). The
    paper's matrix values at MLC=4 cannot be reconstructed from binary raw data, so regenerate_interaction_matrix.py
    caps MLC=4 cells at MLC=3 in data/interactions/interaction_matrix.csv (1288 pairs collapsed; MLC=0/1/2 counts
    preserved exactly; total pair count 38592 unchanged). MLC=1 stays as-is (NOT suspect — standard low-potency
    interaction; the paper's LFW/Abi caveat is a biological note, not an exclusion criterion). USE: MLC 0-3 provides
    graded relevance weights for nDCG evaluation (SX01). Training labels remain binary (any_lysis) per
    label-policy-binary.*
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
- **`omp-score-homogeneity`**: All 369 clinical E. coli hosts express all 12 core OMP receptors at nearly identical HMM
  scores (CV 0.01-0.17), making receptor × OMP cross-terms collapse to phage-only features with no host-side
  discrimination. [validated; source: GT03, GT06, 2026-04-12 pair-level analysis; see also: receptor-variant-richness,
  receptor-specificity-solved, three-layer-hypothesis]
  - *LptD CV=0.01, OmpA CV=0.02, Tsx CV=0.03, OmpC CV=0.04. The cross-term predicted_is_OmpC × host_OmpC_score ≈
    predicted_is_OmpC × constant. Direct test: OmpC phages lyse high-OmpC hosts at 37.0% vs low-OmpC at 33.4% (Cohen
    d=0.086, negligible). Moriniere's receptor classifiers were trained on K-12 (BW25113/BL21) which lack
    O-antigen/capsule — in clinical strains, polysaccharide barriers mask OMP differences. OMP variant/allele-level
    features (extracellular loop regions) might restore host-side variation.*
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
- **`spandex-final-baseline`**: SPANDEX final baseline (Track SPANDEX closing configuration, 2026-04-13): GT03
  all_gates_rfe + AX02 per-phage blending, trained on SX05-corrected MLC 0-3 labels, with SX06 real TL17 family
  projection for BASEL phages. 10-fold CV on our 369×96 panel: nDCG 0.7958 [0.7877, 0.8124], mAP 0.7111 [0.6925,
  0.7290], AUC 0.8699 [0.8570, 0.8819], Brier 0.1248 [0.1187, 0.1309]. Cross-panel Arm C (train Guelin → predict BASEL ×
  ECOR): nDCG 0.7619 [0.7219, 0.8207], mAP 0.5186 [0.4591, 0.5780], AUC 0.7607 [0.6886, 0.8307], Brier 0.1844 [0.1426,
  0.2213]. SX07 and SX09 skipped (plm-rbp-redundant); SX08 continuous depolymerase bitscore validated as null
  (bit-identical Arm C metrics). [validated; source: SX05, SX06, SX08, SX10; see also: autoresearch-baseline,
  mlc-dilution-potency, new-phage-generalization, plm-rbp-redundant, panel-size-ceiling]
  - *The 10.9 pp within-panel AUC (0.87) vs cross-panel AUC (0.76) gap is the main SPANDEX-era unresolved item. Closing
    it requires panel expansion rather than richer phage-side features. Use this record as the reference point for
    future tracks; the canonical artifacts live at lyzortx/generated_outputs/sx05_sx01_eval/ (within-panel) and
    lyzortx/generated_outputs/sx06_sx03_eval/ (cross-panel Arm C).*
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
- **`per-phage-blending-dominant`**: Per-phage LightGBM sub-models blended with all-pairs predictions are the dominant
  AUTORESEARCH architectural gain: +2.0pp AUC (0.810->0.830) on ST03 holdout, +3.1pp top-3 (90.8%->93.8%), and -2.3pp
  Brier (0.167->0.144). Surpasses TL18 on AUC (+0.7pp) and matches top-3 (93.8% vs 93.7%). [validated; source:
  2026-04-09 APEX ablation, 2026-04-09 APEX holdout; see also: adsorption-dominates-paper, family-bias-straboviridae,
  autoresearch-baseline, per-phage-not-deployable, deployment-goal]
  - *Each phage gets its own 32-tree LightGBM on host-only features (surface + typing + stats), blended 50/50 with
    all-pairs predictions. Bootstrap CIs overlap with TL18 — differences not statistically significant on 65-bacteria
    holdout.*
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
- **`top3-metric-retired`**: Top-3 hit rate is retired as an evaluation metric starting with Track SPANDEX. It collapses
  the entire ranking into a binary signal, discards potency information, and cannot handle sparse ground truth from
  multi-source panels. Replaced by nDCG (graded relevance using MLC 0-3 post-SX05) and mAP (binary retrieval quality).
  [validated; source: 2026-04-12 SPANDEX design; see also: mlc-dilution-potency, st03-canonical-benchmark,
  bootstrap-strain-level]
  - *Top-3 assigns identical scores to a model ranking a true positive 4th vs 148th. With mixed-source data (our MLC 0-3
    post-SX05 + BASEL binary), nDCG naturally handles graded relevance while mAP handles partial ground truth (score
    only observed pairs per bacterium). AUC and Brier retained as secondary metrics.*
- **`kfold-cv-replaces-fixed-holdout`**: 10-fold bacteria-stratified cross-validation replaces the single fixed ST03
  holdout (65 bacteria) as the primary evaluation protocol in Track SPANDEX. Every bacterium is evaluated exactly once;
  performance estimates are robust to holdout composition. [preliminary; source: 2026-04-12 SPANDEX design; see also:
  st03-canonical-benchmark, split-contract, top3-metric-retired]
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
  deployment-goal, autoresearch-baseline]
  - *The +2.0pp AUC gain over all-pairs is real on the fixed-panel holdout but non-transferable. All future track
    comparisons should use the all-pairs 0.810 baseline, not per-phage 0.830.*
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
- **`ordinal-regression-not-better`**: Five loss formulations for ordinal potency prediction (vanilla regression/SX04,
  hurdle two-stage, LambdaRank, ordinal all-threshold/SX11) all fail the +2 pp nDCG gate over binary classification. The
  binary model's implicit potency signal (Spearman 0.246 among positives) captures most of the available signal. Best
  alternative is ordinal all-threshold at +1.33 pp nDCG — detectable but sub-threshold, with a trade-off: MLC=1 pairs
  get demoted. [validated; source: SX04, SX11; see also: mlc-dilution-potency, top3-metric-retired]
  - *79% zero-inflation (MLC=0) dominates all metrics. Binary P(lysis) implicitly ranks potency because MLC=3 pairs have
    3 concordant positive training rows vs MLC=1 pairs with 1 positive and 2 negatives. SX11 ordinal all-threshold shows
    a real within-positive reranking effect (Kendall tau 0.208→0.290, 62% of bacteria improve) but equal-weighted
    threshold combination suppresses MLC=1 pairs relative to binary — 55% of MLC=1 pairs drop in rank. The mechanism is
    sound; the data's potency resolution (3 dilution levels) cannot support enough signal to clear the 2 pp gate. Richer
    potency labels (quantitative EOP) could revisit.*
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
- **`kmer-receptor-expansion-neutral`**: Expanding Gate 2 receptor coverage from 8/96 (genus-level) to 39/96
  (k-mer-based) OMP phages produces zero AUC improvement (0.824 vs 0.823, delta CI [-0.005, +0.005]). [validated;
  source: GT06; see also: omp-score-homogeneity, pairwise-cross-terms-dead-end, receptor-specificity-solved]
  - *Used GenoPHI's 815 receptor-predictive k-mers from Moriniere 2026 Dataset S6. The simplified k-mer max-vote
    predictor assigns 39 OMP, 33 LPS, 22 NGR, 2 unknown. Despite 5x more OMP-assigned phages, the cross-terms still
    collapse because the host-side OMP scores have CV 0.01-0.17. More accurate receptor predictions (full GenoPHI
    classifier) would not help because the host side — not the phage side — is the bottleneck.*
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
  adding new discriminative signal. [validated; source: 2026-04-10 AX08 holdout; see also:
  physicochemical-rbp-insufficient, receptor-specificity-solved]
  - *PLM embeddings encode protein-level similarity that correlates with genome-level family membership. Mean-pooling
    across RBPs further dilutes binding-specific signal. Moriniere 2026 confirms the signal is in localized
    hypervariable sequence motifs at RBP tip domains, not global protein embeddings — k=5 amino acid k-mers predict
    receptor class at AUROC 0.99.*

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
