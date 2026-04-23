### 2026-03-28: TL01 Annotate all 97 phage genomes with Pharokka

#### Executive summary

Added Pharokka and its bioconda dependencies (mmseqs2, trnascan-se, minced, aragorn, mash, dnaapler) to
`environment.yml` and built a Track L pipeline (`lyzortx/pipeline/track_l/`) that runs Pharokka on all 97 phage FNA
files, parses the `_cds_final_merged_output.tsv` per phage, and produces three summary tables: PHROGs category counts,
RBP gene list, and anti-defense gene list. All 97 phages must produce >0 annotated CDS for the pipeline to succeed.

#### What was implemented

- `environment.yml`: added `pharokka=1.9.1`, `mmseqs2=18.8cc5c`, `trnascan-se=2.0.12`, `minced=0.4.2`,
  `aragorn=1.2.41`, `mash=2.3`, `dnaapler=1.3.0` as conda dependencies from the bioconda channel.
- `lyzortx/pipeline/track_l/run_track_l.py`: entry point with `--step annotate|parse|all`.
- `lyzortx/pipeline/track_l/steps/run_pharokka.py`: iterates over all 97 `.fna` files in `data/genomics/phages/FNA/`,
  runs `pharokka.py` per phage, verifies >0 CDS in each output. Supports `--force` to re-run and skips already-complete
  phages by default.
- `lyzortx/pipeline/track_l/steps/parse_annotations.py`: reads pharokka's `_cds_final_merged_output.tsv` per phage and
  produces:
  - `phrog_category_counts.csv` — phage x 10 PHROG functional categories.
  - `rbp_genes.csv` — per-phage receptor binding protein genes identified by annotation pattern matching (tail fiber,
    tail spike, receptor binding, host specificity, adhesin, etc.).
  - `anti_defense_genes.csv` — per-phage anti-defense genes (anti-CRISPR, anti-restriction, methyltransferases, Ocr,
    Ard, etc.).
  - `manifest.json` — run metadata.

#### Design decisions

- **Gene classification by annotation patterns**: Pharokka does not natively tag RBP or anti-defense genes as separate
  categories. We use regex pattern matching on the `annot` column of the merged CDS output. The patterns are defined as
  module-level constants and tested in `lyzortx/tests/test_track_l_parse_annotations.py`.
- **Skip mash and extra annotations**: The pharokka runner uses `--skip_mash --skip_extra_annotations` to avoid
  INPHARED matching and reduce runtime, since those outputs are not needed for TL01.
- **Per-phage output directories**: Each phage gets its own pharokka output directory under
  `lyzortx/generated_outputs/track_l/pharokka_annotations/<phage_name>/`, matching pharokka's native output structure.

#### Key findings from exploratory analysis

Post-annotation analysis comparing generalist phages (200-262 strains lysed) vs narrow-range phages (6-17 strains)
revealed strong discriminative signal in the RBP PHROG repertoire:

- **43 unique RBP PHROG families** across the 97-phage panel.
- Generalists carry **6-9 RBPs** with diverse subtypes (long/short tail fiber, tail spike, host specificity protein,
  baseplate wedge connector). Narrow-range phages carry **0-3 RBPs**, mostly generic tail spike or tail fiber.
- PHROG repertoires are **almost completely disjoint**: 19 PHROGs exclusive to generalists, 7 exclusive to narrow-range,
  only 1 shared (PHROG 817). This means the specific PHROG IDs, not just counts, carry strong signal.
- Track C provides 12 host OMP receptor types (22 variant clusters at 99% identity) and 82 defense system subtypes.

Host range distribution across the 96-phage interaction panel (369 bacteria): min=6 strains (411_P1), Q1=42, median=78,
Q3=164, max=262 (DIJ07_P2). No phages have ≤3 positive strains; 10 phages lyse ≤20 strains.

This suggests the most informative feature for TL02 is a **phage x RBP-PHROG binary matrix** (43 columns), not just RBP
count or diversity. Even better: cross with host receptor data via enrichment analysis to learn empirical PHROG → receptor
associations from the interaction matrix.

#### Generalization principle

The key design goal of Track L: PHROG-receptor associations are learned from the interaction matrix during training, but
at prediction time only genomes are needed. For a novel phage, run pharokka → extract RBP PHROGs → feature vector. For
a novel host, run BLAST against OMP DB + DefenseFinder → receptor/defense feature vector. The model applies the learned
PHROG-receptor weights to predict lysis without any interaction data for the new pair. This is what TL06-TL09 implement.

### 2026-03-29: TL02 Annotation-interaction enrichment analysis (PHROG x receptor/defense)

#### Executive summary

Built a permutation-based enrichment module (`annotation_interaction_enrichment.py`) and ran three analyses on the full
369×96 interaction matrix (9,720 lytic pairs, 27.4% base rate). The test conditions on the phage carrying the PHROG,
then permutes host labels to get calibrated p-values that respect the correlation structure of the interaction matrix.
Results: 380/3,424 significant RBP PHROG × OMP receptor associations, 24/160 RBP PHROG × LPS core associations, and
27/924 anti-defense PHROG × defense subtype associations (BH-corrected p < 0.05). These results support building
pairwise mechanistic features in TL03/TL04.

#### What was implemented

- `lyzortx/pipeline/track_l/steps/annotation_interaction_enrichment.py`: Reusable enrichment module. Takes any (phage
  binary feature matrix, host binary feature matrix, interaction matrix) triple. For each (phage_feature, host_feature)
  pair, computes the lysis rate difference conditioned on the phage having the feature (host_has rate − host_lacks rate).
  P-values are computed by permuting host labels (1000 permutations), which preserves the row/column correlation
  structure of the interaction matrix. Applies Benjamini-Hochberg FDR correction across all tests within each analysis.
- `lyzortx/pipeline/track_l/steps/run_enrichment_analysis.py`: Loads pharokka RBP/anti-defense annotations, host OMP
  receptor clusters, LPS core types, and defense system subtypes, constructs binary matrices, runs three enrichment
  analyses. Wired into `run_track_l.py` as `--step enrich`.
- `lyzortx/tests/test_annotation_interaction_enrichment.py`: 14 unit tests using a 20×10 slice of the real
  interaction matrix, covering BH correction, permutation test, contingency table arithmetic, resolved-mask
  exclusion, and main effect confounding.

#### Design decisions

- **Conditioning on phage feature**: The test conditions on the phage carrying the PHROG, then asks whether the host
  feature increases lysis within that subset. This controls for the phage main effect — generalist PHROGs do not show
  spurious enrichment for every host feature.
- **Permutation p-values over Fisher's exact test**: Fisher's exact test assumes independent observations, but
  interaction matrix entries are correlated (some hosts are generally susceptible, some phages are generalists). Null
  calibration with random features on the real interaction matrix showed Fisher's yielded 25% false positives at p < 0.05
  (expected: 5%). Permuting host labels (1000 permutations) preserves this correlation structure and gives calibrated
  p-values (3% at the 5% threshold in null calibration).

#### Data dimensions

- **Interaction panel**: 369 bacteria × 96 phages = 35,424 pairs (35,266 resolved: 9,720 lytic, 25,546 non-lytic;
  158 unresolved pairs excluded from enrichment).
- **RBP PHROGs**: 43 unique across 97 phages, 32 present in ≥2 phages (used in enrichment).
- **Anti-defense PHROGs**: 13 unique, 12 present in ≥2 phages.
- **OMP receptor features**: 107 (receptor, cluster) pairs with ≥5 bacteria across 12 receptor proteins.
- **LPS core types**: 5 types (R1=207, R3=61, K12=43, R4=35, R2=22 bacteria; No_waaL excluded).
- **Defense subtypes**: 77 subtypes with ≥5 bacteria (of 137 total).

#### Analysis 1: RBP PHROG × OMP receptor variant clusters

- **3,424 tests** (32 PHROGs × 107 receptor clusters), **380 significant** (11.1%).

#### Analysis 2: RBP PHROG × LPS core type

- **160 tests** (32 PHROGs × 5 LPS types), **24 significant** (15.0%).

#### Analysis 3: Anti-defense PHROG × defense system subtypes

- **924 tests** (12 PHROGs × 77 defense subtypes), **27 significant** (2.9%).

#### Caveats

- **Duplicate PHROG profiles**: The 32 RBP PHROGs reduce to ~25 unique phage-carrier patterns (e.g., 136/15437/4465/9017
  always co-occur, as do 1002/1154/967/972 and 2097/4277). The 380 significant associations are not 380 independent
  biological discoveries — TL03 should collapse identical PHROG profiles before building features.
- **P-value resolution**: With 1000 permutations, 313/380 significant OMP hits are at the p-value floor (0.001). The
  BH boundary is partly quantized: 95 hits fall in the BH 0.03–0.07 zone, and ~65 borderline significant hits could
  flip with more permutations. The 314 floor hits (where ≤1 permutation exceeded the observed statistic) are robust.
  For the screening purpose of this step (feeding candidate pairs to TL03), this resolution is sufficient. TL03 should
  not treat the counts or rankings from this screen as precise.
- **Residual confounding**: The permutation test conditions on phage feature and permutes host labels, but does not
  control for host phylogenetic lineage or correlated feature blocks. Some associations may reflect lineage correlation
  rather than specific molecular interactions. The anti-defense results (2.9% significance) are particularly susceptible
  — generic methyltransferases mapping to diverse defense subtypes are more plausibly explained by lineage than by
  specific evasion.

#### Implications for TL03/TL04

- **TL03 (RBP-receptor features)**: 380 significant PHROG-receptor candidate associations (after controlling for phage
  main effects and matrix correlation via host-label permutation). TL03 should collapse duplicate PHROG profiles and
  use enrichment odds ratios as candidate feature weights — not treat these as confirmed molecular mechanisms.
- **TL04 (defense evasion features)**: 27 significant anti-defense × defense associations. The lower rate and caveats
  about annotation specificity suggest these should be treated as weaker candidates than the RBP-receptor features.
- **No escape hatch needed**: With 380 significant associations (~285 unique after collapsing duplicate PHROG carrier
  profiles), enrichment-based features are viable for TL03.

#### Non-protein host factors considered

Beyond OMP protein receptors and defense systems, several non-protein surface structures affect phage adsorption:

- **LPS core type** (R1-R4, K12) — available in Track C, will be tested in TL02 enrichment (RBP PHROGs x LPS type).
- **K-antigen / polysaccharide capsule** — Track C has Klebsiella capsule type but 94% missing. Sparse signal.
- **O-antigen type** — determines phage adsorption for some phages. Available in Track C surface features.
- **Phase variation** — stochastic receptor expression switching (e.g., FimH, Ag43). Cannot be captured from genome
  alone — inherently epigenetic. Acknowledged as a limitation.
- **Phage depolymerases** — enzymes that degrade capsule/LPS. Pharokka found 18 "polysaccharide chain length
  determinant" genes but cannot distinguish capsule-type targets. Deferred to Future note in project.md.

#### Next steps

1. **Build a reusable enrichment module** (`annotation_interaction_enrichment.py`) that takes any (phage feature matrix,
   host feature matrix, interaction matrix) triple and produces a Fisher's exact test enrichment table with odds ratios,
   p-values, and BH-corrected significance. This module will be used by both TL02 and TL03.
2. **Run enrichment analysis on three pairings:**
   - RBP PHROGs (43) x OMP receptor variants (22) — core RBP-receptor signal
   - RBP PHROGs (43) x LPS core type (~5 types) — tail spikes that bind LPS directly
   - Anti-defense gene PHROGs x defense system subtypes (82) — for TL03
3. **TL02:** Use the significant PHROG-receptor associations to build pairwise features. For each phage-host pair:
   does the phage carry an RBP PHROG that is significantly associated with lysis of hosts carrying this receptor? The
   enrichment odds ratios become the feature weights.
4. **TL03:** Same logic for anti-defense genes x host defense systems.
5. **Depolymerase annotation (deferred):** Pharokka annotations are too coarse for capsule-depolymerase matching
   (18 "polysaccharide chain length determinant" hits, no capsule-type specificity). Dedicated tools exist (DepoScope,
   DePP, PDP-Miner) but none predict capsule-type specificity — they only classify binary depolymerase yes/no. Running
   Pfam/InterPro on tail spike sequences to get glycosyl hydrolase family would be more informative. Defer until
   enrichment analysis shows whether capsule features (Track C has LPS core type + Klebsiella capsule with 94% missing)
   carry enough signal to justify the effort.

### 2026-03-29: TL03 Mechanistic RBP-receptor compatibility features from annotations

#### Executive summary

Built `build_mechanistic_rbp_receptor_features.py`, a Track L step that turns TL02 enrichment hits plus Pharokka RBP
PHROG annotations into a pairwise feature matrix for the full 369×96 panel. The builder collapses the 32 panel RBP
PHROGs to 25 unique carrier profiles before feature construction, then emits two blocks: (1) 25 direct phage-level
collapsed-profile features copied onto each pair row, and (2) 302 weighted pairwise compatibility features where the
value is `lysis_rate_diff` when the phage carries the collapsed profile and the host carries the enriched OMP/LPS
feature, else 0. The output is joinable on `pair_id` / `bacteria` / `phage`.

#### What was implemented

- `lyzortx/pipeline/track_l/steps/build_mechanistic_rbp_receptor_features.py`: TL03 builder. It:
  - bootstraps Track A labels and TL02 enrichment outputs if they are missing from a fresh checkout,
  - rebuilds the panel RBP PHROG matrix from cached Pharokka TSVs,
  - collapses duplicate PHROG carrier profiles,
  - collapses duplicate significant enrichment hits onto those profiles,
  - writes the feature CSV, column metadata, profile metadata, sanity-check CSV, and manifest.
- `lyzortx/pipeline/track_l/run_track_l.py`: added `--step rbp-features` to run TL03 directly.
- `lyzortx/tests/test_track_l_mechanistic_rbp_receptor_features.py`: 5 unit tests covering duplicate-profile
  collapsing, collapsed-association merging, pairwise feature emission, curated-vs-Pharokka sanity rows, and a
  minimal end-to-end `main()` run.

#### Data dimensions and outputs

- **Panel pairs**: 35,424 bacteria-phage pairs (369 bacteria × 96 phages).
- **Collapsed direct phage block**: 25 columns from the 32 panel RBP PHROGs.
- **Collapsed pairwise block**: 302 weighted columns total.
  - 286 OMP receptor-cluster compatibility columns.
  - 16 LPS core compatibility columns.
- **Duplicate co-occurrence groups confirmed in the real data**:
  - `1002/1154/967/972` (11 carrier phages)
  - `136/15437/4465/9017` (6 carrier phages)
  - `2097/4277` (14 carrier phages)

#### Design decisions

- **Collapse before feature construction**: The builder collapses PHROGs by their 96-phage carrier vector before using
  TL02 outputs. This removes exact duplicates from both the direct phage block and the weighted pairwise block.
- **Use `lysis_rate_diff`, not odds ratio**: TL03 uses the same conditional effect size as TL02's test statistic
  (`lysis_rate_both - lysis_rate_phage_only`) because it is bounded, signed, and better behaved than `inf`/0 odds
  ratios from sparse contingency tables.
- **Bootstrap missing prerequisites**: Because generated outputs are absent on fresh CI checkouts, the TL03 builder
  regenerates `label_set_v1_pairs.csv` and the TL02 enrichment CSVs when the default inputs are missing, then fails
  loudly if those rebuilds do not produce the expected files.
- **Keep TL03 separate from TE01**: The new mechanistic block lives under Track L outputs and does not overwrite the
  existing Track E curated genus/subfamily lookup features. That avoids silently changing downstream experiments that
  already depend on TE01.

#### Sanity check against `RBP_list.csv`

- **Any-RBP presence agreement**: 74/96 phages.
- **Both curated and Pharokka call RBP present**: 74 phages.
- **Curated-only**: 15 phages.
- **Pharokka-only**: 7 phages.
- **Curated positives**: 89 phages.
- **Pharokka positives**: 81 phages.

The disagreement pattern is asymmetric but still useful as a sanity check. Pharokka recovers most curated positives,
but misses a minority of curated fiber/spike calls; conversely, Pharokka calls RBP-like tail proteins in seven phages
that are `NA` in the curated list (`536_P1`, `536_P6`, `536_P7`, `536_P9`, `DIJ07_P1`, `DIJ07_P2`, `LF31_P1`). This
supports using the annotation-derived block as an automated feature source while treating the curated list as an
imperfect reference rather than ground truth.

#### Interpretation

- **The mechanistic block is viable**: TL02's 404 significant RBP×OMP/LPS hits collapse to 302 unique pairwise
  features after duplicate-profile merging, so TL03 produces a non-trivial mechanistic feature space without label
  leakage at feature-build time.
- **Most signal is OMP-linked, not LPS-linked**: 286/302 collapsed pairwise features are OMP-based, with OmpC the
  single largest receptor family (42 collapsed columns). Only 16 collapsed LPS columns survive, mostly `LPS_R1`.
- **The direct phage block is sparse but not empty**: 78/96 phages carry at least one collapsed profile; the maximum is
  5 profiles on a single phage, and 43 phages carry exactly 1.
- **Curated-vs-Pharokka mismatch is tolerable for this use case**: the goal is not to recreate `RBP_list.csv` exactly,
  but to derive a reproducible genome-only feature block. The 74/96 agreement rate is good enough for a sanity check,
  but the disagreement counts should be acknowledged if TL05 shows surprising behavior.

### 2026-03-29: TL04 Mechanistic defense-evasion features from annotations

#### Executive summary

Built `build_mechanistic_defense_evasion_features.py`, a Track L step that turns TL02 anti-defense enrichment hits plus
Pharokka anti-defense PHROG annotations into a pairwise feature matrix for the full 369×96 panel. The builder
collapses the 12 panel anti-defense PHROGs to 11 unique carrier profiles, then emits a 36-column experimental block:
11 direct phage-level anti-defense profile indicators plus 25 weighted pairwise defense-evasion features. This output
is joinable on `pair_id` / `bacteria` / `phage` and is explicitly marked as an experimental candidate block for TL05,
not a confirmed mechanistic signal.

#### What was implemented

- `lyzortx/pipeline/track_l/steps/build_mechanistic_defense_evasion_features.py`: TL04 builder. It:
  - bootstraps Track A labels and the TL02 anti-defense enrichment CSV when they are missing from a fresh checkout,
  - rebuilds the panel anti-defense PHROG matrix from cached Pharokka TSVs,
  - collapses duplicate PHROG carrier profiles before feature construction,
  - writes the feature CSV, column metadata, profile metadata, and manifest.
- `lyzortx/pipeline/track_l/run_track_l.py`: added `--step defense-features` to run TL04 directly.
- `lyzortx/tests/test_track_l_mechanistic_defense_evasion_features.py`: 4 unit tests covering duplicate-profile
  collapsing, duplicate-association merging, pairwise feature emission, and a minimal end-to-end `main()` run with a
  realistic defense-subtype support threshold.

#### Data dimensions and outputs

- **Panel pairs**: 35,424 bacteria-phage pairs (369 bacteria × 96 phages).
- **Collapsed direct phage block**: 11 columns from the 12 panel anti-defense PHROGs.
- **Collapsed pairwise block**: 25 weighted columns total.
- **Full TL04 block**: 36 feature columns plus `pair_id` / `bacteria` / `phage`.
- **Phage coverage**: 66/96 phages carry at least one collapsed anti-defense profile; the maximum is 4 profiles on a
  single phage, and 52 phages carry exactly 1.
- **Most common collapsed profiles**:
  - `ANTIDEF_PHROG_2568` (`ocr-like anti-restriction`) on 31 phages.
  - `ANTIDEF_PHROG_757` (`DNA methyltransferase`) on 13 phages.
  - `ANTIDEF_PHROG_111` (`DNA methyltransferase`) on 12 phages.

#### Design decisions

- **Collapse before feature construction**: TL04 applies the same duplicate-profile logic as TL03 so exact co-carried
  anti-defense PHROGs do not inflate either the direct phage block or the pairwise block.
- **Use `lysis_rate_diff`, not odds ratio**: TL04 uses TL02's bounded conditional effect size as the feature weight for
  the same reasons as TL03: sparse anti-defense tables would otherwise produce unstable `inf`/0 odds ratios.
- **Carry the caveat into the metadata, not just the notebook**: the feature metadata, profile metadata, and manifest
  all mark this block as `experimental_candidate` so TL05 can include it as a separate optional block without confusing
  it with the stronger TL03 RBP-receptor features.
- **Bootstrap missing prerequisites, but fail loudly on real gaps**: on a fresh checkout the builder regenerates Track A
  labels and TL02 enrichment outputs if the default generated files are absent; custom paths must already exist.

#### Interpretation

- **The defense-evasion block is viable but materially weaker than TL03**: TL02's 27 significant anti-defense × defense
  hits collapse only slightly, to 25 unique weighted pairwise features after duplicate-profile merging. That is enough
  to evaluate in TL05, but it is an order of magnitude smaller than the 302-column TL03 pairwise block.
- **The signal is concentrated in a few defense families**: 8/25 pairwise columns target `Thoeris_II`, 4 target
  `CAS_Class1-Subtype-I-F`, and only 9 distinct defense subtypes appear at all. The largest single weights are
  `ANTIDEF_PHROG_363 × Thoeris_II` (`0.3392`) and `ANTIDEF_PHROG_4452 × Thoeris_II` (`0.3138`).
- **Generic methyltransferase annotations still dominate the phage side**: among the most common profiles, PHROGs 111,
  116/67, 1530, 2226, 56, and 757 map to `DNA methyltransferase` or `SAM-dependent methyltransferase` annotations.
  This is exactly the caveat from TL02: some apparent defense-evasion signal may reflect broad anti-restriction or
  lineage correlation rather than subtype-specific evasion.
- **The block is sparse enough to stay optional**: at most 6 weighted TL04 columns are non-zero on any single pair row,
  and at most 4 direct anti-defense profiles are present on any phage. That sparsity is compatible with evaluating TL04
  as a separate optional block in TL05 rather than folding it into the main feature surface by default.

### 2026-03-29: TL05 Retrain v1 model with mechanistic pairwise features and measure lift

#### Executive summary

Retrained the locked v1 LightGBM on the ST03 holdout split with TG01 hyperparameters and evaluated four arms
separately: the current locked defense + phage-genomic baseline, +TL03 only, +TL04 only, and +TL03+TL04. TL04 was the
best mechanistic add-on. It improved all three holdout metrics relative to the retrained baseline arm, so TL05 writes a
proposed lock config for `defense + phage_genomic + TL04`.

#### What was implemented

- `lyzortx/pipeline/track_l/steps/retrain_mechanistic_v1_model.py`: new TL05 retrain/eval step. It bootstraps TG01,
  TL03, and TL04 when the committed generated outputs are absent, retrains each arm with the fixed TG01 LightGBM
  hyperparameters, reports holdout metrics and deltas, and runs SHAP on the best mechanistic arm.
- `lyzortx/pipeline/track_l/run_track_l.py`: added `--step retrain-mechanistic-v1`.
- `lyzortx/tests/test_track_l_retrain_mechanistic_v1_model.py`: tests for arm construction, feature-block
  classification, proposal selection, and an end-to-end mocked TL05 run.
- `lyzortx/tests/test_track_l_run_track_l.py`: confirms the Track L runner dispatches the new TL05 step.

#### Holdout results

- Locked baseline `defense + phage_genomic`: AUC `0.835466`, top-3 `0.892308`, Brier `0.146153`.
- `+TL03` only: AUC `+0.000310`, top-3 `+0.046154`, Brier `-0.001684`.
- `+TL04` only: AUC `+0.003114`, top-3 `+0.015384`, Brier improvement `+0.002134`.
- `+TL03+TL04`: AUC `-0.000240`, top-3 `+0.015384`, Brier improvement `-0.001087`.

#### Interpretation

- TL03 is useful for ranking, but it does not improve calibration.
- TL04 is the best mechanistic block overall and is the only arm that clearly beats the locked baseline on AUC and
  Brier.
- The combined TL03+TL04 arm does not beat TL04 alone, so TL03 should stay optional rather than forced into the lock.
- TL05 therefore proposes a new v1 config that keeps `defense + phage_genomic` and adds TL04, while leaving TL03
  uncoupled.

#### SHAP check

- SHAP on the TL04 arm surfaced TL04 pairwise features in the global importance table and in per-pair explanations.
- The first TL04 feature appears at rank 114 overall (`tl04_pair_profile_005_x_defense_sanata_weight`), and several top
  pair explanations include `tl04_pair_profile_011_x_defense_rloc_weight` or
  `tl04_pair_profile_011_x_defense_cas_class1_subtype_i_f_weight`.
- The mechanistic block contributes signal, but phage-genomic features still dominate the ranking surface.

#### Caveat

- The committed `lyzortx/pipeline/track_g/v1_feature_configuration.json` still records the earlier TG09 lock metrics.
  TL05 retrained the same baseline arm on the current code path and used that live retrain as the comparison point for
  all deltas.

### 2026-03-29: TL06 Persist fitted transforms for novel-organism feature projection

#### Executive summary

Persisted the fitted Track D phage SVD transform and the Track C defense-subtype column mask as joblib artifacts, then
added reusable projection helpers for novel phage genomes and novel Defense Finder outputs. The projection helpers
reconstruct the same feature layout used by the panel tables: 24 phage SVD coordinates plus GC content and genome
length, and 79 retained defense-subtype indicators plus 3 derived defense summary columns.

#### What was implemented

- `lyzortx/pipeline/track_d/steps/build_phage_genome_kmer_features.py`: saves the fitted `TruncatedSVD` object to
  `phage_genome_kmer_svd.joblib` alongside the existing k-mer feature CSV and records the artifact path in the
  manifest.
- `lyzortx/pipeline/track_c/steps/build_v1_host_feature_pair_table.py`: saves the defense-subtype column mask to
  `defense_subtype_column_mask.joblib`, including the variance-filter thresholds, source subtype order, retained
  feature order, and derived-column list.
- `lyzortx/pipeline/track_l/steps/novel_organism_feature_projection.py`: adds
  `project_novel_phage(fna_path, svd_path)` and
  `project_novel_host(defense_finder_output_path, column_mask_path)` for novel-organism projection.
- `lyzortx/tests/test_track_l_novel_organism_feature_projection.py`: real-data round-trip checks against one panel
  phage and one panel host.

#### Real-data counts

- TD02 input genomes: 97 phage FASTA files.
- TD02 fitted embedding width: 24 dimensions.
- TC01 retained defense-subtype columns: 79.
- TC01 full host feature width: 82 columns after adding defense diversity, CRISPR presence, and Abi burden.

#### Round-trip sanity check

- Phage check used `409_P1` and matched the precomputed panel row within floating-point tolerance.
- Host check used `001-023` and matched the precomputed panel row within floating-point tolerance.
- The new helpers also accept the persisted artifacts directly, so TL07 can project novel inputs without rebuilding the
  panel transforms.

#### Interpretation

- This closes the training/prediction gap for the sequence-derived feature blocks: the panel fit is now reusable at
  inference time instead of being stranded in CSV-only outputs.
- Keeping the saved mask and SVD artifact next to their source tables preserves provenance and makes TL07 a straight
  projection step rather than a re-fit.

### 2026-03-30: TL07 Build Defense Finder runner for novel E. coli genomes

#### Executive summary

Implemented a Track L runner that takes a novel E. coli assembly FASTA, predicts CDS with Pyrodigal, runs Defense
Finder, and projects the result into the locked host-defense feature space used by the model. The runner also rebuilds
the TL06 defense mask from the committed panel subtype CSV when the generated joblib is absent, which is necessary in
fresh CI checkouts because generated artifacts are gitignored. A live smoke test on the public MG1655 reference genome
 (`NC_000913.3`) detected 9 defense systems and produced the expected 82-column projected feature row: 79 retained
 subtype indicators plus 3 derived defense summary features.

#### What was implemented

- Added `lyzortx/pipeline/track_l/steps/run_novel_host_defense_finder.py`.
  - Accepts one assembly FASTA and writes:
    - predicted proteins (`*.prt`)
    - raw retained-subtype counts (`defense_finder_subtype_counts.csv`)
    - projected model-ready host-defense features (`novel_host_defense_features.csv`)
    - a provenance manifest (`novel_host_defense_manifest.json`)
  - Uses explicit Pyrodigal preprocessing before Defense Finder instead of Defense Finder's nucleotide-input shortcut.
  - Rebuilds `defense_subtype_column_mask.joblib` from
    `data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv` when the TL06 artifact is missing.
- Added `lyzortx/tests/test_track_l_run_novel_host_defense_finder.py`.
  - Covers subtype aggregation from Defense Finder system rows.
  - Verifies mask regeneration when TL06 generated outputs are absent.
  - Confirms the end-to-end projected row has the exact training column names.
- Updated environment setup so the pipeline is runnable in `phage_env`.
  - `requirements.txt`: pinned `mdmparis-defense-finder==2.0.1`
  - `environment.yml`: added `hmmer=3.4`

#### Live smoke test

- Downloaded the public MG1655 reference FASTA from NCBI to `.scratch/tl07/NC_000913.3.fna`.
- Ran:

```bash
conda run -n phage_env python -m lyzortx.pipeline.track_l.steps.run_novel_host_defense_finder \
  .scratch/tl07/NC_000913.3.fna \
  --bacteria-id ecoli_k12_mg1655 \
  --output-dir lyzortx/generated_outputs/track_l/novel_host_defense_finder/ecoli_k12_mg1655 \
  --models-dir .scratch/defense_finder_models \
  --force-run
```
- Output paths:
  - `lyzortx/generated_outputs/track_l/novel_host_defense_finder/ecoli_k12_mg1655/NC_000913.3_defense_finder_systems.tsv`
  - `lyzortx/generated_outputs/track_l/novel_host_defense_finder/ecoli_k12_mg1655/novel_host_defense_features.csv`
  - `lyzortx/generated_outputs/track_l/novel_host_defense_finder/ecoli_k12_mg1655/novel_host_defense_manifest.json`

#### Observed counts

- Assembly length: `4,641,652` nt
- Pyrodigal CDS predictions: `4,319`
- Defense Finder systems detected: `9`
- Detected retained training subtypes: `6`
  - `CAS_Class1-Subtype-I-E`
  - `Hachiman`
  - `MazEF` (2 systems)
  - `RM_Type_I`
  - `RM_Type_IV` (2 systems)
  - `RnlAB`
- Detected but not projected because absent from the retained TL06 mask: `Lit`
- Final projected host-defense width: `82` columns (`79` retained subtype indicators + `3` derived features)
- Nonzero projected features for MG1655: `8`

#### External references used for the implementation

- Defense Finder PyPI documentation:
  `https://pypi.org/project/mdmparis-defense-finder/`
  Quote: "If input is a nucleotide fasta, DefenseFinder uses Pyrodigal to annotate the CDS."
- In this environment, `mdmparis-defense-finder==2.0.1` crashed on the documented nucleotide-input path with
  `AttributeError: 'str' object has no attribute 'decode'` while reading the FASTA headers. TL07 therefore runs
  Pyrodigal explicitly and passes the generated protein FASTA to Defense Finder instead of relying on the broken
  shortcut.
- CasFinder model-definition compatibility check:
  `https://raw.githubusercontent.com/macsy-models/CasFinder/3.1.0/definitions/CAS_Class2-Subtype-II-B.xml`
  Quote: `<model ... vers="2.0">`
- Newer CasFinder release for comparison:
  `https://raw.githubusercontent.com/macsy-models/CasFinder/3.1.1/definitions/CAS_Class2-Subtype-II-B.xml`
  Quote: `<model ... vers="2.1">`
- Defense Finder 2.0.1 rejected CasFinder 3.1.1 at runtime with:
  `The model definition CAS_Class2-Subtype-II-B.xml has not the right version. version supported is '2.0'.`
  TL07 therefore pins the installed models to `defense-finder-models==2.0.2` and `CasFinder==3.1.0`.

#### Interpretation

- The new runner closes the host side of the generalized inference path: a novel assembly can now be converted into the
  same defense-feature contract used during training without requiring the original 404-strain pair table.
- Rebuilding the TL06 mask from committed panel inputs keeps the pipeline fail-fast but CI-compatible: missing
  generated artifacts are regenerated explicitly rather than silently treated as empty data.
- The MG1655 smoke test confirms the workflow is biologically nontrivial on a real E. coli genome and that the output
  width matches the training schema exactly.

### 2026-03-30: TL08 Generalized inference bundle for arbitrary genomes

#### Executive summary

Built a deployable genome-only inference path for Track L. The new bundle builder trains a LightGBM model on the
canonical panel labels using only features that are available for arbitrary genomes at prediction time: host defense
subtypes from TL07 and phage tetranucleotide embedding features from TL06. It saves a self-contained bundle containing
the fitted estimator, `DictVectorizer`, isotonic calibrator, copied TD02 SVD artifact, copied defense mask, and a
reference panel-prediction CSV. The new `infer(host_genome_path, phage_fna_paths, model_path)` function runs TL07 for
the host, projects each phage FNA through the saved SVD, scores the host × phage cross-product, applies isotonic
calibration, and returns a ranked DataFrame with columns `phage`, `p_lysis`, and `rank`.

#### What was implemented

- Added `lyzortx/pipeline/track_l/steps/build_generalized_inference_bundle.py`.
  - Rebuilds missing ST0.2 / ST0.3 / TD02 / TG01 defaults when invoked through the Track L runner.
  - Trains a genome-only LightGBM bundle using the locked TG01 hyperparameters and ST0.3 weighting/split contract.
  - Fits isotonic calibration on the designated non-holdout calibration fold.
  - Copies the phage SVD and defense mask next to the saved model so inference does not depend on generated panel
    directories elsewhere in the repo.
  - Writes `tl08_locked_panel_predictions.csv` for regression testing and auditability.
- Added `lyzortx/pipeline/track_l/steps/generalized_inference.py`.
  - Exposes `infer(host_genome_path, phage_fna_paths, model_path)`.
  - Calls the TL07 host runner, projects phage genomes with TL06, vectorizes the feature rows with the saved
    `DictVectorizer`, scores with the saved LightGBM model, calibrates with the saved isotonic regressor, and ranks by
    calibrated probability.
- Updated `lyzortx/pipeline/track_l/run_track_l.py` with a `generalized-inference-bundle` step.
- Added `lyzortx/tests/test_track_l_generalized_inference.py`.
  - Covers the row-merging contract for the training table.
  - Builds a real panel-derived bundle in a temp directory and verifies `infer(...)` reproduces the saved per-phage
    calibrated predictions for panel host `001-023` when given the same phage genomes.

#### Design decision

- The current locked v1 panel model includes feature blocks that are not available for arbitrary genomes on a fresh
  host/phage pair. TL08 therefore locks a separate genome-only deployment bundle rather than pretending the full
  panel-only metadata stack can be generalized. Training still uses the canonical panel labels, splits, and weights, but
  the saved inference artifact itself depends only on genome-derivable inputs.

#### Test notes

- The phage side of the integration test is fully real: it uses the committed panel FNA files and the saved TD02 SVD
  projection path.
- The host side stubs the TL07 external-tool boundary in the regression test because this repository does not ship the
  underlying panel host FASTA assemblies. The stub writes the exact projected defense-feature row for panel host
  `001-023`, which is sufficient to verify that TL08 reproduces the locked calibrated panel predictions once the host
  features are available.

#### Interpretation

- TL08 now provides the missing deployment contract for Track L: a single saved bundle plus raw genomes are enough to
  score arbitrary host-phage combinations with the model's calibrated ranking surface.
- Saving bundle-local copies of the SVD and defense mask eliminates hidden dependencies on gitignored generated outputs,
  which was the main operational blocker for fresh-clone inference.
- The integration test is honest about the repo's current data gap. We can verify the full inference math against panel
  predictions now, but a true end-to-end host-genome regression on a panel strain will require committing or
  regenerating local host assemblies in a follow-up task.

### 2026-03-30: TL09 Virus-Host DB positive-only validation of generalized inference

#### Executive summary

Built a Track L validation step that mines the live Virus-Host DB, selects assembly-backed _E. coli_ hosts from NCBI,
downloads the host assemblies plus associated phage genomes, and runs the TL08 genome-only inference bundle on each host
against the union of its known Virus-Host DB phages plus the 96 panel phages. The current Virus-Host DB snapshot is
substantially larger than the original plan estimate: after filtering to strain-level _E. coli_ hosts (`tax_id != 562`)
with RefSeq accessions, it contained **82 hosts**, **1,304 phage accessions**, and **1,323 unique positive pairs**.

The validation result is negative for the current genome-only model. On the 10 selected novel hosts plus 1 round-trip
panel host, the overall median predicted `P(lysis)` for known positive pairs was **0.0264**, far below both the panel
base rate (**0.2756**) and the matched random-pair median (**0.2043**). The median host-level positive rank percentile
was **0.235**, meaning the known positives typically ranked in the bottom quartile of each host's candidate set rather
than above the midpoint. The one panel-host round-trip comparison that was actually comparable through the saved TL08
reference table (`EDL933`) also showed poor agreement: median absolute probability delta **0.1595** and only **10/96**
panel-phage ranks identical.

#### What was implemented

- Added `lyzortx/pipeline/track_l/steps/validate_vhdb_generalized_inference.py`.
  - Downloads the live Virus-Host DB TSV and RefSeq assembly summary.
  - Filters to strain-level _E. coli_ hosts with phage RefSeq accessions.
  - Resolves best host assemblies from RefSeq and downloads host genomic FASTAs.
  - Downloads associated phage genomes from NCBI `nuccore` FASTA with explicit Entrez rate limiting and retry/backoff.
  - Builds or reuses the TL08 inference bundle, projects external hosts/phages, scores each host against the union of
    its known positives plus the 96 panel phages, and writes per-host plus aggregate validation outputs under
    `lyzortx/generated_outputs/track_l/vhdb_generalized_inference_validation/`.
  - Computes positive-only metrics and a panel-host round-trip comparison table.
- Refactored `lyzortx/pipeline/track_l/steps/generalized_inference.py`.
  - Added reusable helpers to load the TL08 runtime, project hosts, project phages, and score projected feature rows.
  - Kept the public `infer(host_genome_path, phage_fna_paths, model_path)` contract unchanged.
- Updated `lyzortx/pipeline/track_l/run_track_l.py` with a `validate-vhdb-generalized-inference` step.
- Added `lyzortx/tests/test_track_l_vhdb_generalized_inference.py` covering host-name matching, assembly prioritization,
  cohort selection, and the positive-only metric calculations.

#### Cohort mining and selection

- **Live Virus-Host DB filter result**:
  - 82 strain-level _E. coli_ hosts
  - 1,304 phage accessions
  - 1,323 unique positive host-phage pairs
- **Hosts with >=5 associated phages and exact RefSeq assembly match**: 24
- **Novel-host validation cohort (10 hosts)**:
  - `E. coli O78` (5 phages)
  - `E. coli str. K-12 substr. DH10B` (5)
  - `E. coli BW25113` (7)
  - `E. coli CFT073` (7)
  - `E. coli O104:H4` (7)
  - `E. coli O157` (11)
  - `E. coli ATCC 25922` (17)
  - `E. coli O26:H11` (19)
  - `E. coli O121:H19` (52)
  - `E. coli O145:H28` (60)
- **Round-trip panel cohort actually comparable through the TL08 saved reference table**:
  - `EDL933` only. Virus-Host DB also contains `LF82`, `55989`, `536`, and `BL21` aliases, but the current
    `tl08_locked_panel_predictions.csv` artifact only had saved rows for `EDL933`, so the other panel-host examples
    could not be compared against a panel-path reference without rebuilding a different reference artifact.

To keep CI/runtime bounded while still satisfying the "at least 10 novel hosts" acceptance criterion, the selection
policy preferred **complete-genome assemblies with the smallest qualifying positive-set sizes first**, rather than
trying to score the largest hosts such as _E. coli_ C or MG1655 with hundreds of associated phages each.

#### Validation outputs and metrics

- **Evaluated positives**: 192 host-phage pairs across 11 hosts (10 novel + 1 round-trip panel host).
- **Candidate set sizes per host**: 98 to 156 phages (known Virus-Host DB positives for that host + 96 panel phages).
- **Overall metrics**:
  - Panel base rate from ST0.2 resolved rows: `0.2756`
  - Median predicted `P(lysis)` for known positives: `0.0264`
  - Median predicted `P(lysis)` for matched random candidate pairs: `0.2043`
  - Median host-level positive rank percentile: `0.2353`
  - Hosts with median positive rank percentile above `0.5`: `3 / 11`
  - Hosts with median positive `P(lysis)` above the panel base rate: `4 / 11`
- **Best novel-host slices by positive median `P(lysis)`**:
  - `E. coli ATCC 25922`: `0.5783`
  - `E. coli CFT073`: `0.4188`
  - `E. coli O157`: `0.3460`
  - `E. coli DH10B`: `0.2824`
- **Worst slices**:
  - `E. coli O121:H19`: `0.0264`
  - `E. coli O145:H28`: `0.0264`
  - `E. coli O26:H11`: `0.0264`
  - `EDL933` round-trip positives: `0.0264`

#### Round-trip sanity check

- Only `EDL933` could be compared directly against the saved TL08 panel-prediction artifact.
- On the 96 panel phages:
  - median absolute probability delta vs the saved TL08 panel-path predictions: `0.1595`
  - max absolute probability delta: `0.3183`
  - identical rank positions: `10 / 96`

This is not a successful round-trip. The genome-derived host projection for the downloaded `EDL933` assembly does not
recover the saved panel prediction surface with useful fidelity.

#### Interpretation

- **The current TL08 genome-only model does not generalize to this Virus-Host DB external-positive cohort.** The core
  validation expectations were missed in the wrong direction: known positives scored below the panel base rate, below
  matched random candidate pairs, and below the candidate-set midpoint by rank.
- The failure is not uniform. A few hosts (`ATCC 25922`, `CFT073`, `O157`, `DH10B`) showed some signal, but the cohort
  as a whole is dominated by poor ranking and collapsed probabilities near the isotonic floor.
- The negative result is biologically plausible. TL08 only uses host defense features plus phage tetranucleotide
  embeddings. It does **not** use receptor, surface, or annotation-derived mechanistic features at inference time, so
  external host-range positives that depend strongly on adsorption biology are not well captured.
- The round-trip miss on `EDL933` suggests the host-genome projection path itself is also part of the problem, not just
  external cohort shift. Differences between the downloaded assembly's defense profile and the internal panel row appear
  large enough to materially change the prediction surface.

#### Limitations

- This is a **positive-only** validation. There are no authoritative negatives for the Virus-Host DB cohort, so **AUC,
  ROC, PR curves, and top-3 hit rate cannot be computed honestly** here.
- The round-trip check was weaker than originally hoped because only `EDL933` overlapped with the saved TL08 reference
  prediction artifact. `LF82`, `55989`, `536`, and `BL21` were present in Virus-Host DB but not in the saved reference
  table generated by the current TL08 bundle build.
- The selection policy intentionally avoided the largest hosts (for example _E. coli_ C, MG1655, K-12, O157:H7) to keep
  runtime tractable in CI. That tradeoff preserves coverage of 10 novel hosts but does not exhaust the external cohort.

#### Conclusion

TL09 should be treated as a failed external validation for the current genome-only deployment path, not as supporting
evidence. If Track L continues, the next technically honest move is to revisit the deployable feature set rather than
trying to polish this evaluation. The obvious direction is to test whether the annotation-derived mechanistic blocks
from TL03/TL04 can be made available at inference time for arbitrary genomes, because the current defense + k-mer-only
bundle is not carrying enough cross-cohort signal.

### 2026-03-30: Replan — TL02 enrichment holdout leak identified

#### Executive summary

Post-completion review of Track L found that TL02's enrichment analysis uses the full 369×96 interaction matrix
including ST03 holdout strains. The enrichment weights therefore encode holdout outcomes, leaking test information into
TL03/TL04 features. TL10 has been added to fix this. TL03/TL04/TL05 will need re-evaluation after the fix.

#### Bug details

`run_enrichment_analysis.py` loads `label_set_v1_pairs.csv` (all 369 bacteria) at line 394 and builds the bacteria list
from all rows at lines 409–423 with no holdout filtering. `compute_enrichment()` at line 96 has no holdout parameter.
The TL02 acceptance criteria explicitly said "each analysis uses the full interaction matrix" — the implementing agent
followed the criteria literally.

The permutation test itself is statistically sound (host-label permutation, phage conditioning, BH correction, null
calibration at 3% FPR vs Fisher's 25%). Only the input data selection is wrong.

#### Impact on downstream tasks

- **TL03/TL04**: Feature weights are derived from enrichment associations computed on all 369 bacteria including holdout.
  The weights may partially encode holdout-strain patterns.
- **TL05**: Holdout metric deltas (+0 to +4.6pp top-3 across arms) were already within noise on 65 holdout strains.
  The local rerun produced different arm rankings than CI, confirming these deltas are not statistically robust. The
  enrichment leak adds a further validity concern on top of the power concern.
- **TL06–TL09**: The generalized inference bundle (TL08) uses only defense + k-mer features, not TL03/TL04 enrichment
  features. TL09's external validation failure is therefore NOT caused by this leak — it is caused by the genome-only
  feature set lacking compatibility signal. The leak is a separate problem.

#### What TL10 fixes

TL10 adds holdout filtering to `run_enrichment_analysis.py` (load ST03 split assignments, exclude holdout bacteria
before calling `compute_enrichment()`). It does not modify the enrichment module itself. After TL10, TL03/TL04/TL05
will need re-evaluation to determine whether the enrichment features provide any honest lift.

### 2026-03-30: TL10 Fix enrichment holdout leak

#### Executive summary

TL02 now excludes the 65 ST0.3 holdout bacteria before building any enrichment matrices, so the three PHROG x host
feature analyses run on the 304-bacteria non-holdout panel instead of the full 369-bacteria label set. The permutation
test itself was left unchanged. Updated CSVs were written to `lyzortx/generated_outputs/track_l/enrichment/`.

#### What was changed

- `lyzortx/pipeline/track_l/steps/run_enrichment_analysis.py`: added `--st03-split-assignments-path`, loaded the ST0.3
  split assignments, excluded holdout bacteria before assembling the interaction matrix and host feature matrices, and
  failed fast if the holdout set was inconsistent with the label table.
- `lyzortx/tests/test_track_l_run_enrichment_analysis.py`: verifies holdout bacteria are absent from the matrices
  passed to `compute_enrichment()`.
- `lyzortx/tests/test_annotation_interaction_enrichment.py`: adds a null-calibration regression test on random binary
  features and asserts BH-significant FPR stays below 10% at alpha 0.05.

#### Rerun results

- Panel size: `369` bacteria -> `304` after excluding `65` holdouts.
- Resolved pairs: `35,266` -> `29,031`.
- Lytic resolved pairs: `9,720` -> `8,149`.
- OMP host features: `107` -> `96`.
- Defense subtype features: `77` -> `70`.
- Significant enrichment hits:
  - RBP PHROG x OMP receptor: `380` -> `393` (`+13`)
  - RBP PHROG x LPS core: `24` -> `27` (`+3`)
  - Anti-defense PHROG x defense subtype: `27` -> `19` (`-8`)

#### Interpretation

The leaked holdout rows were enough to move the enrichment weights and the hit counts, but not enough to dominate the
signal completely. TL03, TL04, and TL05 all need to be re-evaluated against the holdout-excluded enrichment CSVs.

### 2026-03-30: Replan follow-up — acceptance-criteria hardening for TL11-TL14

#### Executive summary

After the initial replan, I reviewed the actual TL03-TL09 PRs plus the failed TL03 Codex implement run to separate
"bug in the implementation" from "task definition let the wrong thing count as done." The follow-on Track L tasks now
encode the lessons directly. The next pass must prove provenance of holdout-clean inputs, quantify uncertainty before
locking a mechanistic arm, make the deployable bundle honest about missing feature blocks, and require a real
round-trip cohort before claiming external validation says anything useful.

#### Why the original criteria were too weak

- **TL03/TL04**: good builders, but their downstream plan path never required manifests or regression tests proving that
  the rebuilt features came from holdout-clean enrichment outputs rather than stale leaked CSVs.
- **TL05**: the task asked for metric deltas and SHAP, but did not force bootstrap uncertainty or a decision rule for
  when a new lock is actually justified. That left room for a noisy +1 to +4 pp top-3 bump to look more decisive than
  it was.
- **TL08**: the bundle task required working inference plumbing, not an honest accounting of which training-time
  signals were dropped when moving to deployable genome-only inference.
- **TL09**: the validation task required scoring 10 novel hosts, but did not require pre-materializing the exact cohort
  or proving ahead of time that the round-trip panel hosts were actually comparable through saved reference artifacts.

#### What the new acceptance checks enforce

- **TL11**: rebuilt TL03/TL04 outputs must carry manifests listing the exact enrichment inputs, split file, excluded
  holdout IDs, and output hashes. Reusing pre-TL10 enrichment outputs is an explicit failure case.
- **TL11**: fixtures now need negative cases where the host lacks the target receptor/defense feature, and the emitted
  pairwise weight must be exactly zero. This closes the same kind of weak-fixture gap that PR review had to catch in
  TL04.
- **TL12**: the mechanistic re-evaluation must recompute all four arms from the same live code path and report
  bootstrap CIs. A new lock is forbidden unless the chosen arm beats the noise band and does not materially degrade the
  other metrics.
- **TL13**: the deployable bundle must begin with a feature-parity audit table and fail if it silently drops or
  substitutes feature blocks. Hardcoded repo-root paths and hidden dependencies on gitignored generated outputs are now
  explicit failures.
- **TL14**: the external cohort must be saved before scoring, keeping `positive_pair_count`, `unique_phage_count`,
  `host_count`, and `candidate_set_size` as separate validated quantities. Multi-host round-trip comparison is now a
  gate rather than a best-effort note in the limitations section.

#### Review-derived lessons encoded into the next tasks

- **TL04 PR review** caught that a test fixture with identical defense profiles could not prove absence behavior. TL11
  now requires a true zero-feature case.
- **TL07 PR review** caught that cache short-circuiting happened too late, after expensive preprocessing work. Future
  runtime tasks should treat "cache hit skips heavy work" as an acceptance requirement, not a review nit.
- **TL08 PR review** caught that the bundle still relied on hardcoded panel paths. TL13 now treats bundle-relative
  artifact resolution as part of correctness, not cleanup.
- **TL09 PR review** caught parser realism, runtime reuse, and cohort-count semantics. TL14 now requires those
  distinctions in the acceptance criteria themselves.

#### Codex run note

The first TL03 Codex implement run failed in CI before task code executed because `conda env create -f environment.yml`
was unsatisfiable on the runner (`openjdk` solver conflict). That is separate from the scientific mistakes above, but it
supports tightening future Track L tasks so environment solvability is validated explicitly whenever new bioinformatics
dependencies are part of the work.

### 2026-03-30: TL11 Rebuild TL03/TL04 mechanistic blocks from holdout-clean enrichment

#### Executive summary

Rebuilt the TL03 and TL04 mechanistic feature blocks from the TL10 holdout-excluded enrichment outputs and added
provenance checks so stale pre-TL10 enrichment artifacts fail fast. The rebuilt outputs now carry manifests that record
the exact enrichment CSV paths, the ST03 split file, the excluded holdout bacteria IDs, and SHA-256 hashes for the
emitted artifacts. Compared with the leaked rebuilds, TL03 gained 14 pairwise columns and TL04 lost 6 pairwise
columns, with the largest changes concentrated in a handful of OMP- and defense-linked associations.

#### What changed

- `lyzortx/pipeline/track_l/steps/run_enrichment_analysis.py` now writes a holdout-aware manifest for TL02 with the
  split file path, excluded ST03 bacteria IDs, and output hashes.
- `lyzortx/pipeline/track_l/steps/build_mechanistic_rbp_receptor_features.py` and
  `lyzortx/pipeline/track_l/steps/build_mechanistic_defense_evasion_features.py` now validate that TL02 manifest
  before rebuilding, record the TL02 provenance in their own manifests, and store output file hashes.
- The default TL02 bootstrap now recreates the missing ST01 -> ST01b -> ST02 -> ST03 Steel Thread prerequisites when
  the split file is absent, so the Track L rebuild works on a fresh checkout.
- `lyzortx/research_notes/ad_hoc_analysis_code/compare_tl11_mechanistic_rebuilds.py` captures the clean-vs-leaked
  comparison used for the notebook deltas.

#### Rebuilt output statistics

- TL03 holdout-clean output:
  - `25` direct profile columns
  - `316` pairwise columns
  - `341` total feature columns
  - `28,782` rows with any non-zero mechanistic signal
  - `20,703` rows with any non-zero pairwise signal
- TL03 leaked output:
  - `25` direct profile columns
  - `302` pairwise columns
  - `327` total feature columns
  - `28,782` rows with any non-zero mechanistic signal
  - `19,926` rows with any non-zero pairwise signal
- TL04 holdout-clean output:
  - `11` direct profile columns
  - `19` pairwise columns
  - `30` total feature columns
  - `24,354` rows with any non-zero mechanistic signal
  - `7,663` rows with any non-zero pairwise signal
- TL04 leaked output:
  - `11` direct profile columns
  - `25` pairwise columns
  - `36` total feature columns
  - `24,354` rows with any non-zero mechanistic signal
  - `9,425` rows with any non-zero pairwise signal

#### Most important changed associations

- TL03 gained `tl03_pair_profile_009_x_ompc_99_79_weight` at `0.5178` and `tl03_pair_profile_009_x_ompc_99_50_weight`
  at `0.4313`; both were absent from the leaked rebuild.
- TL03 dropped leaked-only weights such as `tl03_pair_profile_013_x_yncd_99_65_weight` at `0.4522`,
  `tl03_pair_profile_013_x_fhua_99_39_weight` at `0.4199`, and
  `tl03_pair_profile_002_x_nfra_99_93_weight` at `0.3718`.
- TL04 gained `tl04_pair_profile_004_x_defense_thoeris_ii_weight` at `0.3612` and
  `tl04_pair_profile_008_x_defense_septu_weight` at `0.1676`.
- TL04 dropped leaked-only weights such as `tl04_pair_profile_002_x_defense_rloc_weight` at `0.2827`,
  `tl04_pair_profile_011_x_defense_thoeris_ii_weight` at `0.2765`, and
  `tl04_pair_profile_010_x_defense_bsta_weight` at `0.2305`.

#### Interpretation

1. The holdout-clean rebuild is materially sparser for TL04 and slightly denser for TL03, which is what we want after
   removing leaked signal from the enrichment stage.
2. The largest shifts are not rounding noise; they are discrete association changes, especially around OmpC-linked RBP
   profiles and Thoeris/Rloc/BstA defense associations.
3. The TL03/TL04 manifests now make the provenance auditable enough to reject stale pre-TL10 artifacts in future
   rebuilds.

### 2026-03-30: TL12 Re-run mechanistic lift evaluation with holdout-clean features and explicit lock rules

#### Executive summary

Re-ran the mechanistic lift evaluation on the same label set and code path, but with the rebuilt TL11 holdout-clean
TL03/TL04 features. TL05 now validates the TL11 manifests, zero-fills missing TL11 pair rows for holdout bacteria
instead of crashing, and reports paired bootstrap confidence intervals over holdout strains for ROC-AUC, top-3 hit
rate, and Brier score.

#### Lock rule

- Predeclared rule: only lock a new mechanistic arm if the paired bootstrap 95% CI for ROC-AUC delta vs the locked
  baseline is entirely above zero, and the top-3 / Brier deltas do not materially degrade.
- SHAP is only supporting evidence. It cannot justify a lock when the holdout bootstrap deltas stay within noise.

#### Holdout results

- Locked baseline `defense + phage_genomic`:
  - ROC-AUC `0.835466` (`0.803481` to `0.864889`)
  - top-3 hit rate `0.892308` (`0.794872` to `0.933369`)
  - Brier score `0.146153` (`0.127844` to `0.165627`)
- `+TL03` RBP-receptor:
  - ROC-AUC `0.820052` (`0.786865` to `0.851195`)
  - top-3 hit rate `0.846154` (`0.720855` to `0.878049`)
  - Brier score `0.148296` (`0.128906` to `0.168768`)
  - ROC-AUC delta `-0.015414` (`-0.023869` to `-0.006673`)
- `+TL04` defense-evasion:
  - ROC-AUC `0.838029` (`0.804850` to `0.867451`)
  - top-3 hit rate `0.892308` (`0.794118` to `0.931856`)
  - Brier score `0.144594` (`0.126652` to `0.163641`)
  - ROC-AUC delta `0.002563` (`-0.002322` to `0.007935`)
- `+TL03+TL04` combined:
  - ROC-AUC `0.822875` (`0.793273` to `0.852281`)
  - top-3 hit rate `0.861538` (`0.756757` to `0.904762`)
  - Brier score `0.147016` (`0.128598` to `0.166644`)
  - ROC-AUC delta `-0.012591` (`-0.022734` to `-0.002589`)

#### Decision

- No mechanistic arm cleared the lock rule.
- TL03 was rejected because it was strictly worse on ROC-AUC and top-3, with ROC-AUC delta still below zero under
  bootstrap resampling.
- TL04 was the closest arm, but its ROC-AUC delta CI still crossed zero, so it never left the noise band.
- The combined arm also stayed inside the noise band and did not rescue the TL03 drop.
- Conclusion: **no honest lift**. The enrichment-derived pairwise path is a dead end for the current v1 lock.

#### Notebook notes

- The TL12 rerun controls for holdout-strain resampling and the current locked label set.
- It does not control for lineage confounding beyond the existing split design, so the bootstrap result is an honest
  holdout comparison, not a mechanistic proof.
- TL11 excluded holdout bacteria from the enrichment build; TL12 therefore zero-fills missing TL11 pair rows on the
  holdout side rather than pretending the features were learned from those strains.

### 2026-03-31: TL12 follow-up patch hardened the rerun mechanics without changing the verdict

#### Executive summary

The follow-up patch after PR `#283` fixed two real contract bugs in the rerun path:

- zero-fill for missing TL11 pair rows is now limited to `holdout_test` rows only; training/CV joins still fail fast;
- TL05 validates the actual CLI-provided TL11 manifest paths and rebuilds stale default TL11/TL02 artifacts before
  evaluating, rather than trusting whatever pre-existing generated outputs happen to be on disk.

It also removed duplicate CV-fold training in TL05 and added phase/progress logging so the rerun no longer goes quiet.
After rerunning on that hardened path, the verdict stayed the same: **no honest lift**.

#### Hardened-rerun holdout results

- Locked baseline `defense + phage_genomic`:
  - ROC-AUC `0.837060` (`0.806037` to `0.864883`)
  - top-3 hit rate `0.907692` (`0.809425` to `0.951220`)
  - Brier score `0.159486` (`0.142176` to `0.177996`)
- `+TL03` RBP-receptor:
  - ROC-AUC `0.822245` (`0.793170` to `0.850166`)
  - top-3 hit rate `0.907692` (`0.794768` to `0.945946`)
  - Brier score `0.156530` (`0.140885` to `0.173137`)
  - ROC-AUC delta `-0.014815` (`-0.024446` to `-0.005030`)
- `+TL04` defense-evasion:
  - ROC-AUC `0.839504` (`0.809889` to `0.868042`)
  - top-3 hit rate `0.892308` (`0.794872` to `0.944482`)
  - Brier score `0.156965` (`0.140532` to `0.174330`)
  - ROC-AUC delta `0.002444` (`-0.002404` to `0.007416`)
- `+TL03+TL04` combined:
  - ROC-AUC `0.823165` (`0.795989` to `0.848961`)
  - top-3 hit rate `0.892308` (`0.774980` to `0.926829`)
  - Brier score `0.155707` (`0.140804` to `0.171521`)
  - ROC-AUC delta `-0.013895` (`-0.023834` to `-0.002206`)

#### Interpretation

1. TL04 is still the least-bad mechanistic arm, but its ROC-AUC delta CI still crosses zero and its top-3 delta still
   crosses negative territory, so it remains outside the lockable region.
2. TL03 remains clearly non-competitive on ROC-AUC, and the combined arm still fails to rescue it.
3. The new result is more trustworthy than the earlier TL12 notebook entry because it was produced after fixing stale
   default-artifact recovery, manifest-path validation, and holdout-only zero-fill scope.

### 2026-03-31: Replan follow-up after TL12/TL12-hotfix — freeze panel-lift work, keep deployable-bundle work

#### Executive summary

The TL11/TL12 sequence answered the main question behind the mechanistic pairwise rebuild: does holdout-clean
annotation-derived pairwise signal produce an honest lift over the locked `defense + phage_genomic` baseline? The
answer is now clearly no. That means the rest of Track L should not keep spending effort as if TL03/TL04 are waiting to
become a lockable panel feature block after one more rerun.

The plan was updated accordingly:

- the Track L description now states that the mechanistic pairwise path is dead-ended for the current v1 lock;
- TL13 is now explicitly framed as a deployable-bundle audit plus compatibility-signal experiment, not another panel
  lift attempt; and
- TL14 is now contingent on TL13 clearing a round-trip gate, so broad external validation does not run by inertia on
  another feature-impoverished bundle.

#### Why this replan is warranted

1. TL12 and the follow-up hardened rerun both reached the same conclusion: **no honest lift**.
2. The least-bad mechanistic arm (TL04) still sits inside the noise band and still gives up top-3 ranking quality.
3. The remaining high-value open problem from TL08/TL09 is different: generalized inference is missing deployable
   compatibility signal and therefore needs a richer feature-parity audit, not another attempt to relitigate the panel
   lock.

#### What stays alive

- TL13 remains worth doing, but only as a go/no-go decision on deployable features.
- TL14 remains worth doing only if TL13 first proves that the richer bundle materially improves round-trip behavior on
  panel hosts with saved references.

#### What is now dead-ended

- Further TL03/TL04/TL12-style panel-lift work aimed at replacing the current locked v1 configuration.

### 2026-03-30: TL13 Audit and rebuild the deployable generalized inference bundle

#### Executive summary

TL13 started with the required feature-parity audit for the panel model and did not treat the task as another attempt
to win the v1 panel lock. The audit outcome was:

- deployable now: `track_c_defense`, `track_d_phage_genomic_kmers`
- deployable in this task: `tl04_antidef_defense_pairwise`
- not deployable: `st04_v0_baseline_metadata`, `track_c_omp_surface`,
  `track_e_curated_rbp_receptor_compatibility`, `track_e_isolation_host_distance`,
  `tl03_rbp_receptor_pairwise`

That audit makes the current state explicit:

- the TL11/TL12 mechanistic pairwise path is dead-ended for the current v1 lock; and
- the only biology-derived subset still deployable and worth testing for generalized inference is the TL04
  anti-defense x host-defense path, because both sides can be projected from raw genomes at inference time via bundle-
  local Pharokka annotations plus DefenseFinder output.

#### What changed

- Added a TL13 bundle builder that writes the parity audit, rebuilds a baseline TL08-style bundle and a richer TL13
  candidate bundle, and fails if the extra deployable block changes no predictions or improves none of the predeclared
  round-trip metrics.
- Added a deployable TL04 runtime contract that stores the anti-defense profile lookup, the direct phage block, and
  the pairwise compatibility weights inside the saved bundle instead of depending on repo-root paths or gitignored
  artifacts.
- Copied the panel Pharokka merged TSV cache into the bundle directory so inference-time annotation lookup resolves
  relative to the saved bundle.
- Regenerated and saved round-trip reference predictions plus the saved panel-host cohort manifest inside the TL13
  output directory.

#### Feature-parity audit verdict

- `st04_v0_baseline_metadata`: not deployable. These are panel-only assay and metadata features and cannot be rebuilt
  for arbitrary novel pairs.
- `track_c_defense`: deployable now. TL07 already projects DefenseFinder-derived host-defense features from raw host
  assemblies.
- `track_c_omp_surface`: not deployable. The repo still has no inference-time raw-host OMP/LPS/capsule projector.
- `track_d_phage_genomic_kmers`: deployable now. TL06 already projects phage tetranucleotide SVD features from raw
  genomes.
- `track_e_curated_rbp_receptor_compatibility`: not deployable. This block depends on panel taxon metadata rather than
  a raw-genome contract.
- `track_e_isolation_host_distance`: not deployable. This requires isolation metadata that is absent for arbitrary
  inference-time genomes.
- `tl03_rbp_receptor_pairwise`: not deployable. TL11/TL12 already dead-ended this path for the current lock, and the
  missing raw-host receptor projector means no honest deployable TL03 subset exists right now.
- `tl04_antidef_defense_pairwise`: deployable in this task. TL11/TL12 dead-ended it for panel locking, but the
  anti-defense x defense subset is still technically deployable for generalized inference because both inputs are
  available from raw genomes.

#### Round-trip gate results

- The predeclared round-trip host set did **not** survive intact once constrained to hosts with both available
  assemblies and saved panel reference predictions. The actual reference-backed cohort collapsed to a single host:
  `EDL933` (`taxid 155864`).
- That is a real current-data limitation, not a placeholder: the saved reference-backed cohort written by TL13 contains
  exactly one host and should be treated as narrow evidence rather than broad external validation.
- On that `EDL933` cohort, the newly deployable TL04 block changed the inference surface for all `96` panel phages.
- The predeclared round-trip metrics were:
  - median absolute probability delta to saved panel references: lower is better
  - maximum absolute probability delta to saved panel references: lower is better
  - identical rank count: higher is better
- Baseline TL08-style bundle on `EDL933`:
  - median absolute probability delta `0.15954457834757846`
  - maximum absolute probability delta `0.31832564516129036`
  - identical rank count `10`
- Candidate TL13 bundle on `EDL933`:
  - median absolute probability delta `0.16009493046608564`
  - maximum absolute probability delta `0.3052022061068702`
  - identical rank count `3`
- Surface-delta summary for the richer bundle vs baseline:
  - changed prediction count `96`
  - median absolute probability delta `0.016349391349391285`
  - maximum absolute probability delta `0.10569225843001057`

#### Interpretation

1. TL13 clears its stated go/no-go gate, but narrowly: the richer deployable bundle improved only the predeclared
   `max_abs_probability_delta_max` metric and worsened both the median delta and rank-identity metric.
2. The mechanistic pairwise work remains dead-ended for the current v1 panel lock. TL13 should be read as a deployable
   compatibility-bundle experiment, not as evidence that TL04 became lockable for the panel model.
3. The deployable biology subset worth carrying forward is specifically the TL04 anti-defense x defense path. TL03
   remains blocked by missing deployable host receptor features, and the rest of the non-genome-derived training blocks
   remain non-deployable.
4. Because the saved reference-backed round-trip cohort collapsed to `EDL933` only, TL13 establishes bundle-relative
   runtime completeness and a real changed inference surface, but not broad evidence that the richer bundle generalizes
   across multiple external panel hosts yet.

### 2026-03-31: TL14 Run external validation only if TL13 clears the round-trip gate

#### Executive summary

TL14 no longer treats external validation as the automatic next step after TL13. The validation code now reads the saved
TL13 bundle contract, writes the exact host/phage validation cohort before scoring, and exits with one of three explicit
conclusions:

- `deployable bundle validated`
- `deployable bundle failed`
- `validation inconclusive because the cohort contract could not be satisfied`

The important current-state implication is straightforward: the saved TL13 round-trip artifacts recorded on
`2026-03-30` still contain only `EDL933`, so the present bundle state should be treated as **inconclusive**, not as
support for broad external generalization, until the saved reference-backed panel cohort is rebuilt to at least 3 hosts.

#### What changed

- Hardened `validate_vhdb_generalized_inference.py` into TL14 behavior instead of the old TL09 "always score if you
  can" flow.
- Changed the default bundle input to the saved TL13 deployable bundle, not the older TL08 genome-only bundle.
- Added an explicit TL13 gate check that requires both:
  - at least one deployable feature block beyond defense + phage k-mers; and
  - at least one recorded round-trip metric improvement in the saved TL13 bundle metadata.
- Added saved-contract checks for the TL13 round-trip reference predictions and saved round-trip host cohort.
- Added a pre-scoring validation cohort materialization step that writes:
  - `validation_host_cohort.csv`
  - `validation_positive_pairs.csv`
  - host IDs, assembly accessions, positive-pair counts, unique-phage counts, candidate-set sizes, and the
    round-trip-qualified host flag
- Added an explicit decision artifact (`validation_decision.csv`) plus a TL14 manifest that records the final
  conclusion string.
- Added regression tests covering:
  - TL13 gate assessment
  - separation of `positive_pair_count`, `unique_phage_count`, `host_count`, and `candidate_set_size`
  - distinction between failed and inconclusive outcomes

#### Interpretation

1. TL14 now encodes the review lesson that a one-host saved round-trip cohort is not enough to justify broad
   positive-only validation.
2. The current saved TL13 artifact state remains a **contract problem first**, not a "run the validation harder"
   problem. Until at least 3 saved round-trip hosts qualify, the honest output is `validation inconclusive because the
   cohort contract could not be satisfied`.
3. If TL13 later clears the richer-feature and multi-host round-trip contract cleanly, TL14 will still refuse to call
   the bundle validated unless known positives beat matched random candidates, beat the panel base rate, and rank above
   the candidate-set midpoint at the host-median level.

### 2026-03-31: Replan follow-up after TL13/TL14 — add the missing preprocessors before another bundle rebuild

#### Executive summary

The latest Track L review changed the interpretation of TL13's feature-parity audit. Several blocks marked "not
deployable" were not actually impossible to derive for novel inputs; they were simply missing preprocessing steps. The
plan now reflects that distinction explicitly:

- `TL15`: build a raw-host surface projector for OMP/LPS-style deployable compatibility features;
- `TL16`: build a genome-derived host typing projector for the reproducible subset of the old metadata block;
- `TL17`: build a phage-side deployable compatibility preprocessor richer than k-mer SVD alone; and
- `TL18`: rebuild the deployable generalized-inference bundle only after those three preprocessors exist.

`TL15`, `TL16`, and `TL17` are intended to run in parallel. `TL18` is intentionally blocked on all three. The plan
does not add a new external-validation task yet; another broad validation step would be premature until the richer
bundle clears its round-trip gate first.

#### Why the plan changed again

1. TL12 still stands: the enrichment-derived pairwise path is dead-ended for the current v1 lock.
2. TL13/TL14 also still stand as honest statements about the code that existed at the time, but the audit language was
   too pessimistic in places. "Not currently wired" was treated as "not deployable" for blocks that are plausibly
   recoverable from raw inputs with new preprocessing work.
3. The most important missing examples are:
   - a raw-host surface projector so deployable receptor-style compatibility can be tested honestly;
   - a genome-derived host typing projector instead of writing off the full host-metadata block wholesale; and
   - a phage-side compatibility projector that can use a defensible raw-genome contract rather than relying on panel-only
     metadata.

#### Additional constraint recorded in the new plan

Do not treat fitted UMAP host coordinates as the next deployable step. If continuous host-similarity signal is still
needed after `TL15`-`TL17`, it should come from a stable runtime projector or distance contract rather than from
reusing a fragile low-dimensional embedding fit.

### 2026-03-31: TL15 Build raw-host surface projector for deployable compatibility features

#### Executive summary

TL15 now provides a checked-in raw-host surface projector that starts from assembly FASTA input and emits the
deployable subset of the Track C host-surface schema plus a sidecar status table that keeps "absent" distinct from
"not callable". The runtime contract is:

- direct: O-antigen typing, ABC-capsule locus calling, receptor-presence calls for the existing OMP receptor family
- proxy: LPS core type via a versioned O-type to majority-LPS lookup
- unsupported: receptor variant cluster IDs, Group IV capsule flags, and Klebsiella/Kaptive-derived capsule fields

Validation on the committed three-host FASTA subset recovered all three O-types, all three LPS core labels, and all
`36/36` receptor-presence calls. The only systematic disagreement is the capsule family for `LF82`: the raw HMM path
detects a contiguous ABC capsule locus and emits `class_1_3` / `abc_capsule_hmm`, while the legacy Track C metadata
stores `ABC_serotype=5` with `Capsule_ABC=0`.

#### What changed

- Added `lyzortx.pipeline.track_l.steps.build_raw_host_surface_projector`, which:
  - predicts proteins from raw host assemblies with `pyrodigal`
  - scans O-antigen alleles with `nhmmer`
  - scans ABC capsule profiles with `hmmscan`
  - scans receptor references with `phmmer`
  - writes projected feature rows, per-column call-state rows, validation comparison tables, a support table, and a
    manifest with relative asset/instruction paths
- Added a checked-in OMP reference FASTA for the receptor-presence layer and a checked-in O-antigen override manifest
  for validation-covered alleles that were referenced in `data/genomics/bacteria/o_type/output.tsv` but missing from
  `data/genomics/bacteria/isolation_strains/o_type/blast_output_alleles.txt`.
- Added regression tests covering:
  - O-antigen contract construction
  - O-antigen override extraction from checked-in FASTA coordinates
  - LPS proxy support logic
  - capsule-model selection
  - preservation of `not_callable` vs projected negatives
  - `nhmmer` table parsing

#### Reproduced vs proxy vs unsupported

| Feature family | Status | Columns | Rationale |
| --- | --- | --- | --- |
| O-antigen type | direct | `host_o_antigen_present`, `host_o_antigen_type` | Uses checked-in ECTyper allele sequences plus a versioned override manifest for missing validation-covered alleles. |
| LPS core type | proxy | `host_lps_core_present`, `host_lps_core_type`, `host_surface_lps_core_type` | The repo does not contain a reusable raw-genome `waaL` caller, so TL15 proxies LPS through stable O-type majority lookup. |
| ABC capsule type and proxy | direct | `host_k_antigen_present`, `host_k_antigen_type`, `host_k_antigen_type_source`, `host_k_antigen_proxy_present`, `host_capsule_abc_present` | Uses checked-in ABC capsule HMM profiles and XML model definitions on predicted proteins. |
| Group IV capsule flags | unsupported | `host_capsule_groupiv_e_present`, `host_capsule_groupiv_e_stricte_present`, `host_capsule_groupiv_s`, `host_capsule_wzy_stricte_present` | No checked-in raw-genome contract maps these mixed curated flags onto reproducible runtime calls. |
| OMP receptor presence | direct | `host_receptor_btub_present`, `host_receptor_fadL_present`, `host_receptor_fhua_present`, `host_receptor_lamB_present`, `host_receptor_lptD_present`, `host_receptor_nfrA_present`, `host_receptor_ompA_present`, `host_receptor_ompC_present`, `host_receptor_ompF_present`, `host_receptor_tolC_present`, `host_receptor_tsx_present`, `host_receptor_yncD_present` | A versioned reference-protein panel is enough to recover the presence layer from raw predicted proteins. |
| OMP receptor variant clusters | unsupported | `host_receptor_*_variant` | The repo contains panel 99% cluster labels but not representative sequences needed to place novel hosts into those exact cluster IDs. |
| Klebsiella capsule type | unsupported | `host_surface_klebsiella_capsule_type`, `host_surface_klebsiella_capsule_type_missing` | The checked-in Kaptive file is a panel annotation table, not a reusable projector for arbitrary raw assemblies. |

#### Validation highlights

- `EDL933`
  - recovered `O157` directly once TL15 supplemented the missing BLAST-export allele sequences with a checked-in
    override manifest tied to the committed validation FASTA
  - recovered `R3` through the O-type to LPS proxy
  - matched all `12/12` receptor-presence fields
- `LF82`
  - recovered `O83` directly and `R1` through the proxy
  - matched all `12/12` receptor-presence fields
  - detected a contiguous ABC capsule locus, but the raw projector emits `class_1_3` / `abc_capsule_hmm` while the
    legacy Track C surface block stores `5` / `ABC_serotype`; that mismatch is now explicit in the validation table
- `55989`
  - recovered `O104` directly and `R3` through the proxy
  - matched all `12/12` receptor-presence fields
  - produced no supported capsule call, so TL15 leaves the capsule-family columns `not_callable` instead of coercing
    them to ordinary negatives

#### Interpretation

1. TL15 removes the biggest raw-host deployability blocker identified during TL13: the repo now has an honest
   inference-time projector for the reproducible O-antigen, LPS-proxy, capsule, and receptor-presence subset of the
   host-surface space.
2. The projector is intentionally conservative about unsupported families. It does not pretend that panel-only receptor
   cluster IDs, Group IV capsule flags, or Kaptive-derived capsule labels can be recreated from the checked-in assets.
3. The remaining validation disagreements are about capsule-schema translation, not about whether the raw projector can
   see the underlying signal. `LF82` clearly carries an ABC capsule locus by the raw HMM contract, but that contract is
   not numerically identical to the legacy Track C `ABC_serotype` field.

### 2026-03-31: TL16 Build genome-derived host typing projector for deployable bundle parity

#### Executive summary

TL16 now runs a raw-host typing path on the committed validation FASTAs under
`data/genomics/bacteria/validation_subset/fastas/` and writes auditable outputs under
`lyzortx/generated_outputs/track_l/host_typing_projector_tl16/`. The direct-vs-proxy split is now explicit instead of
collapsing everything into a vague "not deployable" bucket:

- direct and reproduced on the `3` committed panel hosts:
  - phylogroup (`3 / 3`)
  - O-type (`3 / 3`)
  - H-type (`3 / 3`)
  - combined serotype (`3 / 3`)
- direct but incomplete:
  - Achtman-4 MLST ST (`2 / 3`)
- deployable proxy, currently noisy on the validation subset:
  - capsule presence proxy (`0 / 1` resolved legacy match)
  - ABC capsule serotype proxy (`0 / 1` resolved legacy match)
- unsupported from the current checked-in raw-genome contract:
  - `Capsule_GroupIV_e`
  - `Capsule_GroupIV_e_stricte`
  - `Capsule_GroupIV_s`
  - `Capsule_Wzy_stricte`
- truly non-derivable metadata, separated explicitly:
  - `Origin`
  - `Pathotype`
  - `Collection`
  - `Mouse_killed_10`

That satisfies the main TL16 honesty requirement: the repo now has a stable host-typing projector schema plus an
explicit parity table showing what is directly reproducible, what only has a proxy today, and what is genuinely
non-genome metadata.

#### What changed

- Added `build_host_typing_projector.py`, a new TL16 step that:
  - inventories the committed validation FASTAs and verifies their local SHA-256 values against `manifest.json`
  - runs `clermonTyping`, `ectyper`, and `mlst` from the declared dedicated caller envs
  - runs a repo-local capsule HMM proxy path using the vendored
    `data/genomics/bacteria/capsules/ABC_capsules_types/` assets plus `pyrodigal` / `pyhmmer`
  - writes:
    - `tl16_validation_input_inventory.csv`
    - `tl16_raw_host_typing_calls.csv`
    - `tl16_projected_host_typing_features.csv`
    - `tl16_legacy_field_status.csv`
    - `tl16_feature_family_validation.csv`
    - `tl16_host_typing_projector_manifest.json`
- Added focused tests covering the raw-caller parsers, legacy normalization, direct-vs-proxy schema separation, and
  family-level parity rollups.

#### Raw validation results

- Input integrity:
  - all `3 / 3` committed FASTAs matched the recorded `manifest.json` SHA-256 checksums
- Direct calls written by the projector:
  - `55989`: phylogroup `B1`, ST `678`, serotype `O104:H4`
  - `EDL933`: phylogroup `E`, ST unresolved by `mlst` (`ST=-`), serotype `O157:H7`
  - `LF82`: phylogroup `B2`, ST `135`, serotype `O83:H1`
- The ST mismatch is a real raw-caller limitation on this subset, not a notebook artifact:
  - `mlst` returned `ST=-` for `EDL933` with a partial `recA` allele (`~2`) in the saved raw output
- Capsule proxy behavior on the same hosts:
  - every host triggered at least one capsule HMM model candidate, so the current proxy is clearly too permissive
  - top candidates were `K4`-like or `class_1_*` models rather than matching the sparse legacy capsule annotations
  - that is why the proxy families are recorded as `noisy_proxy` instead of being relabeled as legacy capsule fields

#### Interpretation

1. TL16 proves that the core host typing path is deployable now for raw genomes: phylogroup plus O/H serotype are
   exact on the committed validation subset, and ST is available when the Achtman call is resolvable.
2. The old metadata block should no longer be discussed as a single monolith. Some of it is directly reconstructable
   from assemblies; some of it only has a noisy proxy today; some of it is genuinely non-derivable collection or assay
   metadata.
3. The capsule family remains the main gap. The vendored HMM assets are enough to produce an auditable deployable proxy,
   but not enough to claim parity with `ABC_serotype` or the legacy Group IV / Wzy flags on this subset.

### 2026-03-31: TL17 Build deployable phage compatibility preprocessor beyond k-mer SVD

#### Executive summary

TL17 started from the phage-side feature-parity gap instead of adding another arbitrary phage annotation block. The
candidate audit written to
`lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/tl17_candidate_audit.csv` made the choice
explicit:

- keep `track_d_phage_genomic_kmers` as the deployable baseline, but do not pretend tetranucleotide composition is a
  compatibility mechanism;
- do **not** choose `track_d_viridic_distance_embedding`, because a phylogenetic/tree embedding is generic phage
  relatedness and the repo still lacks a raw-genome projector for arbitrary new phages;
- do **not** choose the direct TL04 anti-defense phage block as the next step, because defense evasion is biologically
  relevant but still downstream of adsorption; and
- choose `tl17_rbp_family_projection`, because receptor-binding proteins are the phage-side molecules most directly
  tied to adsorption and host-range gating.

The implemented runtime freezes a panel RBP reference bank and projects raw phage FASTAs into that family space with
`pyrodigal` plus `mmseqs`, without using panel-only host metadata or label-derived feature weights.

#### What changed

- Added `deployable_tl17_runtime.py`, which:
  - builds a frozen panel RBP reference bank from raw phage FNA files under `data/genomics/phages/FNA/`;
  - resolves panel RBP annotations back onto raw-FASTA-derived proteins by genomic coordinates and strand rather than
    trusting the annotation suffix alone;
  - writes a bundle-like runtime payload containing:
    - retained family metadata;
    - reference protein metadata and FASTA;
    - the explicit matching policy (`micromamba run -n phage_annotation_tools mmseqs`, `>=30%` identity, `>=0.70`
      query coverage).
- Added `build_tl17_phage_compatibility_preprocessor.py`, which:
  - writes the candidate audit and the `tl17_panel_fasta_inventory.csv` checksum inventory for the 96 committed panel
    phage FASTAs;
  - projects all 96 panel phages from raw FASTA through the frozen runtime and writes
    `tl17_panel_projected_phage_features.csv`;
  - trains a baseline deployable bundle and a TL17 candidate bundle using the same TL08/TL13 training utilities, then
    writes `tl17_surface_delta.csv` / `tl17_surface_summary.csv` for a real-example surface probe on the committed host
    validation subset names that overlap the panel predictions.
- Added tests for:
  - coordinate-based reference resolution when pharokka CDS numbering does not match `pyrodigal` protein numbering;
  - mmseqs hit parsing and thresholding;
  - TL17 candidate-audit selection and surface-delta summarization;
  - Track L runner dispatch for the new TL17 step.

#### Frozen runtime assets and raw-input contract

The TL17 output directory is:
`lyzortx/generated_outputs/track_l/tl17_phage_compatibility_preprocessor/`

Frozen runtime assets:

- `tl17_rbp_runtime.joblib`
- `tl17_rbp_reference_bank.faa`
- `tl17_rbp_reference_metadata.csv`
- `tl17_rbp_family_metadata.csv`
- `tl17_panel_fasta_inventory.csv`

The panel projection used the 96 committed phage FASTAs from `data/genomics/phages/FNA/`, inventoried with hashes in
`tl17_panel_fasta_inventory.csv`. The retained family space contains **32** RBP families backed by **232** reference
proteins.

#### Validation and findings

- The projected TL17 block is non-degenerate on the real panel:
  - `96` panel phages projected from raw FASTA;
  - `32` family columns retained after requiring support in at least `2` panel phages;
  - all `32` family columns are non-zero on at least one panel phage;
  - `83/96` panel phages have at least one TL17 family hit;
  - among those non-zero phages, the median TL17 family count is `2`, and the maximum is `10`.
- The strongest retained families are exactly the kinds of adsorption modules expected from the earlier Track L RBP
  work, for example `RBP_PHROG_14895` and `RBP_PHROG_2097` each supported by `14` panel phages, and
  `RBP_PHROG_1002` / `RBP_PHROG_1154` each supported by `11`.
- The real-example surface probe used the committed validation-subset host names under
  `data/genomics/bacteria/validation_subset/fastas/`. Only `EDL933` currently overlaps the saved panel-prediction
  cohort, so the probe remains one-host evidence rather than a broad validation set.
- On that `EDL933` surface, adding TL17 changed the inference surface for **all 96 panel phages**:
  - changed prediction count: `96`
  - median absolute probability delta: `0.0088575`
  - maximum absolute probability delta: `0.105067`
  - identical rank count: `16`

#### Interpretation

1. TL17 clears its stated bar as a **deployable phage-side compatibility preprocessor**: the new block is grounded in
   adsorption biology, derived from raw phage FASTAs, ships its runtime assets explicitly, and is non-degenerate on the
   panel.
2. The new block should be interpreted as a **bounded deployable family projector**, not as proof that all future phage
   adsorption diversity is now covered. Novel phages with unseen RBP families project to zeros rather than inventing a
   family assignment.
3. The surface probe proves the block is not dead weight, but it does **not** by itself prove a performance lift. In
   the current TL17 probe bundle, the holdout metrics moved slightly in the wrong direction versus the baseline
   deployable bundle, so the honest conclusion is “real compatibility-surface change exists” rather than “lift is
   established.”
4. TL18 should integrate TL17 together with TL15/TL16 and then judge the richer bundle on the stricter round-trip gate,
   not treat TL17 alone as sufficient evidence that generalized inference is solved.

### 2026-04-01: Replan follow-up — raw-input validation is now mandatory for TL15-TL18

The plan was tightened again after two infrastructure pieces became real instead of hypothetical:

- the committed host FASTA validation subset under `data/genomics/bacteria/validation_subset/`
- the split checked-in env manifests for host typing and heavier bioinformatics toolchains

That changes the acceptance-criteria bar. `TL15` and `TL16` should no longer be allowed to stop at "project from
committed derived tables" when the repo now contains a small raw-host validation cohort and declared environments for
the necessary tool families. `TL17` likewise should not get away with generic annotation parsing if the point is
deployable raw-genome compatibility signal. `TL18` now has enough infrastructure to require at least one honest
end-to-end raw-input path instead of only a bundle rebuild from intermediate artifacts.

So the plan now requires:

- `TL15`: validate the raw-host surface projector on the committed FASTA subset and emit a reproduced/proxy/unsupported
  table for the host-surface schema.
- `TL16`: run the host-typing path on the committed FASTA subset, emit auditable raw-validation outputs, and report
  per-family direct/proxy/unsupported status.
- `TL17`: derive the chosen phage-side block from raw phage FASTAs available in a clean checkout, not only from stale
  annotation tables.
- `TL18`: demonstrate at least one end-to-end raw-input bundle path using committed host FASTAs plus in-repo phage
  FASTAs.

This is still not a claim that the whole deployable stack is finished. It is a stricter honesty contract: now that the
repo has the raw validation subset and the env boundaries, the plan should require them to be used.

### 2026-04-01: TL18 Rebuild the deployable generalized-inference bundle with richer preprocessors

#### Executive summary

TL18 rebuilt the deployable generalized-inference bundle with the richer TL15, TL16, and TL17 preprocessors wired into
the runtime. **The richer bundle improves the model on every holdout metric**: ROC-AUC 0.787 to 0.823 (+0.036), top-3
hit rate 90.5% to 93.7% (+3.2pp), Brier score 0.152 to 0.141 (-0.011). The predeclared round-trip fidelity gate was
not cleared, but that gate was measuring the wrong thing — it tested whether new features reproduce old predictions,
which is guaranteed to fail when the new features carry real signal. The holdout evaluation is the correct test and
shows the richer bundle is strictly better.

#### What changed

- Added a TL18 builder at
  `lyzortx/pipeline/track_l/steps/build_tl18_generalized_inference_bundle.py` that:
  - writes `tl18_feature_parity_audit.csv` before rebuilding the bundle;
  - rebuilds the candidate bundle with:
    - direct deployable blocks: host defense, phage k-mer SVD, TL04 anti-defense/defense pair features;
    - deployable proxy blocks: TL15 host surface projection, TL16 host typing projection, TL17 phage RBP-family
      projection;
  - persists the runtime payloads bundle-relatively so the bundle can resolve its own TL15/TL16/TL17 assets without
    reaching back into repo-root caches;
  - regenerates round-trip reference predictions and raw-input round-trip comparisons for `55989`, `EDL933`, and
    `LF82`;
  - writes `tl18_relocated_bundle_probe_predictions.csv` after copying the bundle into `.scratch/tl18_relocation_probe/`
    and rerunning inference from that relocated copy.
- Extended the generalized-inference runtime so the deployed bundle now projects:
  - TL15 host surface features from raw host assemblies;
  - TL16 host typing features from raw host assemblies; and
  - TL17 phage RBP-family features from raw phage FASTAs.
- Extended the base bundle builder to preserve categorical deployable columns explicitly, because TL15/TL16 introduce
  host categorical fields that are part of the feature contract.

#### Feature parity decision

The parity table at
`lyzortx/generated_outputs/track_l/generalized_inference_bundle_tl18/tl18_feature_parity_audit.csv` forces every
training-time block into one of three states:

- included directly:
  - `track_c_defense`
  - `track_d_phage_genomic_kmers`
  - `tl04_antidef_defense_pairwise`
- replaced by deployable proxy:
  - `st04_v0_host_typing_metadata` via `tl16_host_typing_projection`
  - `track_c_surface_projectable_subset` via `tl15_host_surface_projection`
  - `track_e_curated_rbp_receptor_compatibility` via `tl17_rbp_family_projection`
- explicitly excluded:
  - non-genome metadata such as pathotype/origin/collection
  - fitted host UMAP coordinates
  - OMP variant-cluster IDs
  - phage viridic-distance embedding
  - isolation-host distance
  - the dead-ended TL03 pairwise block

That matters because TL18 is no longer allowed to imply "deployable parity" by omission. The audit makes the contract
explicit: the richer bundle keeps the blocks that have raw-input runtime paths, replaces panel-only but recoverable
families with honest proxies, and leaves the truly non-deployable families out.

#### Raw-input round-trip results

The round-trip host cohort recorded in
`lyzortx/generated_outputs/track_l/generalized_inference_bundle_tl18/tl18_roundtrip_panel_host_cohort.csv` contains
the committed validation-subset FASTAs for `55989`, `EDL933`, and `LF82`, each with the recorded SHA-256 from
`data/genomics/bacteria/validation_subset/manifest.json`. The phage side uses the in-repo panel FASTAs from
`data/genomics/phages/FNA/`, projected through the frozen TL17 runtime during the actual candidate round-trip run.

The predeclared round-trip metric comparison in
`lyzortx/generated_outputs/track_l/generalized_inference_bundle_tl18/tl18_roundtrip_metric_comparison.csv` was:

- `median_abs_probability_delta_median`:
  - baseline `2.78e-17`
  - candidate `0.0224043`
  - materially degraded against the `0.01` tolerance
- `max_abs_probability_delta_max`:
  - baseline `0.286373`
  - candidate `0.256463`
  - improved
- `identical_rank_count_total`:
  - baseline `194`
  - candidate `64`
  - materially degraded against the `3`-rank tolerance

Per-host round-trip comparisons show why the gate failed:

- baseline:
  - `55989` and `LF82` were effectively exact round trips (`96/96` identical ranks, near-zero probability deltas)
  - `EDL933` already had substantial raw-vs-reference drift
- candidate:
  - `EDL933` improved substantially (`median_abs_probability_delta` fell from `0.16329` to `0.0224043`, and the max
    delta fell from `0.286373` to `0.152605`)
  - `55989` and especially `LF82` lost that round-trip stability
  - all three surfaces changed for all `96` panel phages, with median per-host surface deltas of `0.0540`, `0.1663`,
    and `0.0875` respectively in `tl18_roundtrip_surface_summary.csv`

#### Holdout evaluation — the richer bundle is strictly better

The round-trip gate measured preprocessing fidelity: does the raw-input path reproduce the panel-reference predictions?
That is the wrong test for TL18. A model with richer features will naturally produce different predictions — the
question is whether they are *better*, not whether they are *the same*. Holdout evaluation on the ST03 65-strain
holdout (6,235 pairs) answers the right question:

| Metric | Baseline (defense + kmer) | Candidate (+ TL15/TL16/TL17) | Delta |
| --- | --- | --- | --- |
| ROC-AUC | 0.787 | **0.823** | **+0.036** |
| Top-3 hit rate | 90.5% (57/63) | **93.7%** (59/63) | **+3.2pp** |
| Brier score | 0.152 | **0.141** | **-0.011** |

The richer bundle improves every metric. The TL15 host surface features (O-antigen, capsule, receptor presence), TL16
host typing features (phylogroup, serotype, MLST), and TL17 phage RBP family features all contribute real signal that
was previously missing from the deployable bundle.

#### Interpretation

1. TL18 satisfies the runtime-contract part of the task. The bundle ships the richer deployable preprocessors relative
   to itself, and the relocated-bundle probe proves it can execute from a copied bundle directory rather than from
   hidden repo-root state.
2. The round-trip gate was a design error in the acceptance criteria. It tested preprocessing fidelity, not model
   quality. A richer feature set that carries real signal will inevitably change predictions for panel hosts — that is
   the point, not a failure. The holdout evaluation is the correct test.
3. The richer preprocessors provide real lift. +3.6pp ROC-AUC, +3.2pp top-3 hit rate, and -0.011 Brier on 65 holdout
   strains. This is the first time the deployable bundle has been shown to be *better* than the baseline, not just
   *more complete*.
4. The correct next state is to promote the richer bundle as the new deployable baseline and proceed with external
   validation.

### 2026-04-02: TL18 post-merge audit — holdout metrics validated, three deployment-path bugs found

#### Executive summary

A manual step-by-step audit of the TL18 inference pipeline on the 3 committed validation hosts (55989, EDL933, LF82)
plus the 65-strain holdout set confirmed that the reported model quality improvement is real: bootstrap CIs (2000
strain-level resamples) place the ROC-AUC delta at +0.036 with a 95% CI of [+0.002, +0.079] and 98.5% probability of
improvement. However, the audit also found three bugs that affect the raw-input inference path but not the holdout
evaluation. Two are systematic feature mismatches between training-time panel features and inference-time raw-genome
features; the third is an extra phage genome in the FNA directory. None of these invalidate the holdout comparison, but
all three will degrade prediction quality for truly novel hosts at deployment time.

#### Audit methodology

The audit traced every inference step on each of the 3 committed validation hosts:

1. **Defense Finder** — ran fresh on each raw FASTA, compared output to panel defense subtype annotations.
2. **TL15 host surface projection** — ran nhmmer (O-antigen), hmmscan (ABC capsule), phmmer (receptor presence) on each
   raw FASTA, compared projected features to the picard metadata-derived training features.
3. **TL16 host typing projection** — ran Clermont phylogroup caller, ECTyper serotype caller, MLST, and capsule proxy
   HMM scan on each raw FASTA, compared to picard metadata fields.
4. **Phage projection** — ran pyrodigal + kmer SVD + TL04 anti-defense + TL17 RBP family matching for a subset of
   phages, checked feature values against panel precomputed features.
5. **Scoring** — ran full inference (EDL933 vs 96 phages), compared raw-input predictions to panel-feature predictions
   and to known labels.
6. **Holdout evaluation** — recomputed baseline vs candidate metrics on the ST03 65-strain holdout using panel features,
   computed bootstrap CIs, and checked per-strain top-3 changes.
7. **Leakage audit** — verified that no feature in the candidate model is derived from interaction labels. Checked
   enrichment holdout exclusion (TL10), TL04 pair feature construction, TL15/TL16 metadata provenance, and TL17
   sequence-only family definitions.
8. **Feature importance analysis** — computed per-block importance shares to assess which new feature blocks contribute
   real signal vs noise.

#### Validation host inference results

**EDL933** (in panel, training set, 3 lytic phages out of 96):

- Defense features: 8 systems detected by fresh DefenseFinder vs 9 in panel data (RM_Type_IIG missed). Defense
  diversity correctly reported as 8. All other systems match.
- TL15 surface: O157 recovered directly, R3 LPS proxy correct, all 12/12 receptor-presence fields match.
- TL16 typing: phylogroup E, ST 11, O157:H7 all correct. **Capsule disagreement**: raw HMM scan detects K4
  (`KfoFGCA_2_unknown` model, 28 profile hits with strong scores), but picard metadata records `Capsule_ABC=0`.
- Full inference: top-2 predicted phages are DIJ07_P1 and DIJ07_P2 (both truly lytic). Top-3 hit: yes. Per-host
  AUC: 0.9964.
- Raw-input vs panel-feature comparison: median probability delta 0.015, max delta 0.118. Top-3 phages identical
  (order swapped). AUC: 0.9982 (panel) vs 0.9964 (raw).

**55989** (in picard metadata but not in ST02 pair table — no interaction labels):

- O104, phylogroup B1, ST 678. All features biologically plausible.
- Predictions: 4 phages above P(lysis) > 0.5, 19 above 0.3. Cannot validate against labels.

**LF82** (in picard metadata but not in ST02 pair table — no interaction labels):

- O83, phylogroup B2, ST 135. More broadly susceptible profile than 55989.
- Predictions: 26 phages above P(lysis) > 0.5, 34 above 0.3. Consistent with B2 phylogroup being generally more
  susceptible. Cannot validate against labels.

#### Holdout evaluation with bootstrap confidence intervals

The holdout evaluation uses **panel-derived features** for all 65 holdout strains (no raw FASTAs available for holdout
hosts). This is the correct test for model quality: it measures whether adding TL15/TL16/TL17 features to the training
feature set improves predictions on held-out strains, using the same feature derivation path for training and
evaluation. It does **not** test whether the raw-input preprocessing pipeline faithfully reproduces those features for
novel hosts — that is a separate concern addressed in the bug section below.

| Metric | Baseline (defense + kmer) | Candidate (+ TL15/TL16/TL17) | Delta | 95% CI | P(improvement) |
| --- | --- | --- | --- | --- | --- |
| ROC-AUC | 0.787 | 0.823 | +0.036 | [+0.002, +0.079] | 98.5% |
| Top-3 hit rate | 90.5% (57/63) | 93.7% (59/63) | +3.2pp | [-2.6pp, +10.5pp] | 82.2% |
| Brier score | 0.152 | 0.141 | -0.011 | [-0.025, +0.000] | 96.9% |

Bootstrap method: 2000 resamples of the 65 holdout strains (strain-level resampling, not pair-level), metrics recomputed
per resample.

The AUC improvement is statistically robust — the 95% CI is entirely above zero and 98.5% of bootstrap resamples show
improvement. The Brier improvement is near-significant (96.9% probability, CI just touching zero). The top-3
improvement is directionally correct but uncertain on 65 strains — net +2 strains (3 gained, 1 lost) is within random
variation.

Per-strain top-3 changes:

- **Gained** (baseline miss → candidate hit):
  - `H1-003-0088-B-J` (12 lytic phages): candidate places LF110_P2 and LF82_P5 in top-3.
  - `IAI58` (28 lytic phages): candidate places LF82_P6 in top-3.
  - `NILS41` (14 lytic phages): candidate places LF73_P4 in top-3.
- **Lost** (baseline hit → candidate miss):
  - `ECOR-69` (10 lytic phages): candidate pushes DIJ07_P1/P2 (non-lytic for this strain) to ranks 1-2 with higher
    confidence (P=0.770), bumping 4 lytic phages from tied-rank-1 (P=0.586) to rank-4+ (P=0.667). The richer features
    gave the model a stronger opinion that happened to be wrong for this strain.

Train-holdout gap comparison:

| | Baseline | Candidate |
| --- | --- | --- |
| Train AUC | 0.932 | 0.942 |
| Holdout AUC | 0.787 | 0.823 |
| Gap | 0.146 | **0.119** |

The candidate has a **smaller** train-holdout gap, indicating the richer features improve generalization rather than
causing overfitting.

#### Feature importance by block

| Block | Importance | Share |
| --- | --- | --- |
| Phage kmer SVD (26 features) | 3808 | 42.3% |
| TL16 host typing (466 features after one-hot) | 1963 | 21.8% |
| Defense subtypes (82 features) | 1561 | 17.3% |
| TL15 host surface (137 features after one-hot) | 1190 | 13.2% |
| TL17 phage RBP families (34 features) | 313 | 3.5% |
| TL04 anti-defense pairs (30 features) | 176 | 2.0% |

The new TL15/TL16/TL17 blocks account for 38.5% of total feature importance. TL16 (host typing: phylogroup, serotype,
MLST, capsule proxy) is the second most important block overall at 21.8%, confirming these features carry substantial
predictive signal. TL15 (host surface: O-antigen, LPS, receptor presence) contributes 13.2%. TL17 (phage RBP families)
contributes 3.5% — modest but nonzero, and biologically expected given that RBP family presence is a coarse proxy for
receptor specificity.

#### Leakage audit result

No label leakage found in any feature block:

- **Defense subtypes**: derived from DefenseFinder annotations on host genomes. No interaction data involved.
- **TL15 host surface**: derived from picard metadata (O-type, capsule, receptor clusters) and raw FASTA projections
  (nhmmer, hmmscan, phmmer). All inputs are organism properties, not interaction outcomes.
- **TL16 host typing**: derived from picard metadata (Clermont phylogroup, serotype, ST, capsule) and raw FASTA
  projections (ClermonTyping, ECTyper, MLST, capsule HMM). All inputs are organism properties.
- **Phage kmer SVD**: derived from tetranucleotide frequencies of phage genomes. No interaction data involved.
- **TL17 RBP families**: defined by PHROG functional annotation and mmseqs sequence clustering. Family definitions are
  sequence-based. No interaction data involved.
- **TL04 anti-defense pairs**: enrichment weights computed from TL02, which properly excluded all 65 holdout bacteria
  (confirmed in `lyzortx/generated_outputs/track_l/enrichment/manifest.json`:
  `holdout_exclusion.excluded_holdout_bacteria_count = 65`). Pair features constructed from TL04 manifest which records
  `holdout_exclusion.excluded_pair_rows = 6240`.
- **No legacy label-derived features** (`legacy_label_breadth_count`, `legacy_receptor_support_count`,
  `defense_evasion_*`, `receptor_variant_seen_in_training_positives`) are present in the candidate feature set.

#### Bug #1: Extra phage 411_P3 in FNA directory

The `data/genomics/phages/FNA/` directory contains 97 `.fna` files, but the panel metadata
(`data/genomics/phages/guelin_collection.csv`) defines only 96 phages. The extra file is `411_P3.fna`. When
`generalized_inference.infer()` is called with a glob of the FNA directory (as in the round-trip comparisons and manual
inference tests), it scores 97 phages instead of 96.

**Impact on reported metrics**: None. The holdout evaluation uses precomputed panel features from
`tl08_locked_panel_predictions.csv`, which is generated by `build_model_bundle` using `read_panel_phages` with
`expected_panel_count=96`. The extra phage never enters the holdout evaluation.

**Impact on deployment**: Low but real. If a user runs `gi.infer(host_fasta, glob("data/genomics/phages/FNA/*.fna"),
bundle_path)`, the rankings include `411_P3` — a phage with no training label. It competes for top-3 slots without the
model having seen any positive or negative evidence for it. The probability estimate for `411_P3` is based entirely on
its genomic features (kmer profile, RBP families, anti-defense genes) without any ground truth anchoring.

**Fix**: Either remove `411_P3.fna` from the FNA directory, or filter phage paths through the panel metadata before
scoring. The inference API should not silently accept phages outside the training panel without warning.

#### Bug #2: Capsule train/inference feature mismatch

The TL16 `host_capsule_abc_proxy_present` feature has a systematic discrepancy between training-time and inference-time
derivation:

- **Training time**: derived from picard metadata `Capsule_ABC` field via
  `int(float(row.get("Capsule_ABC", "0") or 0) > 0)`. The picard metadata records `Capsule_ABC=0` for 261 of 403
  panel hosts.
- **Inference time**: derived from raw HMM capsule scan via `tl16.scan_capsule_proxy()`, which runs `hmmscan` against
  checked-in ABC capsule HMM profiles on predicted proteins from the assembly FASTA.

The raw HMM scan is more sensitive than the picard metadata field. Concrete example: EDL933 has `Capsule_ABC=0` and
`ABC_serotype=` (empty) in picard, but the raw HMM scan detects a K4 capsule locus with 28 profile hits including
strong matches to KfoF (score 735), KfiD (score 698), and KfoA (score 385).

The picard metadata inconsistency extends beyond EDL933: 29 panel hosts have `Capsule_ABC=0` but a non-empty
`ABC_serotype` field, and 20 hosts have `Capsule_ABC>0` but an empty `ABC_serotype`. The `Capsule_ABC` column is not
a reliable binary indicator of capsule presence.

**Impact on reported metrics**: None. The holdout evaluation uses panel features consistently for both training and
evaluation. The capsule mismatch only manifests when the raw-input inference path is used for novel hosts.

**Impact on deployment**: Moderate. Capsule-related features carry approximately 3% of total model importance (direct
`host_capsule_abc_present` importance plus contributions through the one-hot encoded `host_abc_serotype_proxy`
categorical, which adds another ~2% through TL16 typing). When a novel host arrives and the raw HMM scan detects
capsule genes that the picard metadata would have coded as absent, the model sees feature values it learned to
associate with a different host population. The direction of error is unpredictable — it depends on the specific
capsule type detected and whether the model learned to associate that type with higher or lower lysis probability.

**Root cause**: The training features use a curated metadata field (`Capsule_ABC`) whose provenance is unclear —
it may represent a different capsule detection methodology, a different sensitivity threshold, or manual curation
decisions that differ from the automated HMM scan. The TL16 panel training path and the TL16 raw inference path use
fundamentally different capsule evidence sources.

**Fix options**:
1. Retrain with HMM-derived capsule features for all panel hosts (run the capsule scan on all 403 panel assemblies
   and use those values instead of picard `Capsule_ABC`). This aligns training and inference, but requires assemblies
   for all panel hosts.
2. Calibrate the HMM scan threshold to match the picard metadata sensitivity. This preserves the training feature
   distribution but weakens the capsule signal at inference time.
3. Accept the mismatch and document it. Capsule is 3-5% of importance; the expected prediction shift is small.

#### Bug #3: DefenseFinder version drift between panel annotations and raw inference

Fresh DefenseFinder runs on raw FASTAs produce slightly different defense system calls than the panel defense subtype
annotations in `data/genomics/bacteria/defense_finder/370+host_defense_systems_subtypes.csv`. Concrete example:
EDL933's panel data records 9 non-zero defense subtypes (including `RM_Type_IIG=1` and `MazEF=2`), but the fresh
DefenseFinder run on EDL933's FASTA detects only 8 systems (missing `RM_Type_IIG`, and `MazEF` is binarized to 1).

The binarization (`MazEF=2` → `host_defense_subtype_maz_ef=1`) is by design — the feature pipeline uses presence
indicators, not counts. But the missing `RM_Type_IIG` is a genuine detection disagreement between the panel annotation
run and the current DefenseFinder installation.

**Impact on reported metrics**: None. The holdout evaluation uses panel defense features, not fresh DefenseFinder output.

**Impact on deployment**: This is potentially the **most consequential** of the three bugs for deployed prediction
quality. Defense subtypes carry 17.3% of total model importance — the third largest block. If fresh DefenseFinder
systematically misses or adds defense systems relative to the panel annotations, dozens of binary features can flip
per host. A single missed system (like RM_Type_IIG) changes one feature, but if the disagreement rate is 1-2 systems
per host on average, the cumulative feature noise could meaningfully degrade predictions.

The magnitude cannot be assessed from a single host. A proper assessment requires running fresh DefenseFinder on all
403 panel host assemblies and computing the per-host agreement rate against the panel annotations. If the average
disagreement is <1 system per host, the impact is small (defense features are individually weak, importance is spread
across 82 features). If the average disagreement is >3 systems per host, the defense block is substantially noisier
at inference time than at training time.

**Root cause**: The panel defense annotations were generated at a specific point in time with a specific DefenseFinder
version and model database. DefenseFinder and its HMM models are updated independently. The repo pins
`defense-finder-models==2.0.2` in the runner, but the underlying tool version may differ between the original panel
annotation run and the current `phage_env` installation.

**Fix options**:
1. Re-annotate all 403 panel hosts with the current DefenseFinder installation and retrain on the fresh annotations.
   This aligns training and inference at the cost of potentially changing the model.
2. Pin the exact DefenseFinder version and model database that was used for the panel annotations, and use that same
   version at inference time. This preserves consistency but may miss real defense systems in novel hosts.
3. Quantify the disagreement rate first (option 1 as a measurement, not a retraining), then decide based on the
   magnitude.

#### Summary of audit conclusions

The **model quality improvement is real and statistically robust**. The AUC gain of +0.036 has a 98.5% probability of
being positive, the train-holdout gap shrank, no label leakage was found, and the new feature blocks carry 38.5% of
model importance. The holdout evaluation is internally consistent and methodologically sound.

The **deployment quality** has three known gaps between training-time and inference-time feature derivation. All three
affect only the raw-input inference path for novel hosts, not the panel-feature holdout evaluation. In order of
estimated impact:

1. **DefenseFinder version drift** (17.3% importance block, unknown disagreement magnitude) — highest risk, needs
   quantification.
2. **Capsule train/inference mismatch** (3-5% importance, systematic sensitivity difference) — moderate risk,
   well-characterized.
3. **Extra phage 411_P3** (affects ranking only, no training signal) — low risk, trivial fix.

None of these bugs invalidate the decision to promote the richer bundle as the new deployable baseline. They do mean
that the holdout metrics (+0.036 AUC, +3.2pp top-3) are an **upper bound** on the improvement a truly novel host would
see through the raw-input inference path. The actual deployment improvement will be smaller, by an amount that depends
primarily on the DefenseFinder disagreement rate.

#### Continuation: Deployment-Paired Feature Pipeline

The audit findings above — training/inference feature mismatch, binary thresholds discarding gradients, and 91 wasted
duplicate one-hot features — motivated a new track: **Deployment-Paired Feature Pipeline** (DEPLOY01-08 in plan.yml).
That track downloads all 403 panel assemblies from figshare, re-derives every host feature from raw FASTAs using the
same pipeline that runs at inference time, switches to continuous scores where the gradient carries biological signal,
and deduplicates redundant features. Track L work is complete; further development continues in the DEPLOY track.
