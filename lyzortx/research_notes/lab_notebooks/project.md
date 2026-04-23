### 2026-02-15: PLAN update after external-data literature deep dive

#### Summary

We updated `PLAN.md` to make Track I execution source-prioritized and evaluation source-aware. The main principle
remains unchanged: internal paper data is the baseline foundation, and external data is integrated incrementally as a
measured enhancer.

#### What changed in the plan

1. Track I now has explicit ingestion order for supervised external datasets: `VHRdb -> BASEL -> KlebPhaCol -> GPB`.
2. Track I now separates Tier A supervised matrices from Tier B weak-label sources (`Virus-Host DB`,
   `NCBI Virus/BioSample`).
3. Track I now requires a `source_registry.csv` capturing source metadata, label type, host resolution, assay type,
   license, and access path.
4. Track F now includes source-aware evaluation and leakage controls: leave-one-datasource-out and cross-source transfer
   checks.
5. Track G now includes explicit training sequence: internal-only baseline first, then Tier A, then optional Tier B.
6. Track J now includes external data-license and use-restriction tracking in manifests.

#### Why this improves execution quality

1. Prevents mixing noisy and high-quality sources too early.
2. Makes performance gains attributable to specific sources.
3. Reduces hidden leakage risk across merged external datasets.
4. Preserves reproducibility and compliance when datasets have different use terms.

### 2026-02-15: Defining a Meaningful Model and Product-Oriented Benchmarks

Based on a review of the project plan and the available data, we've established a strategic framework for what
constitutes a "meaningful model" and how to evaluate its utility for a clinical product that recommends phage cocktails
for new _E. coli_ strains.

#### 1. Data Sufficiency Assessment

We have a sufficient foundation of raw data (`raw_interactions.csv`, genomic assemblies) to build a predictive model.
However, success is contingent on executing the data integrity and feature engineering tracks outlined in `PLAN.md`
(Tracks A, C, D, E). The raw data is not yet model-ready, and a significant effort is required to create clean labels
and potent predictive features.

#### 2. Defining a Meaningful Model

A "meaningful model" in this context is not a simple classifier. It is a system that produces a **calibrated probability
of lysis** for any given phage-bacterium pair.

- **Input:** The genomic assembly of a new _E. coli_ strain.
- **Output:** A ranked list of all phages in our panel, sorted by their predicted probability of successfully lysing the
  new strain, e.g., `P(Lysis | new_strain, phage_i)`.

This probabilistic output is essential for the downstream cocktail recommendation algorithm (Track H), which needs to
weigh evidence and confidence, not just binary predictions.

#### 3. Proposed Product-Oriented Evaluation Benchmarks

To be considered "very useful" for a clinical product, the model must meet a set of rigorous benchmarks that go beyond
standard academic metrics.

- **Level 1: Foundational Model Performance (Offline)**

  - **Top-3 Lytic Hit Rate > 95%:** The 3-phage cocktail recommended by the model must contain at least one truly lytic
    phage for over 95% of new strains.
  - **Precision @ P(Lysis) > 0.9 must be > 99%:** When the model is highly confident, it must be extremely reliable to
    build clinical trust.

- **Level 2: Cocktail Recommendation Utility (Simulation)**

  - **Simulated 3-Phage Cocktail Coverage > 98%:** An end-to-end simulation must show that the final recommended
    cocktails are effective against over 98% of strains in the held-out test set.

- **Level 3: Safety and Robustness (Negative Controls)**
  - **Sentinel Strain Recovery = 100%:** The model must correctly identify known effective phages for a predefined set
    of challenging "sentinel" strains, ensuring it doesn't fail on critical edge cases.

#### 4. Rationale for 3-Phage Cocktails

The choice of a 3-phage cocktail is a deliberate trade-off between efficacy and complexity.

- **Why Not Fewer (1-2 Phages)?**

  - **Higher Risk of Resistance:** Bacteria can rapidly evolve resistance to a single phage. Using multiple phages
    targeting different receptors or using different lytic mechanisms makes this much harder.
  - **Lower Coverage:** A single phage is less likely to be effective against a new, unknown strain.

- **Why Not More (4+ Phages)?**
  - **Manufacturing & Regulatory Complexity:** Each phage in a cocktail must be individually produced, characterized,
    and proven safe and stable. The regulatory burden and cost increase exponentially with each added phage.
  - **Risk of Antagonistic Interference:** Phages can compete for host cell resources or even block each other's entry,
    reducing the overall efficacy of the cocktail.
  - **Diminishing Returns:** The marginal benefit of adding a fourth or fifth phage is often small and does not justify
    the significant increase in complexity and cost.

#### 3-phage cocktail rationale note

A 3-phage cocktail is considered the current industry "sweet spot," providing robust defense against bacterial
resistance while remaining tractable from a manufacturing and regulatory perspective.

### 2026-02-15: Research on External Phage-Host Interaction Datasets

#### External Data Summary

To improve model performance, we must expand our training data beyond the internal `raw_interactions.csv`. Research was
conducted to identify public databases that contain phage-host interaction data. The findings indicate that while a
direct equivalent to our raw interaction matrix is rare, a significant volume of valuable data can be compiled from
various sources.

#### Key Data Sources Identified

1. **Dedicated Phage-Host Databases:** These are the highest-value targets for finding experimentally confirmed
   interaction pairs.

   - **Relevant Databases:** `Virus-Host DB`, `ViralHostRangeDB`.
   - **Content:** Aggregate interaction data from literature, providing "Phage X infects Bacterium Y" pairs. They
     typically lack the fine-grained dilution/replicate data but are excellent for expanding our core training set.

2. **Public Sequence Archives (NCBI, CNCB):** These contain the raw data from sequencing experiments and associated
   metadata.

   - **NCBI BioSample Database:** A good source for "known positive" interactions. The `isolation_host` field in a
     phage's BioSample record provides a confirmed host. This is a low-effort way to augment our dataset with thousands
     of positive examples.
   - **NCBI Sequence Read Archive (SRA) / China National Center for Bioinformation (CNCB) Genome Warehouse (GWH):** A
     potential goldmine, but high-effort. These archives contain raw data from high-throughput screening projects.
     Extracting an interaction matrix would require a significant bioinformatics effort (re-processing raw reads,
     mapping metadata), making it a long-term strategic goal.

3. **International Databases (Focus on China):**
   - **China National Center for Bioinformation (CNCB):** This is a key resource, hosting databases like the **Virus &
     Host (V&H) database**, which aggregates data from Chinese research and is a priority for exploration.

#### Next Steps & Strategy

The most pragmatic approach is to start with the lowest-effort, highest-value tasks.

1. **Short-Term:** Systematically query the APIs of **Virus-Host DB** and the **NCBI BioSample database**. The goal is
   to script a process to download all available phage-host pairs and integrate them as "known positive" in our dataset.
2. **Mid-Term:** Investigate the curated databases at the CNCB.
3. **Long-Term:** Scope the effort required to re-process a full screening project from the SRA or GWH.

### 2026-03-21: Strategic plan revision — v1 push

#### Context

Steel thread v0 is complete and all go/no-go gates pass. The v0 logistic regression on metadata features achieves 84.6%
top-3 hit rate (AUC 0.827) on the holdout set. However, the model has a "popular phage" bias: it recommends the same
broad-range Myoviridae for almost every strain because it has zero compatibility signal between specific phage RBPs and
specific host receptors. The strict-confidence slice drops to 62.5%.

Meanwhile, the repo contains rich genomic data that is completely unused: 138 defense system subtypes
(`370+host_defense_systems_subtypes.csv`), 12 outer-membrane receptor variant clusters
(`blast_results_cured_clusters=99_wide.tsv`), per-phage RBP annotations (`RBP_list.csv`), 97 complete phage genomes
(`FNA/`), UMAP phylogenomic embeddings (`coli_umap_8_dims.tsv`), and a pangenome matrix (`unique_host_genes.csv`).

This revision refocuses the plan to maximize prediction accuracy by exploiting this untapped data, targeting mid-May 2026
as an aspirational deadline for discussion-ready results.

#### What changed in the plan

**Tracks kept as-is (complete):** ST (11/11), A (10/10).

**Tracks cut or radically downsized:**

- **B** (EDA): marked done. TB06 (uncertainty map) and TB07 (mechanistic hypotheses) cut — TB07 is subsumed by actually
  building features in C/D/E; TB06 is nice-to-have, not blocking.
- **F** (Splits/Eval): 12 tasks → 2 tasks. ST03 already provides leakage-safe host-group and phage-family holdouts. Keep
  only: lock existing split as v1 benchmark + bootstrap CIs.
- **G** (Modeling): 14 tasks → 4 tasks. Two-stage mechanistic decomposition (P(adsorption) × P(lysis|adsorption))
  requires labeled adsorption outcomes that don't exist. Multi-task learning for strength/potency is lower-ROI than
  fixing the core binary prediction. Keep: LightGBM model + calibration + ablation suite + SHAP explanations.
- **H** (Cocktail): 8 tasks → 2 tasks. The heuristic recommender works at 84.6%. Optimization-based cocktail design is a
  later concern. Keep: existing top-3 + explained recommendations with SHAP features.
- **J** (Reproducibility): 7 tasks → 2 tasks. Keep: one-command regeneration + environment freeze.
- **K** (Wet-Lab): eliminated. No wet-lab access exists. Held-out strain evaluation serves the same credibility purpose.

**Tracks refocused (critical path):**

- **C** (Host Features): 6 tasks → 4 tasks. Defense subtypes (138 binary cols from defense_finder), OMP receptor variants
  (12 proteins with cluster IDs), capsule/LPS detail, and UMAP embeddings.
- **D** (Phage Features): 7 tasks → 3 tasks. RBP features from `RBP_list.csv`, k-mer embeddings from 97 FNA genomes
  (tetranucleotide SVD), phage distance embedding from VIRIDIC tree.
- **E** (Pairwise Features): 4 tasks → 3 tasks, moved to stage 1. RBP×receptor compatibility, defense evasion proxy
  (collaborative filtering from training data), phylogenetic distance to isolation host.

**Tracks kept at full scope:** I (External Data) — full Tier A pipeline: VHRdb + BASEL + KlebPhaCol + GPB ingestion with
strict ablation sequence.

**New track added:** P (Presentation) — 3 tasks: digital phagogram visualization, panel coverage heatmap, feature lift
visualization.

**Net effect:** 101 tasks → ~37 tasks. The cut tasks are done, deferred, or eliminated as low-ROI.

#### Execution timeline (aspirational, ~8 weeks)

- **Weeks 1–3 (Phase 1):** Feature engineering sprint. Expand pair table from ~28 metadata features to ~160–200 genomic
  features across C, D, E.
- **Weeks 2–6 (Parallel):** External data integration (Track I). Full Tier A ingestion with ablations.
- **Weeks 3–5 (Phase 2):** Model upgrade. LightGBM replacing logistic regression, calibration, ablation suite, SHAP.
- **Weeks 5–7 (Phase 3):** Evaluation and presentation artifacts. Bootstrap CIs, before/after comparison, explained
  recommendations, visualizations.
- **Week 8 (Phase 4):** Buffer. One-command reproducibility, environment freeze.

#### Expected performance targets

| Metric | v0 (current) | v1 (target) | Source of lift |
|--------|-------------|-------------|----------------|
| Top-3 hit rate (full-label) | 84.6% | 90–93% | OMP receptor + RBP compatibility resolves "popular phage" bias |
| Top-3 hit rate (strict-conf) | 62.5% | 72–78% | Defense subtypes provide discriminative signal |
| AUC | 0.827 | 0.87–0.90 | GBM captures nonlinear defense×phage interactions |
| Brier score | 0.171 | 0.12–0.15 | Better feature set + GBM calibration |

#### Risk factors

1. RBP data is sparse — not all phages have annotations. Handle with indicator features.
2. Defense subtype sparsity — many subtypes in <5 strains. Aggressive variance filtering needed.
3. E2 (defense evasion proxy) leakage risk — must compute strictly on training fold per CV split.
4. RBP-receptor lookup curation — requires ~2–3 days of manual literature work.
5. GBM overfitting — with ~200 features and 29K training pairs, need careful regularization.
6. Pangenome data deferred — `unique_host_genes.csv` (7,511 records) parked unless defense+receptor features plateau.

### 2026-03-22: v1 model results — ablation paradox and feature-selection plan adjustment

#### What happened

Tracks C, D, E, and most of G landed in a single sprint. The v1 LightGBM model on 191 features (21 categorical + 170
numeric) was trained, calibrated, and ablated against the locked ST03 holdout (65 strains, 6,235 pairs).

#### Headline v0 → v1 comparison

| Metric | v0 (logreg, metadata) | v1 (LightGBM, all features) | Delta |
|--------|----------------------|---------------------------|-------|
| ROC-AUC | 0.827 | **0.910** | +0.083 |
| Top-3 hit rate (all strains) | 84.6% | **89.2%** | +4.6% |
| Top-3 hit rate (susceptible only) | 87.3% | **92.1%** | +4.8% |
| Brier score | 0.171 | **0.113** | -0.058 |
| ECE (isotonic, full-label) | 0.032 | **0.020** | -0.012 |
| ECE (isotonic, strict-conf) | 0.124 | **0.094** | -0.030 |

AUC exceeded the 0.87–0.90 target at 0.910. Calibration is excellent (ECE 0.020 isotonic on full-label). Top-3 hit rate
at 89.2% is close to the 90% target but does not clear it.

#### The ablation paradox

The TG03 ablation suite revealed an unexpected pattern: individual feature blocks outperform the combined model on the
top-3 ranking metric.

| Arm (v0 baseline + one block) | Top-3 hit rate | ROC-AUC | Brier |
|-------------------------------|---------------|---------|-------|
| v0 only (metadata baseline) | 86.2% | 0.908 | 0.114 |
| **+defense subtypes** | **90.8%** | 0.907 | 0.114 |
| +OMP receptors | 87.7% | **0.910** | **0.112** |
| **+phage genomic** | **90.8%** | 0.909 | 0.112 |
| +pairwise compatibility | 87.7% | 0.905 | 0.117 |
| all features combined | 87.7% | 0.909 | 0.113 |

Key observations:

1. **Defense subtypes and phage genomic features are the clear winners.** Each independently pushes top-3 hit rate to
   90.8% (+3 holdout strains recovered). These are the features that break the "popular phage" bias — they encode *which
   specific defense systems* a strain carries and *which specific phage genome architecture* can evade them.

2. **OMP receptors win on discrimination but not ranking.** Best AUC (0.910) and Brier (0.112) as a single block, but
   only 87.7% top-3 hit rate. The receptor variants help separate lytic from non-lytic pairs more precisely, but that
   precision does not translate into moving the *right* phages into the top-3 slots.

3. **Pairwise compatibility features hurt.** Adding them alone degrades both AUC (0.905, worst of all arms) and Brier
   (0.117, worst of all arms). The genus-level receptor lookup covers 80% of phages, but the compatibility signal may
   be too coarse, or the defense evasion proxy may be partially redundant with the raw defense subtype block.

4. **The all-features model does not dominate.** At 87.7% top-3 / 0.909 AUC, it underperforms the single-block defense
   and phage-genomic arms on ranking. This means feature interactions are creating noise when all blocks are thrown
   together without selection.

#### Why this happens

The most likely explanation: when all 191 features are present, the GBM splits on pairwise-compatibility and
OMP-receptor features that have good binary discrimination (high AUC) but that *dilute the ranking signal* from defense
and phage-genomic features. Top-3 hit rate is a ranking metric sensitive to the relative ordering of the top few phages
per strain, not just the classification boundary. A feature that improves average-case AUC can still hurt the top-k
ranking for specific strains by pushing a marginally-higher-scoring wrong phage above a correct one.

#### What this means for the plan

The 90.8% top-3 hit rate from defense-only or phage-genomic-only arms proves the target is reachable — we just need
to find the right feature combination that preserves it.

**Plan adjustment (two changes):**

1. **Added TG05: feature-subset sweep.** Train models on all 2-block and 3-block combinations of the 4 new feature
   blocks (defense, OMP, phage-genomic, pairwise). Identify the winning subset that maximizes top-3 hit rate without
   degrading AUC. Lock the final v1 feature configuration for downstream Tracks F, H, and P.

2. **Extended TG04 acceptance criteria.** SHAP analysis must now also produce a concrete recommendation of which feature
   blocks to keep in the final v1 model, informed by both SHAP evidence and TG03 ablation results. This ensures SHAP
   is not just descriptive but prescriptive.

No other plan changes needed. The pairwise compatibility work (Track E) was still worth doing — the features may become
useful after refinement (e.g., finer-grained receptor lookup, per-genus rather than per-family evasion rates) or in the
optimization-based recommender (Track H). We just should not force them into the v1 model if they hurt ranking.

#### Next steps (priority order)

1. **TG04** (SHAP): understand *which* pairwise features cause the ranking degradation and whether any are worth keeping.
2. **TG05** (feature-subset sweep): find the combination that clears 90%+ top-3 on holdout.
3. **TF01** (bootstrap CIs): run on the winning feature configuration, not necessarily the all-features model.
4. **TF02** (before/after error analysis): compare the winning v1 model against v0 on the specific holdout miss strains.

### 2026-03-22: Per-task model selection for Codex CI — cost analysis and assignments

#### Problem

The orchestrator dispatches all tasks to `gpt-5.4` via `codex-implement.yml`, regardless of task complexity. After the
first full day of automated v1 sprint execution (8 implementation runs for TD03, TE01–TE03, TG01–TG03, and TG04
attempt), total token consumption reached **856,463 tokens** with an average of **106,752 tokens per successful run**.

At gpt-5.4 pricing ($2.50/1M input, $15.00/1M output), this represents significant cost. Meanwhile, gpt-5.4-mini
($0.75/1M input, $4.50/1M output — approximately **70% cheaper**) scores 54.4% on SWE-Bench Pro and handles
straightforward coding tasks well. The key trade-off: gpt-5.4-mini has a 400K context window (vs 1.05M for gpt-5.4)
and shows a 5–10% accuracy gap on complex reasoning tasks.

#### Token usage breakdown from today's runs

| Task | Tokens | Status | Complexity |
|------|--------|--------|------------|
| TD03 (phage distance embedding) | 58,488 | success | Low — single MDS embedding from Newick tree |
| TE01 (RBP-receptor compatibility) | 134,626 | success | Medium — curated lookup + pair feature generation |
| TE02 (defense evasion proxy) | 103,687 | success | Medium — collaborative filtering with leakage guard |
| TE03 (phylogenetic distance) | 114,781 | success | Low-Medium — UMAP distance computation |
| TG01 (LightGBM training) | 129,532 | success | High — hyperparameter tuning, CV, comparator model |
| TG02 (calibration) | 85,193 | success | Medium — isotonic/Platt scaling, metric reporting |
| TG03 (ablation suite) | 120,957 | success | Medium — parameterized training loop, metric collection |
| TG04 (SHAP explanations) | 109,199 | failure | High — TreeExplainer, interpretation, cross-referencing |

Observations: TD03 (simplest task) used 58K tokens. TG01 and TE01 (most complex) used 130K+ tokens. The failed TG04
attempt still consumed 109K tokens. There is no correlation between token usage and task complexity strong enough to
predict model needs from token count alone — the decision must be based on task characteristics.

#### Waste patterns observed

Analysis of shell command patterns across today's 8 runs revealed systematic waste:

- **12 environment discovery issues**: Codex attempts `micromamba activate` in CI where micromamba does not exist. This
  is partially addressed in AGENTS.md but the agent still fumbles on first attempts.
- **30 git config attempts**: Git identity is pre-configured by the workflow, but Codex tries to set it again.
- **24 failed commands**: Various command failures across runs (exit codes 1, 127, 128).

These waste patterns account for an estimated 5–10% of total tokens. Fixing them helps but does not change the
fundamental cost picture — the bulk of tokens are legitimate implementation work.

#### Model capability comparison

**gpt-5.4** (full model):
- 1,050,000 token context window
- Highest accuracy on complex reasoning, architectural decisions, multi-file coordination
- Best for tasks requiring deep domain knowledge, novel algorithmic design, or cascading design decisions
- $2.50/1M input, $15.00/1M output

**gpt-5.4-mini**:
- 400,000 token context window
- 54.4% on SWE-Bench Pro (vs gpt-5.4 full score)
- 2x+ faster inference
- Strong on targeted edits, standard library usage, parameterized loops, visualization code
- $0.75/1M input, $4.50/1M output (~70% cheaper)

The 400K context window is sufficient for all tasks in this repo — the largest single-step context needed is the pair
table (~30K pairs × ~200 features), which is loaded from disk, not from the prompt. The deciding factor is reasoning
quality, not context size.

#### Assignment methodology

Each of the 16 pending tasks was evaluated on four axes:

1. **Novelty**: Does the task require implementing a pattern not yet in the codebase, or does it follow established
   patterns (TG01 training loop, TG03 ablation structure, ST05 calibration)?
2. **Domain criticality**: Do design errors cascade to downstream tasks? A wrong harmonization protocol (TI05) poisons
   all external data work; a wrong bar chart color (TP03) is trivially fixable.
3. **Reasoning depth**: Does the task require multi-step logical reasoning (combinatorial search, leakage analysis,
   SHAP interpretation) or straightforward data assembly (bootstrap resampling, metric aggregation)?
4. **Established patterns**: Can the task largely reuse existing code structure? Tasks that follow TG01/TG03 patterns
   (train model, collect metrics, output CSV) are good candidates for gpt-5.4-mini.

#### Assignments: gpt-5.4 (5 tasks)

**TG04 — Compute SHAP explanations for per-pair and global feature importance**
- TreeExplainer integration is new to the codebase (no existing SHAP usage to copy)
- Must cross-reference TG03 ablation results with per-feature SHAP values to produce prescriptive recommendations
- Per-strain narrative synthesis ("what makes strain X hard to predict") requires domain interpretation
- Acceptance criteria explicitly require a concrete recommendation of which feature blocks to keep — this is a design
  decision, not a computation

**TG05 — Run feature-subset sweep to find best block combination for top-3 ranking**
- 10 combinatorial model runs (C(4,2) + C(4,3) = 6 + 4) with careful feature-block bookkeeping
- Must identify the winning subset, compare against the TG01 all-features model, and lock the final v1 configuration
- The "lock" decision affects all downstream tracks (F, H, P) — getting it wrong means rework
- Unlike TG03 (sequential ablation), this requires combinatorial logic and winner-selection heuristics

**TI05 — Define harmonization protocol for Tier A datasets**
- Multi-source schema alignment across VHRdb, BASEL, KlebPhaCol, and GPB — each has different column semantics, label
  types, and confidence levels
- Design decisions here cascade to TI06 (Tier B ingestion), TI07 (confidence tiers), TI08 (integration), and TI09
  (ablation sequence)
- Requires domain understanding of what "lysis" means in each source's assay context
- Errors are silent and expensive to detect (wrong label mapping looks correct until model evaluation)

**TI07 — Define confidence tiers for external labels**
- Subjective tier design: what confidence level does a VHRdb curated record deserve vs an NCBI BioSample isolation_host
  record?
- Weighting strategy for training (how much to down-weight low-confidence labels) requires balancing coverage vs noise
- Cascades to TI08 (integration uses tier weights) and TI09 (ablation sequence is ordered by tier)
- The ST01b strict-confidence concept provides a pattern, but extending it to external sources with heterogeneous assay
  types is a design challenge

**TI08 — Integrate external data as non-blocking enhancer**
- Architectural pattern: external data must be truly optional — the internal-only pipeline must remain runnable
- Leakage prevention: external datasets may contain organisms that overlap with holdout strains — careful ID validation
  needed
- Fallback handling: graceful degradation when external data files are absent (CI runs without them)
- The "non-blocking" constraint means conditional imports, feature-flag-like behavior, and defensive coding — harder to
  get right than it looks

#### Assignments: gpt-5.4-mini (11 tasks)

**TF01 — Lock ST03 split as v1 benchmark and add bootstrap CIs**
- Bootstrap resampling is a standard NumPy pattern (np.random.choice with replacement, 1000 iterations)
- Dual-slice filtering (full-label vs strict-confidence) already implemented in TG02
- Metric functions (AUC, top-3 hit rate, Brier, ECE) already exist in the codebase
- No design decisions — just apply existing patterns with confidence intervals

**TF02 — Before/after comparison of v0 vs v1 with error bucket analysis**
- Side-by-side metric table is a DataFrame join between ST05 and TG02 outputs
- Error bucket identification is algorithmic: for each holdout strain, did v0 miss and v1 hit (or vice versa)?
- The "honest reporting" requirement is met by listing strains that remain unpredictable — no narrative synthesis needed

**TH02 — Add explained recommendations with calibrated P(lysis), CI, and SHAP features**
- Top-3 recommendation assembly reuses ST06 logic with TG02 predictions
- SHAP feature extraction is a DataFrame merge with TG04 output (pair_id → top-3 SHAP features)
- CI computation is percentile-based from bootstrap or calibration data
- Output formatting for "clinician-ready" display is light presentation logic

**TI06 — Tier B weak-label ingestion (Virus-Host DB, NCBI BioSample)**
- ID cross-referencing against canonical maps from Track A (already built)
- Follows TI03/TI04 ingestion patterns (source fidelity preservation, standardized output format)
- Confidence tier assignment uses rules defined in TI07 (upstream dependency)

**TI09 — Run strict ablations in sequence**
- Parameterized loop over TG01 training with progressively more external data
- 6 training runs: internal-only → +VHRdb → +BASEL → +KlebPhaCol → +GPB → +Tier B
- Metric collection identical to TG03 ablation suite structure
- Leakage verification is the only subtle part, but the holdout split is already locked

**TI10 — Track incremental lift and failure modes by datasource and confidence tier**
- GroupBy aggregation of TI09 ablation results by source and tier
- Failure mode detection: which strains regressed when adding each source?
- Follows ST09 error analysis pattern (documented hypotheses per miss strain)

**TJ01 — One command to regenerate all v1 outputs from raw data**
- Orchestration script that calls existing `run_track_*.py` entry points in dependency order
- Dependency validation (check prerequisites before running each track)
- Follows `run_track_a.py` structural pattern

**TJ02 — Freeze environment specs and seeds for v1 benchmark run**
- Documentation task: inventory all `random_state` parameters, export `pip freeze`
- Verification: run pipeline twice with frozen seeds, compare outputs
- No new algorithms or ML logic

**TP01 — Build digital phagogram visualization for per-strain phage ranking**
- Plotly or Matplotlib visualization of ranked phage list with P(lysis) and confidence bands
- Data inputs are well-defined (TG02 calibrated predictions, TG04 SHAP values)
- Standard charting library usage — no novel algorithmic work

**TP02 — Build panel coverage heatmap across strain diversity**
- Seaborn/Matplotlib heatmap with host phylogroup × phage family axes
- Aggregation: mean P(lysis) per cell from TG02 predictions grouped by ST02 taxonomy
- Hard-to-lyse gaps are visually evident from low-probability cells

**TP03 — Build feature lift visualization from ablation suite results**
- Bar chart from TG03 ablation CSV (already generated and available)
- Standard Matplotlib bar chart with v0 baseline reference line
- Simplest visualization task — no data processing beyond reading the CSV

#### Projected cost impact

Assuming average token usage of ~107K per run (today's observed average):

- **All gpt-5.4 (status quo):** 16 tasks × 107K tokens × ~$17.50/1M blended ≈ **$30**
- **With model selection:** (5 × 107K × $17.50/1M) + (11 × 107K × $5.25/1M) ≈ $9.36 + $6.18 ≈ **$15.50**
- **Estimated savings:** ~48% on implementation runs

This does not include lifecycle (review feedback) runs, which use the same model and add 1–3 rounds per task. With
lifecycle runs, the savings scale proportionally.

#### Implementation approach

Two PRs to minimize risk:

1. **PR 1 (data-only):** Add `model:` field to all 16 pending tasks in `plan.yml`. The existing `load_plan()` function
   ignores unknown YAML fields, so this is a no-op until PR 2 lands. This allows model assignments to be reviewed and
   adjusted independently of the code changes.

2. **PR 2 (code):** Wire the model field through `plan_parser.py` → `orchestrator.py` → issue body → workflow YAML.
   The orchestrator emits a `<!-- model: gpt-5.4-mini -->` HTML comment in the issue body. The workflow extracts it
   with a small Python helper (`parse_model_directive.py`) and passes it to `openai/codex-action@v1`. No default — if
   the model directive is missing, the workflow fails with a clear error message.

#### Open questions for future refinement

1. **Should lifecycle (review feedback) runs use the same model as the original implementation?** Currently planned: yes,
   for consistency. But review feedback is often simpler than initial implementation — a cheaper model might suffice.
2. **Should we track per-task model cost to validate these assignments?** The `ci_token_usage.py` tool already reports
   per-run tokens. After a few tasks run with model selection, we can compare actual mini vs full token usage.
3. **When should a task be upgraded from mini to full?** If a gpt-5.4-mini run fails, the model assignment can be
   changed in plan.yml and the orchestrator re-dispatched. No code change needed.

### 2026-03-22: Label leakage concern — deployment-realistic evaluation needed

#### Problem identified

The TG04 SHAP analysis revealed that the two highest-impact features in the v1 model are derived from training labels,
not from genomic content:

| SHAP rank | Feature | Mean |SHAP| | Source |
|-----------|---------|------------|--------|
| 1 | `legacy_label_breadth_count` | 1.183 | Panel metadata: count of lytic phages per strain |
| 2 | `legacy_receptor_support_count` | 0.572 | TE01: training-fold lysis frequency per receptor variant |

These features are powerful within-panel predictors but are **unavailable for truly novel strains** in a deployment
scenario. A new clinical isolate arrives with a genome assembly — you can derive its defense systems, receptor variants,
LPS type, and phylogenomic embedding, but you cannot know `legacy_label_breadth_count` (you haven't screened it yet) or
`legacy_receptor_support_count` (its specific receptor variant may not appear in the training data).

The model leans heavily on these: `legacy_label_breadth_count` alone has 2x the SHAP impact of the next feature.
This means the holdout metrics (AUC 0.910, top-3 89.2%) are **optimistic for real-world deployment** because the
holdout strains are
from the same experimental panel and their receptor variants overlap with training strains.

#### What is genuinely novel genomic signal?

The SHAP ranking after the two panel-artifact features shows real genomic signal:

| SHAP rank | Feature | Mean |SHAP| | Deployment-available? |
|-----------|---------|------------|----------------------|
| 3 | `phage_gc_content` | 0.217 | Yes — from genome |
| 4 | `phage_genome_length_nt` | 0.203 | Yes — from genome |
| 5 | `host_lps_type=R1` | 0.159 | Yes — from genome assembly |
| 6 | `defense_evasion_mean_score` | 0.130 | Partially — rates from training, defense profile from genome |
| 7 | `isolation_host_defense_jaccard_distance` | 0.123 | Yes — computable from genome |

These features have 5–10x smaller SHAP values than the panel artifacts. The model *is* learning genomic signal, but the
panel-memorization features dominate.

#### Why this matters

If we present the 89.2% top-3 hit rate to CDMO partners as "what the model can do for your new clinical isolate," we
are overpromising. The honest number — what the model achieves using only features available at deployment time — will be
lower. That's the number partners need to trust.

This also explains the ablation paradox from the previous note: adding the pairwise compatibility block (which includes
`legacy_receptor_support_count`) shifts the model's attention toward collaborative-filtering signal that is
strong on average but misleading for specific holdout strains.

#### Plan adjustment

Added two acceptance criteria to TG05 (feature-subset sweep):

1. Include a **deployment-realistic arm** that excludes all features derived from training labels
   (`legacy_label_breadth_count`, `legacy_receptor_support_count`) to measure generalization to truly novel strains.
2. Report both **panel-evaluation** and **deployment-realistic** metrics for the winning configuration.

This gives us two numbers to present: the panel metric (what the model does on known strains with full context) and the
deployment metric (what a CDMO partner should expect for a new isolate). Both are valuable — the panel metric validates
the approach, the deployment metric sets honest expectations.

#### No other plan changes needed

The concern does not invalidate Tracks C, D, or E. The genomic features are the right features — they just need to be
evaluated separately from the panel-artifact features. TG05 is the right place to do this, and the task is already
pending.

### 2026-03-22: TG05 results — deployment-realistic model outperforms panel model on ranking

#### TG05 sweep results

The feature-subset sweep (PR #156) evaluated all 10 two-block and three-block combinations of the four new feature
blocks (defense, OMP, phage-genomic, pairwise) with fixed TG01 hyperparameters on the ST03 holdout (65 strains).

**Panel-evaluation results (all arms include v0 baseline + legacy_label_breadth_count):**

| Arm | AUC | Top-3 (all) | Top-3 (susceptible) | Brier |
|-----|-----|-------------|---------------------|-------|
| TG01 all-features reference | 0.9091 | 87.7% | 90.5% | 0.113 |
| defense + phage-genomic | 0.9082 | 89.2% | 92.1% | 0.111 |
| OMP + pairwise | 0.9066 | 89.2% | 92.1% | 0.116 |
| **defense + OMP + phage-genomic (WINNER)** | **0.9108** | 87.7% | 90.5% | **0.110** |
| defense + phage-genomic + pairwise | 0.9082 | 84.6% | 87.3% | 0.113 |

Winner: **defense + OMP + phage-genomic** — best AUC (0.9108), best Brier (0.110), pairwise block excluded. The winner
selection rule required AUC ≥ TG01 all-features (0.9091), then maximized top-3. Only this arm cleared the AUC gate.

**Deployment-realistic result (winner minus legacy_label_breadth_count):**

| | AUC | Top-3 (all) | Top-3 (susceptible) | Brier |
|--|-----|-------------|---------------------|-------|
| Panel model | 0.911 | 87.7% | 90.5% | 0.110 |
| Deployment-realistic | 0.835 | **92.3%** | **95.2%** | 0.158 |

#### Interpretation

1. **Removing `legacy_label_breadth_count` improves ranking.** The deployment-realistic model achieves 92.3% top-3
   hit rate — higher than any panel-evaluation arm. This is counterintuitive but mechanistically sound:
   `legacy_label_breadth_count` tells
   the model "this strain is broadly susceptible" which biases it toward recommending the same popular broad-range
   phages. Without that shortcut, the model is forced to use defense subtypes, OMP receptor variants, and phage k-mer
   profiles to make strain-specific picks. Those features produce better *rankings* even though they produce worse
   *pairwise discrimination* (AUC drops from 0.911 to 0.835).

2. **AUC and top-3 measure fundamentally different things.** AUC measures how well the model separates lytic from
   non-lytic pairs across the entire probability range. Top-3 measures whether the correct phages end up in the top 3
   slots per strain. A feature that improves average-case AUC can hurt top-3 by pushing a marginally-higher-scoring
   wrong phage above a correct one. `legacy_label_breadth_count` is exactly this kind of feature — it improves the
   average but dilutes the per-strain signal.

3. **The pairwise block (Track E) was correctly excluded.** No subset containing pairwise cleared the AUC gate while
   improving top-3. The genus-level receptor lookup (TE01) was too coarse — 80% coverage but no within-genus
   specificity. The defense evasion proxy (TE02) showed up at SHAP rank #6 globally, suggesting the signal exists but
   the current encoding doesn't capture it well enough. Worth revisiting in v2 with finer-grained lookups (per-species
   or per-RBP-family), but not worth forcing into v1.

4. **The winner selection rule was too conservative.** Two 2-block arms (defense + phage-genomic, OMP + pairwise) both
   hit 89.2% top-3 but were disqualified because their AUC was 0.001–0.003 below the all-features reference. On 65
   holdout strains, AUC differences this small are noise. The rule prevented the sweep from selecting these
   higher-ranking arms. For future sweeps, the AUC gate should use a tolerance (e.g., AUC ≥ reference − 0.005) rather
   than a strict ≥.

5. **Two numbers for partners.** The locked v1 model gives CDMO partners two honest benchmarks:
   - Panel evaluation (87.7% top-3, AUC 0.911): what the model does on fully characterized strains
   - Novel-strain prediction (92.3% top-3, AUC 0.835): what to expect for a new clinical isolate with only a genome

   The 92.3% clears our 90% target. The 0.835 AUC means per-pair probability estimates are less reliable for novel
   strains, so recommendations should be presented as a ranked shortlist, not as calibrated P(lysis) values.

#### Track P acceptance criteria updated

All three Track P tasks (TP01, TP02, TP03) now require presenting both the panel model and the deployment-realistic
model side by side. Partners need to see both: the panel model validates the approach, the deployment-realistic model
sets honest expectations for their use case.

### 2026-03-23: Label leakage invalidates v1 model — plan restructured

#### Executive summary

Review of the TG04 SHAP results revealed that the v1 model's two strongest features (`legacy_label_breadth_count` with mean
|SHAP| 1.18 and `legacy_receptor_support_count` at 0.57) are derived from training labels, not independent
inputs. The v1 "panel-default" model is not predicting — it is memorizing. The plan has been restructured to delete these
features, retrain from scratch, and report whatever comes out as the honest v1 baseline. The prior dual-arm
panel/deployment framing is abandoned: there is only one model, the leakage-clean one.

#### What was decided

1. **Label-leaked features are a bug, not a variant.** `legacy_label_breadth_count` is literally "how many phages
   lyse this host" repackaged as a feature. `legacy_receptor_support_count` counts training-positive pairs per receptor
   cluster. Both encode the answer. Keeping them as an "optional panel-only arm" would be dishonest.

2. **Track P deleted.** All three presentation artifacts (digital phagogram, coverage heatmap, feature lift
   visualization) were designed around the dual-arm leaked model. The code, tests, and lab notebook have been removed.
   Visualizations will be rebuilt from scratch after the clean model is established.

3. **Track I made a dead end.** All 10 Track I tasks are done, but no Track I output feeds into any downstream track.
   TI09/TI10 count rows but never retrain a model with external data. Track I remains in the plan as completed work but
   nothing depends on it until external data is actually wired into model training.

4. **Track G extended with four new tasks:**
   - TG06: Delete leaked features from ST02, Track E, and Track G code
   - TG07: Retrain, recalibrate, re-run SHAP and ablation on the clean feature set
   - TG08: Re-run Track F evaluation and Track H recommendations, verify Track J end-to-end
   - TG09: Investigate whether non-leaky features can close the ~7.6pp AUC gap

5. **Track J depends only on G now** (removed F, H, I from depends_on to match what TJ01 actually runs).

6. **Track F and H descriptions updated** to note their current metrics are invalidated and will be re-run as part of
   TG08.

#### Why the prior dual-arm framing was wrong

The 2026-03-22 project entry framed two numbers for partners: panel evaluation (AUC 0.911) and novel-strain prediction
(AUC 0.835). The implicit message was "the model works great on known strains and somewhat worse on novel ones." The
actual message should have been: "the model's best feature is the training labels in disguise, and the 0.911 AUC is
inflated by that leakage." The deployment-realistic arm was the honest model all along — it should have been the only
model from the start.

#### What metrics to expect after cleanup

The deployment-realistic numbers from TG05 (top-3 92.3%, AUC 0.835, Brier 0.158) are the current best estimate for the
clean model, but the actual TG07 retrain may produce different numbers since the feature pipeline itself changes (not
just column exclusion at prediction time). TG09 will investigate whether the AUC gap can be partially closed with
non-leaky features.

### 2026-03-24: Clean v1 model locked — `defense + phage_genomic` is the honest baseline

#### Executive summary

Post-merge review of TG06-TG08 found two additional problems: LightGBM nondeterminism causing the sweep winner to flip
across runs, and 5 out of 13 pairwise features being derived from training labels (soft leakage). The v1 winner is now
locked to `defense + phage_genomic` — the 2-block arm that excludes all label-derived pairwise features. LightGBM
determinism will be fixed and the lock file will be treated as a human decision rather than a regenerated output.

#### Current honest v1 numbers

These are from the TG08 Track J end-to-end regeneration (the most realistic run):
- Winner: `defense + phage_genomic`
- Holdout ROC-AUC: ~0.837
- Holdout top-3 hit rate: 90.8%
- Holdout Brier: ~0.160

The 90.8% top-3 meets the 90%+ target. The AUC is lower than the old leaked model (0.911) but honest.

#### Remaining leakage in the codebase

The pairwise block (Track E) still contains label-derived features that are not yet deleted:
- TE02: all 4 `defense_evasion_*` features (collaborative filtering on training lysis rates)
- TE01: `receptor_variant_seen_in_training_positives` (binary flag from training positives)

These are excluded from the v1 model by locking to `defense + phage_genomic`, but the code still exists. A future task
should evaluate whether the clean pairwise features (TE03 distances, TE01 curated lookups) add value individually
without the label-derived ones.

#### Plan updates

- TG09: Fix LightGBM determinism, lock `defense + phage_genomic`, separate sweep from Track J regeneration
- TG10: Re-run downstream verification on the stable 2-block lock
- TG11 (was TG09): Investigate calibration gap — now aware of pairwise soft leakage, should evaluate clean pairwise
  features individually

#### Future: clean up upstream soft leakage in Track E

After TG09-TG11 are done and the v1 baseline is stable, consider a follow-up pass to delete or gate the training-label-
derived features in Track E itself: TE02's `defense_evasion_*` collaborative filtering features and TE01's
`receptor_variant_seen_in_training_positives`. These are currently excluded from v1 by the 2-block lock, but the code
still produces them. Deleting them would make the feature pipeline honest by construction rather than by configuration,
and would prevent future sweep runs from accidentally including them in a winning arm.

#### Future: workflow tooling (Snakemake / Nextflow) — not yet

Considered reimplementing the Track G pipeline in Snakemake for declarative dependencies, automatic skip of up-to-date
steps, and built-in parallelism. Decision: defer. The pipeline is still in flux (TG09-TG11 pending), the pain points
are scientific (leakage, nondeterminism) not operational, and the pipeline is small enough (~5 steps, ~30 min end-to-end,
single machine) that the linear `run_track_g.py` / `run_track_j.py` runners are sufficient. Revisit when: (a) external
data (Track I) is wired in and multiple data variants need training, (b) the pipeline grows beyond ~10 steps with
nontrivial branching, or (c) cluster execution or robust checkpointing becomes necessary.

### 2026-03-24: External-data decision locked for v1

#### Summary

Track K is now closed at the strategy level: the v1 model remains trained on internal data only. TK01 found no
joinable VHRdb rows in the available TI08 artifact, and the fixture-based TK02-TK05 follow-up arms were all neutral on
ROC-AUC, top-3 hit rate, and Brier score relative to the locked internal-only baseline.

#### Decision

- `lyzortx/pipeline/track_g/v1_feature_configuration.json` now records `external_data_lock_task_id: TK06`,
  `locked_training_data_arm: internal_only`, and `locked_external_source_systems: []`.
- No final-model retrain was performed for this task because the promotion condition was not met.
- Future production reruns can still revisit the decision through TK06 once real Track I / Track K manifests exist, but
  the current repo state does not justify changing the v1 release contract.

### 2026-03-24: Track I and Track K completed on empty data — all reopened

#### Executive summary

Post-merge review of TK01-TK05 revealed that all Track K tasks reported zero deltas because Track I never downloaded
any external data. The entire TI03-TI10 chain built plumbing that reads from nonexistent files. Track K inherited the
emptiness and reported "neutral" lift based on zero external rows. TK06 (PR #217) was rejected because it would have
locked an "internal-only" decision based on zero evidence. All 14 tasks (TI03-TI10, TK01-TK06) have been set back to
pending with acceptance criteria that require >0 real rows at every stage.

#### Root cause

1. **No Track I step downloads external data.** The code reads from local paths that were never populated. No HTTP
   requests, no API calls, no `urllib` — the download step was never implemented.
2. **Track K silently tolerated missing TI08 output.** `build_vhrdb_lift_report.py` line 258:
   `cohort_rows = read_csv_rows(...) if path.exists() else []` — silent fallback to empty list.
3. **Agents marked tasks done despite zero results.** Lab notebooks openly acknowledged "0 joinable rows" and
   "validated on a minimal fixture" but tasks were closed as completed anyway.
4. **CI starts with no generated outputs.** The agents ran in GitHub Actions where `lyzortx/generated_outputs/` does not
   exist. Without fail-fast on missing data, every step that depends on generated outputs silently produces nothing.

#### What changed

- New AGENTS.md rules: fail-fast on missing data, substance over plumbing, CI environment note
- TI03-TI06 now require downloading real data from source URLs and producing >0 rows
- TI07-TI10 now require >0 real external rows at each processing stage
- TK01-TK05 now require >0 external rows in the augmented training set
- TI03-TI07 upgraded to gpt-5.4 (external service integration needs research judgment)
- TK06 (PR #217) rejected — cannot synthesize results that don't exist

#### Lessons

- **Zero results is not a finding, it's a failure.** A task that runs to completion on empty inputs and reports zero
  deltas has failed its acceptance criteria, even if the code ran without error.
- **Silent fallback is not graceful degradation.** Code that returns empty results when inputs are missing is hiding a
  bug, not handling an edge case. Raise on missing data.
- **Acceptance criteria must include data volume assertions.** ">0 rows" is the minimum bar. Without it, an agent can
  build correct plumbing that processes nothing and call it done.

### 2026-03-24: V2 plan restructure — external data dead end, mechanistic features forward

#### Executive summary

Deep analysis of Track I and Track K revealed that external data integration is a dead end with the current pipeline:
only VHRdb overlaps with the internal 404×96 panel (23,885 pairs), and that's the paper's own data uploaded to VHRdb by
the original authors (datasource 257). BASEL, KlebPhaCol, GPB, Virus-Host DB, and NCBI all have zero strain overlap.
Track K was deleted. Track I was trimmed to download-only (TI01-TI06). A new Track L (Mechanistic Phage Features)
replaces the deleted label-derived pairwise features with annotation-based features from Pharokka. A label policy fix
(TA11) captures the +3.1pp top-3 improvement from downweighting borderline `matrix_score=0` noise positives.

#### Key findings from VHRdb analysis

1. **VHRdb overlap is circular.** The 23,885 overlapping pairs are the paper's own experiment (datasource 257: "Host
   range of the 96 coliphages from the Antonina Guelin collection on the 403 natural isolates of Escherichia from the
   Bertrand"). VHRdb compressed the paper's 0-4 score to 0-2 for upload.
2. **Perfect agreement on matrix scores 1-4.** Zero off-diagonal entries in the cross-tabulation. VHRdb's 3-level
   scale is a lossless compression of the 0-4 scale for non-zero scores.
3. **1,737 disagreements are all `matrix_score=0, label=1`.** These are pairs where 1-3 replicates out of 8-9 showed
   lysis (noise), the matrix aggregated to 0, but the "any lysis" label policy called them positive. VHRdb correctly
   reports "No infection" for all of them.
4. **No other VHRdb datasource shares our phages.** Only datasource 153 (3 LF82 phages on ECOR) has any panel phage
   overlap — 3 out of 96 phages, across 73 ECOR strains with no panel overlap.

#### What changed in the plan

- **Deleted Track K** — all code, tests, plan entries. Archived notebook to `archive_v1/track_K.md`.
- **Trimmed Track I** to TI01-TI06 (download infrastructure). Deleted TI07-TI10 code (confidence tiers, training
  cohorts, ablations, lift analysis) — these depended on external data being trainable.
- **Added Track L** (Mechanistic Phage Features) — TL01-TL06: set up bioinformatics env, annotate phages with
  Pharokka, build mechanistic RBP-receptor and defense-evasion features, retrain, validate on external strains.
- **Added TA11** (label policy fix) — downweight `matrix_score=0` noise positives based on the VHRdb finding.
- **Track E kept as-is** — Track G depends on it; Track L will eventually supersede it.

#### Pharokka POC results (local, 2026-03-24)

- Pharokka 1.9.1 on LF82_P8: 276 CDS annotated in 2m 43s
- 29 tail genes (including long tail fiber proximal/distal subunits, baseplate components)
- 7 lysis genes (holin, spanins, lysis inhibitors)
- 1 anti-restriction nuclease
- 140 hypothetical proteins (51% unannotated — typical for phage genomes)
- Extrapolation: ~1 hour for all 97 phages at 4 threads

#### Future: External validation dataset for generalized inference pipeline

**Trigger condition:** Revisit when a full interaction matrix (with both positives and negatives) for novel E. coli
strains becomes available. TL09 covers positive-only validation via VHdb; this note is about the stronger full-matrix
validation that would allow AUC and top-3 computation.

**Context:** The v2 plan adds a generalized inference pipeline (TL06-TL08) that accepts arbitrary E. coli genomes and
arbitrary phage FNAs, computes features from sequence, and predicts lysis. The locked v1 model uses only defense subtypes
(79 features from Defense Finder) and phage k-mer SVD (26 features from tetranucleotide frequencies) — both
sequence-derivable, so the pipeline is architecturally sound. TL09 validates on VHdb positive-only pairs, but
full-matrix out-of-distribution validation (with negatives) requires a dataset that doesn't yet exist in public.

**Why existing sources don't provide a full matrix for novel strains:**

- **ECOR strains** are already in the 404-strain training panel (71 of 404). Not novel.
- **VHRdb pair-table overlap** (26,029 rows) is the paper's own data (datasource 257). Novel VHRdb pairs (30,643) failed
  entity resolution. VHRdb strain-level pairs (~500 positives across ~70 hosts) are usable for positive-only validation
  (TL09) but contain no negatives.
- **BASEL** has labeled E. coli interactions but only 4 E. coli host strains. Phages are novel (78 genomes on NCBI,
  MZ501046-MZ501113) but the host count is too small for meaningful metrics.
- **KlebPhaCol** is Klebsiella, not E. coli. Wrong species.
- **GPB** (Gut Phage Biobank) has only 1 E. coli strain (RTGS0219) out of 40 hosts. Phage genomes on CNGB, not NCBI.
- **SNIPR001** (Nature Biotech 2023) has 429 E. coli strains x 162 phages with interaction matrix on GitHub, but strain
  genomes are commercially restricted.

**What a full-matrix validation dataset requires:**

1. **E. coli strains not in our 404-strain panel** — at least 20-30 strains for meaningful ranking metrics.
2. **Full interaction matrix** (lysis AND no-lysis from spot assay or equivalent) — not just positive pairs.
3. **Genome assemblies for the host strains** — so TL07 can run Defense Finder.
4. **Phage genomes (FNA)** — so TL06 can project k-mer features.
5. **Sufficient phage panel size** — at least 10-20 phages per host for top-3 ranking to be meaningful.

**Partial solution found (2026-03-28): Virus-Host DB positive-only validation.**

Virus-Host DB contains ~70 E. coli strains at strain-level resolution (excluding lab strains) with ~900 unique phage
genome accessions and ~500 positive pairs (phage confirmed to lyse host, from literature). Most hosts have NCBI genome
assemblies. This gives a positive-only validation path: download host assemblies + phage FNAs, run generalized inference,
check that known-positive pairs get high predicted P(lysis). Limitations: no negatives, so AUC and top-3 hit rate cannot
be computed. Metrics are limited to recall-on-positives, calibration-on-positives, and rank-of-positives-vs-random. This
is implemented as TL09 in the plan.

**Full interaction matrix validation still requires:**

- Published phage therapy studies with E. coli host-range matrices and deposited genomes (both host and phage).
- Phage biobanks that publish interaction data alongside NCBI accessions (e.g., future GPB or DSMZ releases).
- Direct collaboration with labs running E. coli phage screening panels.
- The SNIPR001 dataset (Nature Biotech 2023) has 429 E. coli strains x 162 phages with interaction matrix on GitHub,
  but strain genomes are commercially restricted.

**What to do when a full-matrix candidate is found:**

- Verify E. coli strain count, phage count, and label quality before committing to ingestion.
- Download host genome assemblies from NCBI, run TL07/TL08, compare predictions against labels.
- Report AUC, top-3 hit rate, and calibration metrics as the honest out-of-distribution benchmark.
- Compare against in-panel holdout metrics to quantify the generalization gap.

#### Future: Structural RBP-receptor prediction via protein folding

**Trigger condition:** Revisit if TL02 (RBP-receptor compatibility from Pharokka annotations) hits the escape hatch —
i.e., PHROG functional annotations are too coarse to map RBPs to specific host receptor targets (FhuA, BtuB, OmpC,
LPS core, etc.).

**The problem:** Pharokka labels phage genes by PHROG family ("tail fiber protein," "tail spike protein") but doesn't
tell you which receptor the RBP binds. That mapping is determined by the 3D structure of the RBP tip domain and its
binding interface, not by sequence homology alone.

**Potential approach — structure-based RBP clustering:**

1. Predict 3D structures for all phage RBPs using AlphaFold2 or ESMFold.
2. Extract receptor-binding tip domains (C-terminal regions of tail fibers / tail spikes).
3. Cluster RBPs by structural similarity of the tip domain. RBPs with similar tip folds likely target the same receptor
   class.
4. Map each structural cluster to a receptor class using known reference structures from the literature (e.g., T4 long
   tail fiber tip → OmpC, T5 pb5 → FhuA).
5. For each phage-host pair, compute a compatibility feature: does this phage's RBP cluster match a receptor present on
   the host (from Track C OMP data)?

**Precedent:** The BASEL phage collection papers (Dunne et al. 2021, Maffei et al. 2025) used structural analysis to
validate receptor assignments for their phage panels. Structure-based receptor grouping is scientifically established.

**Why not now:**

- Heavy computational dependency (AlphaFold or ESMFold, ~200-300 structure predictions).
- Requires curated reference structures to anchor the cluster→receptor mapping.
- Would be a separate track-level effort, not a subtask.
- TL02's escape hatch already handles the annotation-only failure case gracefully.

**What to do if triggered:** Scope a new track (e.g., Track M: Structural RBP-Receptor Prediction) with explicit
acceptance criteria for structure prediction, tip domain extraction, clustering, and receptor assignment. Evaluate
whether the expected lift justifies the computational cost before committing.

#### Future: Depolymerase domain annotation for capsule-specificity matching

**Trigger:** Enrichment analysis (TL02) shows that host capsule/LPS features carry signal for lysis prediction, but
pharokka's generic annotations ("polysaccharide chain length determinant protein", "tail spike protein") are too coarse
to distinguish which capsule types a phage can degrade.

**Context (2026-03-29):** TL01 found 18 "polysaccharide chain length determinant" genes across the panel, but pharokka
does not annotate capsule-type specificity. Dedicated depolymerase tools (DepoScope, DePP, PDP-Miner) only classify
binary depolymerase yes/no — none predict capsule-type targets. The most promising approach is running Pfam/InterPro on
tail spike protein sequences to identify glycosyl hydrolase families, which can be mapped to polysaccharide substrate
classes via literature. Track C has LPS core type (~5 types, well-covered) and Klebsiella capsule type (94% missing), so
the host-side signal may be sparse.

**What to do if triggered:** Run `hmmscan` against Pfam-A on the pharokka `.faa` tail spike sequences (already available
in the per-phage output). Extract glycosyl hydrolase family annotations and map to substrate specificity using CAZy
database cross-references. If enough phages carry classifiable depolymerase domains and enough hosts have capsule type
data, add depolymerase-capsule compatibility as a pairwise feature alongside RBP-receptor features.

#### Future: External interaction matrices as additional training data

**Trigger:** The model architecture is validated on the internal 97-phage panel and generalizes to external phage/host
genomes via the annotation-based feature pipeline (TL05-TL08 complete).

**Context (2026-03-29):** Because Track L features (RBP PHROGs, anti-defense genes, host receptors, defense systems) are
defined by universal sequence databases (PHROGs, Pfam, OMP BLAST clusters, DefenseFinder), any external phage-host
interaction dataset can be encoded into the same feature space. External phage genomes get pharokka-annotated and their
RBPs assigned to the same PHROG families; external host genomes get receptor-typed and defense-system-annotated with the
same pipelines. The external interaction matrix then provides additional (PHROG, receptor) → lysis observations that
strengthen the learned associations.

**What to do if triggered:** Identify external interaction datasets (e.g., from NCBI, PhageScope, published
supplementary tables). Run the Track L annotation + Track C host typing pipelines on the external genomes. Pool the
encoded interactions with internal training data and retrain. Measure lift from the expanded training set.

### 2026-03-30: Track L replan — enrichment holdout leak and external validation failure

#### Executive summary

Post-completion review of Track L identified two independent problems: (1) TL02's enrichment analysis includes ST03
holdout strains in the interaction matrix, leaking test data into TL03/TL04 feature weights; (2) the TL08 genome-only
inference bundle fails external validation on Virus-Host DB (known positives score below random pairs). These are
separate issues — the leak affects in-panel evaluation trustworthiness, the external failure reflects an architectural
gap in the deployable feature set.

#### What changed in the plan

- **TL10 added**: Fix the enrichment holdout leak in `run_enrichment_analysis.py` by filtering out ST03 holdout bacteria
  before calling `compute_enrichment()`. Mechanical fix — the permutation test is sound, only the input selection is
  wrong.
- **TL03/TL04/TL05**: Remain marked done for now. After TL10 lands, they will need re-evaluation (separate tickets) to
  determine whether enrichment features provide any honest lift with holdout-clean weights.

#### Diagnosis summary

Three compounding issues in Track L:

1. **TL02 enrichment circularity (confirmed, 98% confidence)**: `run_enrichment_analysis.py` loads
   `label_set_v1_pairs.csv` (all 369 bacteria) with no holdout filtering. `compute_enrichment()` has no holdout
   parameter. The TL02 acceptance criteria explicitly said "uses the full interaction matrix." Fixed by TL10.

2. **TL08 genome-only bundle is feature-impoverished (confirmed, 80% confidence)**: TL08 trains on only defense subtypes
   (79 features) + phage k-mer SVD (26 features). It drops the entire v0 metadata block (serotype, phylogroup,
   morphotype, ~20 categorical features), OMP receptor variants, UMAP embeddings, and all phage taxonomy. Many of these
   are genome-derivable but were not wired into TL08. The EDL933 round-trip failure (median P(lysis) delta 0.14, 9/96
   ranks matching) is primarily explained by this feature-set mismatch, not by defense annotation lossiness.

3. **No pairwise compatibility signal in deployable model (moderate confidence, 70%)**: Neither defense subtypes nor
   k-mer SVD encode how a specific phage interacts with a specific host surface. The model memorizes panel-specific
   defense-profile × k-mer-profile combinations that don't transfer to novel genomes. Adding genome-derivable
   compatibility features (OMP receptors, enrichment-weighted PHROG × receptor pairs) to the inference bundle is the
   most promising direction, but requires TL10 first to make the enrichment weights honest.

#### What is still sound

- Pharokka annotation pipeline (TL01): correct and biologically reasonable.
- Enrichment module statistics (TL02 implementation): permutation test, BH correction, phage conditioning all sound.
- Generalized inference architecture (TL06–TL08): plumbing is correct, transform persistence works.
- Novel organism projection helpers (TL06, TL07): tested and correct.
- Leakage diagnosis from Track G (TG04–TG12): rigorous and honest.

### 2026-03-30: Track L replan follow-up — tighten downstream acceptance criteria

#### Executive summary

Reviewed the TL03-TL09 PRs, the failed TL03 Codex implement run, and the Track L notebook entries to identify which
mistakes were implementation accidents versus plan-specification failures. The main issue was under-specified acceptance
criteria: several tasks allowed a technically plausible but strategically wrong completion to count as done. The plan
now adds TL11-TL14 to make the next pass fail fast on the specific mistakes we can already anticipate.

#### What changed in the plan

- **TL11 added**: rebuild TL03/TL04 from TL10's holdout-clean enrichment outputs and emit manifests proving which split
  and which excluded bacteria IDs were used.
- **TL12 added**: rerun mechanistic lift with bootstrap confidence intervals and a predeclared lock rule, so tiny noisy
  deltas can no longer justify a new v1 configuration.
- **TL13 added**: rebuild the deployable bundle under an explicit feature-parity audit and self-contained-artifact
  contract; silently dropping training-time feature blocks no longer counts as success.
- **TL14 added**: rerun external validation under a strict cohort contract with a required multi-host round-trip check,
  so "validation succeeded as a script run" is separated from "bundle actually generalized."

#### Mistakes these new tasks are meant to prevent

1. **Leaked-but-plausible mechanistic rebuilds**: TL03/TL04 were scientifically framed correctly, but nothing in their
   original follow-on criteria would have forced the re-evaluation to prove it was using TL10-clean inputs rather than
   stale leaked enrichment CSVs.

2. **Locking on noise**: TL05 proposed a new mechanistic lock from holdout deltas that were already within noise on a
   65-strain holdout. The next task must report uncertainty and apply a stated decision rule before it is allowed to
   promote any arm.

3. **Plumbing without feature-parity honesty**: TL08 proved that the bundle machinery works, but it was allowed to ship
   a genome-only model that omitted many training-time signals without first surfacing that gap as a formal feature
   parity audit.

4. **Validation with too-weak round-trip guarantees**: TL09 named several panel-host examples, but only `EDL933` was
   actually comparable through the saved reference artifact. Future validation tasks must pre-materialize the cohort and
   treat round-trip host count as a gate, not a hope.

5. **Review-thread fixes discovered too late**: several issues were caught only in PR review rather than by the task
   contract itself, including weak negative fixtures (TL04), late cache short-circuiting (TL07), hardcoded panel paths
   in the bundle (TL08), and under-specified parsing/selection semantics (TL09). TL11-TL14 now encode those lessons
   directly as acceptance checks.

#### Additional note from Codex workflow logs

The first TL03 Codex implement run failed before any repo code ran because `conda env create -f environment.yml` was
unsatisfiable on CI (`openjdk` / `fontconfig` solver conflict). That is not a TL03 logic bug, but it is a reminder that
bioinformatics-heavy tasks should explicitly require CI-compatible environment resolution rather than assuming local and
CI solvability stay aligned.

### 2026-03-30: TL12 mechanistic lift rerun came back as no honest lift

TL12 re-ran the mechanistic lift evaluation on the TL11 holdout-clean feature rebuild and added a stricter lock rule:
only promote a mechanistic arm if the paired bootstrap 95% CI for ROC-AUC delta vs the locked baseline is entirely
above zero, with no material degradation in top-3 hit rate or Brier score. The live rerun used the same label set and
code path as TL05, zero-filled missing TL11 pair rows on the holdout side, and produced `no honest lift`.

Key holdout results:

- Baseline `defense + phage_genomic`: ROC-AUC `0.835466`, top-3 `0.892308`, Brier `0.146153`.
- `+TL03`: ROC-AUC `0.820052`, top-3 `0.846154`, Brier `0.148296`.
- `+TL04`: ROC-AUC `0.838029`, top-3 `0.892308`, Brier `0.144594`.
- `+TL03+TL04`: ROC-AUC `0.822875`, top-3 `0.861538`, Brier `0.147016`.

The ROC-AUC delta CI for TL04 was still `[-0.002322, 0.007935]`, so it never cleared the lock threshold even though
its point estimate was the best of the mechanistic arms. TL03 and the combined arm were clearly worse. The correct
interpretation is that the pairwise enrichment path remains exploratory only; it is not a v1 lock candidate.

### 2026-03-31: TL12 follow-up patch confirmed the same outcome on a hardened rerun path

The follow-up patch for issue `#280` did not change the scientific conclusion, but it did materially harden the rerun
path that produced it. The fixes were:

- restrict TL11 zero-fill semantics to holdout evaluation rows only, while keeping non-holdout pair joins fail-fast;
- validate TL11 provenance from the actual CLI-supplied manifest paths rather than assuming default sibling manifests;
- rebuild stale default TL11/TL02 artifacts automatically before the rerun instead of trusting pre-replan generated
  outputs; and
- remove duplicate CV fold training plus add phase/progress logging so the rerun is observable and faster.

After those fixes, TL05 was rerun end-to-end again and still landed on `no honest lift`.

Key holdout results from the hardened rerun:

- Baseline `defense + phage_genomic`: ROC-AUC `0.837060`, top-3 `0.907692`, Brier `0.159486`.
- `+TL03`: ROC-AUC `0.822245`, top-3 `0.907692`, Brier `0.156530`.
- `+TL04`: ROC-AUC `0.839504`, top-3 `0.892308`, Brier `0.156965`.
- `+TL03+TL04`: ROC-AUC `0.823165`, top-3 `0.892308`, Brier `0.155707`.

The best mechanistic arm was still TL04, but the stricter lock gate still rejected it: ROC-AUC delta vs baseline was
`+0.002444` with paired bootstrap CI `[-0.002404, 0.007416]`, and top-3 delta was `-0.015384` with CI
`[-0.047648, 0.025000]`. So the implementation fixes improved correctness and reproducibility, not the model story.

### 2026-03-31: Plan update after TL12 — stop treating mechanistic pairwise rebuilds as a pending v1 lock path

Track L was replanned again after the hardened TL12 rerun. The key change is strategic, not technical: TL03/TL04 are
no longer treated as a promising near-term route to a better locked panel model. TL12 answered that question twice
already, once on the first clean rerun and once after the hotfixes to the rerun path, and both answers were the same:
`no honest lift`.

So the plan now separates two stories that had been mixed together:

- **Dead-ended for the current v1 lock**: annotation-derived mechanistic pairwise features as a replacement for the
  locked panel model.
- **Still alive**: generalized inference, but only if a richer deployable bundle can add real genome-derivable
  compatibility signal and improve round-trip behavior first.

Concretely, TL13 is now framed as a feature-parity audit plus deployable compatibility experiment, and TL14 is now
gated on TL13 proving that the richer bundle actually improves round-trip behavior on panel hosts. If TL13 cannot
clear that gate, the right next step is another replan, not external validation by inertia.

### 2026-03-31: Track L replan follow-up — new preprocessing tasks before any further deployable-bundle claims

The next Track L change is a planning correction, not a new scientific result. Review of TL13/TL14 concluded that the
feature-parity audit was too absolute about deployability: several blocks were labeled "not deployable" when the more
honest statement was "not currently implemented as a preprocessing step."

The plan now adds four new tasks:

- `TL15`: raw-host surface projector for OMP/LPS-style compatibility features;
- `TL16`: genome-derived host typing projector for the reproducible subset of the old metadata block;
- `TL17`: phage-side deployable compatibility preprocessor beyond k-mer SVD; and
- `TL18`: rebuild the deployable generalized-inference bundle only after `TL15`-`TL17` are all available.

Two strategic decisions matter here:

1. `TL15`-`TL17` are parallel preprocessing tasks, not another serial chain hidden inside one track.
2. The plan deliberately does **not** add a new external-validation ticket yet. Another validation pass is only worth
   discussing after `TL18` proves that the richer deployable bundle improves round-trip behavior on panel hosts.

The updated plan also records one negative lesson explicitly: do not treat fitted UMAP host coordinates as the next
deployable shortcut. If continuous host-similarity signal is still needed, it should come from a stable runtime
distance or projection contract, not from reusing a fragile embedding fit.

### 2026-04-02: Deployment-Paired Feature Pipeline — new track after TL18 audit

#### Executive summary

A step-by-step inference audit of the TL18 richer bundle found three deployment-path bugs and a systematic feature
redundancy problem that together motivate a new track: the Deployment-Paired Feature Pipeline (DEPLOY01-07). The
primary goal is deployment integrity — the model should be trained on exactly the features it will see at inference
time. The secondary goal is richer features that give the model more information to work with.

#### What the audit found

The TL18 audit walked through every inference step on the 3 committed validation hosts (55989, EDL933, LF82) and
the 65-strain holdout set. It confirmed that the model quality improvement from the richer bundle is real (ROC-AUC
+0.036, 95% CI [+0.002, +0.079], 98.5% probability of improvement; top-3 +3.2pp; Brier -0.011). No label leakage was
found. But three deployment-path bugs were identified:

1. **DefenseFinder version drift**: fresh runs on raw FASTAs disagree with panel annotations (EDL933: 8 vs 9 systems).
   Defense subtypes carry 17.3% of model importance.
2. **Capsule train/inference mismatch**: raw HMM scan is more sensitive than the picard `Capsule_ABC` metadata field.
   EDL933 gets K4 from raw HMM but `Capsule_ABC=0` in panel metadata. Capsule features carry 3-5% of importance.
3. **Extra phage 411_P3**: 97 FNA files vs 96 panel phages. Minor ranking impact.

A redundancy analysis found 91 wasted one-hot features from exact duplicates (`host_o_type` duplicates
`host_o_antigen_type`: 84 one-hot columns; `host_surface_lps_core_type` duplicates `host_lps_core_type`: 6 columns;
`host_capsule_abc_present` duplicates `host_capsule_abc_proxy_present`: 1 column) plus 4 derived summary features
computable from their constituents (`defense_diversity`, `has_crispr`, `abi_burden`, `rbp_family_count`).

A feature encoding audit found that 141 of 190 numeric features (74%) are binary thresholds that discard continuous
information. Receptor phmmer bit scores, RBP family mmseqs percent identity, capsule HMM profile scores, and defense
system gene counts all carry biological gradients that the current binary encoding throws away.

#### Strategic decision

The new DEPLOY track addresses all of these by:

1. Downloading all 403 Picard collection assemblies from figshare (doi:10.6084/m9.figshare.25941691.v1)
2. Re-deriving all host features from raw FASTAs using the same pipeline that runs at inference time
3. Switching to continuous scores for receptors (phmmer bit scores), RBP families (mmseqs % identity), and capsule
   (per-profile HMM scores). Defense subtypes switch from binary to integer gene counts.
4. Dropping 91+ redundant features
5. Retraining and evaluating with a 3-way comparison: TL18 panel baseline vs parity-only vs parity+gradients
6. Wiring the re-derived pipeline into the inference runtime with zero-delta parity validation

The DEPLOY track replaces Track L as the active development front. Track L is complete.

#### Why not just retrain on the current features

The parity fix is necessary for deployment integrity regardless of whether it improves metrics. A model trained on
curated panel metadata but scored on fresh bioinformatics tool outputs is predicting on a different feature distribution
than it was trained on. Bugs #1 and #2 prove this causes real divergence. Even if the holdout metrics look the same
after the parity fix (because holdout already used consistent panel features), the deployed predictions for novel hosts
will be more trustworthy.

The gradient features are a bonus — they may or may not improve metrics. DEPLOY06 explicitly separates the two effects
with an ablation: train one model with parity-only (raw-derived binary) and one with parity+gradients (raw-derived
continuous) to measure whether gradients help.

#### Infrastructure findings during DEPLOY planning

The figshare assemblies (doi:10.6084/m9.figshare.25941691.v1, CC BY 4.0) are 403 FASTA files totalling 1.9GB. The
"Download all" zip endpoint downloads in ~7 minutes and unzips in ~7 seconds. The `full-bio` CI image is 9.5GB
compressed, which already strains the free GitHub Actions runner (14GB disk, ~6GB usable). Baking 1.9GB of assemblies
into the image is not viable. Instead, DEPLOY tasks download assemblies on demand at runtime, with a skip-if-present
check so local dev only pays the 7-minute cost once. If the free runner's disk proves insufficient during DEPLOY02
(DefenseFinder on 403 genomes), the fallback is the paid `ubuntu-24.04-4core` runner at ~$0.24-0.48 per task run.

Pharokka meta mode (`--meta --split`) annotates all 97 phage genomes in 3 minutes 13 seconds, compared to 1-2 hours
with the current per-phage `ProcessPoolExecutor` approach. The 40x speedup comes from running mmseqs2 database indexing
and profile search once instead of 97 times. The existing 194 committed pharokka annotation files (4.5MB in
`data/annotations/pharokka/`) are small enough to keep as-is; the meta-mode optimization is deferred until the next
time `run_pharokka.py` is refactored.

#### Fork code placement policy enforcement

During DEPLOY planning we identified that ~1,024 files had been added outside `lyzortx/` by previous work (194
pharokka annotations, 223 capsule/phylogroup reference databases, 98 phage FNAs, 3 validation FASTAs, 6 defense
finder files, plus upstream dev/ artifacts). Root-level files (`AGENTS.md`, `CLAUDE.md`, env manifests) are justified
exceptions. The rest are a mix of legitimate reference data and generated outputs that should ideally live under
`lyzortx/`. Retroactive migration would touch hundreds of `DEFAULT_*_PATH` constants across dozens of files — not
worth the churn now. Going forward, all new data and code goes under `lyzortx/` (the DEPLOY assembly download path is
`lyzortx/data/assemblies/picard/`, gitignored).

### 2026-04-01: Project-level implication of the Track L raw-validation subset

The Track L replan is now stricter in one important way: "deployable" work can no longer hide behind the excuse that
the repo lacks raw-input fixtures or declared tool environments. We now have both:

- a committed host FASTA validation subset with provenance and checksums
- checked-in split env manifests for the host-typing and heavier bioinformatics toolchains

That does not finish generalized inference by itself, but it does change what counts as an honest intermediate result.
From this point on, `TL15`-`TL18` should be expected to validate at least part of their contract from raw FASTA inputs
available in a clean checkout, not only from precomputed intermediate tables. The project-level benefit is clarity:
future "deployable bundle" claims are now easier to falsify early if they still depend on hidden local state.

### 2026-03-31: TL14 implementation hardened the external-validation stop conditions

TL14 is now implemented as a strict gatekeeper instead of a permissive rerun of TL09. The code now reads the saved TL13
bundle contract, writes the exact validation cohort before any scoring, and records one of three explicit outcomes:
`deployable bundle validated`, `deployable bundle failed`, or `validation inconclusive because the cohort contract could
not be satisfied`.

This matters immediately because the saved TL13 round-trip artifacts recorded on `2026-03-30` still cover only
`EDL933`. That means the current repo state does **not** support a broad external-validation claim for the richer
deployable bundle. Until TL13's saved reference-backed panel cohort is rebuilt to at least 3 hosts, the honest project-
level reading is: **validation inconclusive because the cohort contract could not be satisfied**.

#### Future: richer defense features from DefenseFinder raw outputs

**Trigger**: revisit when DEPLOY07 evaluation is complete and the baseline defense-count features are locked.

DEPLOY06 checked in integer gene counts per defense subtype (79 columns), but DefenseFinder produces three additional
output files per host that are preserved locally in
`lyzortx/generated_outputs/deployment_paired_features/host_defense/`:

1. **Gene-level HMM scores** (`*_defense_finder_genes.tsv`): `hit_score`, `i_eval`, `hit_profile_cov`, `hit_seq_cov`,
   `sys_wholeness`, `sys_score`. These could yield per-subtype continuous detection-confidence features (analogous to
   the receptor phmmer scores in DEPLOY03) and system-completeness features.
2. **Sub-threshold HMM hits** (`*_defense_finder_hmmer.tsv`): all hits including those below DefenseFinder's
   co-localization threshold. Degenerate or partial defense systems that the current pipeline ignores entirely.
3. **Total defense gene count**: simple sum of defense-associated genes per host, a proxy for overall defense
   investment that the model currently cannot see.

The current integer-count encoding was chosen because the plan's design principle says "count is real biology, HMM
score is a tool artifact." This is approximately true — but system completeness (`sys_wholeness`) is not a detection
artifact; it reflects whether the full defense operon is intact. And sub-threshold hits may identify hosts with
degraded defense arsenals that look defenseless in the current encoding but still carry partial immunity.

**What to do**: after DEPLOY08 locks the baseline model, run a feature-importance comparison with and without the
richer defense features. If the additional features improve holdout metrics or shift SHAP rankings meaningfully, add
them as a DEPLOY track extension. If not, the current counts are sufficient and the extra complexity is not justified.
The raw outputs are already computed and preserved locally — no re-running of DefenseFinder is needed.

#### 2026-04-04: DEPLOY track renumbered — surface pre-compute inserted as DEPLOY07

The first Codex CI attempt at the full 403-host evaluation (old DEPLOY07, now DEPLOY08) failed because the DEPLOY03
surface derivation runs nhmmer at ~72s/host — infeasible on a 4-core CI runner. This mirrors the DEPLOY06 situation
where DefenseFinder was too slow for CI.

**Decision**: insert a new DEPLOY07 (surface pre-compute, `executor: human`) that runs locally and checks in the
aggregated surface CSV, following the same pattern as DEPLOY06 for defense. Old DEPLOY07 becomes DEPLOY08, old DEPLOY08
becomes DEPLOY09.

The runner (`run_all_host_surface.py`) uses pyhmmer for in-process HMMER and replaces the nhmmer DNA search with protein
phmmer (translating O-antigen alleles), cutting per-host scan time from 72s to 6.2s. Total wall time for 403 hosts on a
10-core Mac: ~10 min.

**Impact on DEPLOY08 (retrain/evaluate)**: DEPLOY08 now loads both pre-computed CSVs (defense + surface) and only needs
to run DEPLOY04 host-typing and DEPLOY05 phage-RBP derivation in CI. Those are fast enough for the Codex runner.

**Note on done-task references**: done tasks DEPLOY02-06 still reference "DEPLOY07" as the full evaluation task. Per
done-task immutability policy, those references are historical. The renumbering is tracked in git history and in the
track_DEPLOY.md notebook entry.

### 2026-04-04 20:45 UTC: Replanned AUTORESEARCH as a raw-input track with frozen featurizers

#### Executive summary

We replanned AUTORESEARCH around the raw interaction table plus host and phage FASTAs. Track A still supplies the
label policy, but DEPLOY outputs no longer define the model-search substrate; preprocessing is frozen in `prepare.py`
and the search loop is restricted to `train.py`. The benchmark contract is now explicit as well: the holdout is
bacteria-disjoint, the comparator artifact is named up front, and CI correctness checks are separated from full-scale
runtime measurement.

#### What changed in the plan

- Renamed the track intent to **Track AUTORESEARCH: Raw-FASTA Autoresearch**.
- Removed the track-level dependency on DEPLOY artifacts and replaced it with a dependency on Track A's label policy.
- `AR01` now freezes the raw corpus, label policy, FASTA inventory, comparator benchmark reference, and a bacteria-disjoint
  sealed split contract.
- `AR02` now freezes the sandbox and cache contract before any feature-family implementation starts.
- `AR03`-`AR06` split the cache build by runtime-risk boundary: host defense, host surface, host typing/stats, and
  phage projection/stats.
- `AR07` keeps the strict one-file search contract, but the model now searches over a raw-input cache instead of a
  DEPLOY-era artifact export.
- `AR08` still owns the dedicated RunPod workflow and environment boundary.
- `AR09` imports winners back through sealed-holdout replication before any promotion decision.

#### Why this is the right cut

The earlier AUTORESEARCH draft still treated DEPLOY-era outputs as the natural substrate for search. That would have
preserved too much of the old abstraction. If the goal is to find a better deployable learner, the substrate has to be
the raw corpus plus feature builders that can be rerun unchanged at inference time.

This replan keeps the useful engineering from DEPLOY without keeping DEPLOY as the scientific starting point:

- Track A supplies a fixed label policy instead of reopening the labeling debate.
- Raw-FASTA helper code survives when it is inference-safe.
- Panel-only metadata fields and metadata-derived proxies are cut, even if they once looked predictive.
- The `autoresearch` loop remains small and comparable run-to-run because expensive preprocessing stays frozen.

We also tightened the AUTORESEARCH acceptance criteria to reflect the runtime lessons already written down in
`track_DEPLOY.md`: DefenseFinder-scale preprocessing does not belong inside the search loop, repeated environment setup
is real overhead, and the wrong algorithmic shape can dominate wall time before the model trains. Those are now
explicit plan constraints instead of implicit tribal knowledge.

We then split the middle of the track further because "implement prepare.py" was still too coarse. Host defense,
surface, typing, and phage projection are not one risk surface; they fail with different dependencies, runtimes, and
optimization levers. The finer split should make orchestrator issues more honest and reviews more surgical.

We then tightened the dispatch contract after review: AR02 now fixes the cache schema composition boundary up front,
AR01 records the exact locked comparator benchmark for AR09 and requires bacteria-disjoint splits, AR07 names ROC-AUC
as the primary search metric, and AR03-AR06 explicitly separate CI-scale correctness checks from full-scale runtime
measurement. That should reduce both timeout risk and benchmark ambiguity before the first AUTORESEARCH issue is even
opened.

### 2026-04-05 09:30 UTC: Antiphage-landscape preprint supports adsorption-first modeling and weaker defense negatives

#### Executive summary

We reviewed the January 8, 2025 preprint _Protein and genomic language models chart a vast landscape of antiphage
defenses_ and its public repo. The paper is important for this project, but not as a new strain-level host-range
matrix. Its main project-level implication is that our adsorption-first strategy still looks right for _E. coli_ lysis
prediction, while current defense annotations should be treated as incomplete positive evidence rather than as
near-complete antiviral catalogs.

#### What changed in our interpretation

The 2024 Nature _Escherichia_ paper and the 2025 antiphage-landscape preprint answer different questions:

- the 2024 paper is about pairwise lysis prediction in _Escherichia_;
- the 2025 preprint is about discovering unknown antiphage systems in bacterial genomes.

These are not contradictory. The 2024 paper can still be right that adsorption factors dominate current supervised
prediction in _E. coli_, while the 2025 preprint is also right that the defense universe is much larger than today's
annotation tools capture.

#### Project-level implications

1. Keep adsorption and RBP features as the main modeling axis.
2. Keep defense features, but stop reading defense-feature absence as strong biological evidence.
3. Treat LM-based defense discovery as a candidate-generation resource, not as a reason to center the main classifier
   on defense burden.
4. The next likely high-value feature increment is finer phage-side RBP resolution from FASTAs, not a broad defense
   feature expansion.

#### Data implications

The preprint does not add a new strain-by-strain host-range matrix that can drop directly into training. Its wet-lab
validation is on _Streptomyces albus_ expressing candidate systems against a small phage panel. For supervised data
priority, the project should still rank:

1. highest-fidelity ingestion of the 2024 _Escherichia_ source package;
2. BASEL plus BASEL completion as same-host-genus external supervision;
3. PhageHostLearn as transfer/robustness data;
4. GPB and broader collections as secondary stress-test cohorts.

#### Recommendations

1. Keep the repo adsorption-first.
2. In analyses and write-ups, treat defense absences as weak evidence.
3. Do not describe our approach as generically better than the preprint's approach; describe it as better aligned to
   the repo's task.
4. Prioritize a future branch for higher-resolution RBP FASTA features.

#### Sources

- Preprint: https://www.biorxiv.org/content/10.1101/2025.01.08.631966v1
- Paper repo: https://github.com/mdmparis/antiphage_landscape_2025
- Core 2024 _Escherichia_ paper: https://doi.org/10.1038/s41564-024-01832-5

### 2026-04-05 11:35 UTC: AUTORESEARCH no longer blocks first search on host defense

#### Executive summary

We refined the AUTORESEARCH plan after the antiphage-landscape read: defense stays in the track, but the first honest
search loop should not be blocked on the slowest defense path. The adsorption-first host/phage cache is now sufficient
for the initial baseline, and host defense becomes an additive block rather than the critical-path gate.

#### What changed

- `AR04` now depends on `AR02` instead of `AR03`, which unblocks host-surface work from DefenseFinder runtime.
- `AR07` now defines its first runnable baseline over the adsorption-first minimum cache:
  `host_surface`, `host_typing`, `host_stats`, `phage_projection`, and `phage_stats`.
- `AR03` now says explicitly that defense-feature absences are annotation-limited evidence, not clean biological
  absence.
- `host_defense` remains part of the frozen cache schema from `AR02`, so this is not a defense deletion or a schema
  contraction.

#### Interpretation

This is the smallest plan change that actually changes the AUTORESEARCH critical path. It preserves defense as possible
additive signal while moving the first serious search onto the feature families that are both more adsorption-aligned
and operationally easier to get running quickly.

### 2026-04-05 22:16 UTC: AUTORESEARCH promotion now requires clean-checkout holdout replay, not RunPod inner-val wins

#### Executive summary

We added the project-level rule for AUTORESEARCH promotion: a RunPod search result is only a candidate, not a winner.
The winning artifact must be imported back into the repo, replayed from a clean checkout on the sealed AR01 holdout,
and compared against the current locked production-intent comparator with the same repeated-seed bootstrap decision
path before any promotion claim is allowed.

#### What changed

- The AR08-to-AR09 handoff is now explicit: the exact `train.py` plus RunPod workflow/pod metadata is the raw input to
  replication, not a screenshot of inner-validation metrics.
- Final replay uses all retained AR01 non-holdout rows for fitting and keeps the sealed holdout only for the last
  comparison. That is the correct thin contract after search because the candidate is frozen and no longer being tuned.
- Promotion is predeclared and metric-based:
  - primary metric: holdout ROC-AUC
  - guardrails: no material regression on holdout top-3 hit rate or Brier score
  - failure mode: emit `no_honest_lift` instead of stretching the evidence

#### Interpretation

This raises the bar in the right place. The expensive RunPod loop is allowed to optimize on inner validation, but it
does not get to self-certify promotion. The only evidence that counts for promotion is the clean-checkout,
sealed-holdout replication bundle produced after the search code is frozen and imported back into the repo.

### 2026-04-08 20:55 UTC: AUTORESEARCH reaches parity with TL18 on honest holdout

#### Executive summary

AUTORESEARCH raw-feature baseline achieves 0.810 ROC-AUC on the ST03 holdout vs TL18's 0.823 — the difference falls
inside the 95% bootstrap CI \[0.765, 0.847\]. The SVD bottleneck was the main obstacle. The remaining gap is
attributable to TL18's defense/kmer/preprocessor features that AUTORESEARCH does not yet include. The flawed AR09
comparator arm (different split, different features) has been scrapped for ST03 evaluations.

#### Strategic decision

Retargeted AUTORESEARCH evaluation from the AR01 holdout (74 bacteria, individual-level) to the ST03 holdout
(65 bacteria, cv_group-disjoint). Rationale:

1. ST03 is the same split TL18 was evaluated on — only honest comparison possible.
2. AR01 holdout doesn't group by genomic similarity; ST03 does (at 1e-4 threshold).
3. The AR09 comparator was meaningless: panel-derived V0 metadata features on a different split produced 0.865 AUC
   that could not be compared to anything.

#### What this means for the project

- AUTORESEARCH is now a viable path to beat TL18. The architecture (raw slots → LightGBM) works; we just need to
  close the feature gap.
- The two most obvious feature additions are defense system annotations and phage genome kmer profiles — exactly the
  features TL18 has that AUTORESEARCH lacks.
- Both can be derived from FASTA inputs, keeping the AUTORESEARCH contract intact (no panel-derived metadata).

### 2026-04-08 22:00 UTC: Knowledge consolidation infrastructure — the `/sleeponit` skill

#### Executive summary

We built a structured knowledge consolidation system inspired by sleep-based memory consolidation from cognitive
science. Lab notebook entries (episodic records) are distilled into a unified `knowledge.yml` (semantic knowledge) that
renders to `KNOWLEDGE.md` and loads into Claude's context for all lyzortx work. This addresses the scaling problem:
~9,500 lines of notebooks are too large to load every session, but the knowledge they contain is essential.

#### Strategic rationale

The 2025-2026 agent memory literature (A-Mem, ICLR 2026 MemAgents workshop, Anthropic's context engineering guide)
converges on episodic-to-semantic transformation as the critical underserved gap in agent memory. Key principles
adopted:

1. **Two tiers, not three.** Always-loaded knowledge map + on-demand raw notebooks. No intermediate digest layer needed
   for file-based context systems.
2. **Structured knowledge units over prose summaries.** Each unit has: statement, sources, status, confidence, and
   cross-references. Machine-manipulable YAML source of truth, not markdown prose.
3. **Unified thematic organization.** Knowledge is grouped by concept (features, model, dead ends), not by source track.
   A finding from Track L and Track E that both relate to features belong together.
4. **Forgetting is a feature.** Implementation details recoverable from code, git history, and exact file paths are
   deliberately excluded.
5. **Incremental updates over full rebuilds.** The `diff_knowledge()` function compares old and new models. Re-running
   `/sleeponit` proposes additions/updates/removals rather than regenerating from scratch.

#### Architecture

Follows the established `plan.yml` → `plan_parser.py` → `render_plan.py` → `PLAN.md` pattern:

- `lyzortx/orchestration/knowledge.yml` — YAML source of truth
- `lyzortx/orchestration/knowledge_parser.py` — frozen dataclasses, loader, validator, diff
- `lyzortx/orchestration/render_knowledge.py` — markdown renderer
- `lyzortx/KNOWLEDGE.md` — rendered output, loaded via `@KNOWLEDGE.md` in `lyzortx/CLAUDE.md`
- `.agents/skills/sleeponit/SKILL.md` — user-invocable skill guiding the 3-phase process

#### Data model

`KnowledgeUnit` fields: `id`, `statement`, `sources`, `status` (active/superseded/dead-end), `confidence`
(validated/preliminary), `context`, `relates_to` (cross-references to other unit IDs). Validation catches duplicate IDs,
empty statements, invalid statuses, and broken cross-references. 15 unit tests cover the parser, validator, diff, and
renderer.

#### Future: First consolidation run

Run `/sleeponit` to produce the initial knowledge model from all existing notebooks. This should happen after the
current AUTORESEARCH baseline work stabilizes, so the knowledge model captures the full state.

### 2026-04-08 22:15 UTC: AUTORESEARCH feature search concluded — base slots are the ceiling

#### Executive summary

After ablating defense features (+0.7pp AUC but -4.6pp top-3) and raw tetranucleotide kmers (zero signal), and
cross-referencing the original paper, we conclude that the base AUTORESEARCH configuration (159 raw slot features,
0.810 AUC on ST03 holdout) is at or near the ceiling for the current architecture. The 1.3pp gap to TL18 is not
statistically significant and is more likely architectural than feature-driven.

#### Key findings

1. **Defense features are net-negative for deployment.** +0.7pp AUC but -4.6pp top-3 due to lineage confounding.
   The paper confirms: only 2 defense traits are significant vs 30 adsorption traits.
2. **Raw kmer features are inert.** 256 dims for 96 phages — dimensionality mismatch. TL18's 24-dim SVD is
   denoising, but the paper achieves 86% AUROC without kmer features entirely, so even the SVD form adds little.
3. **The paper's architecture is different from ours.** Per-phage models using bacterial features only, not all-pairs
   LightGBM. This means TL18's remaining advantages over AUTORESEARCH likely come from its feature engineering
   pipeline (pairwise features, isotonic calibration) rather than from features we haven't tried.

#### Strategic implication

Further feature additions within the current train.py-only surface are unlikely to close the gap. The next
improvement step would be architectural: pairwise interaction features, post-hoc calibration, or per-phage
sub-models. These exceed the current AUTORESEARCH search contract and would need a plan update before proceeding.

### 2026-04-09 01:00 UTC: Track APEX launched — targeting 95% AUC via structural RBP features

#### Executive summary

Launched Track APEX (Adsorption-Prediction EXpansion) with 6 tasks targeting 95% AUC and 95% top-3 on ST03
holdout. This is the first track to incorporate protein structure prediction (via PHIStruct/ESMFold) and per-phage
sub-models into the pipeline.

#### Strategic decision

The AUTORESEARCH feature search is concluded at 0.810 AUC on the current all-pairs LightGBM architecture. Rather
than continuing to search for features within the frozen train.py surface, APEX expands the architecture:

1. **Structural RBP embeddings** — the paper identifies RBPs as the #1 phage-side variable, but AUTORESEARCH only
   has binary family presence. PHIStruct provides structure-aware embeddings that capture binding-interface similarity.
2. **Per-phage sub-models** — the paper achieves 86% AUROC with per-phage models. APEX adds phage-specific
   classifiers that break Straboviridae collapse.
3. **Pairwise compatibility** — label-free RBP-receptor structural matching, replacing the leaky label-derived
   pairwise features that were removed from TL18.

#### What this means for the project

- AUTORESEARCH remains the frozen FASTA-only baseline. APEX extends it with structural and architectural innovations.
- The codex-implement.yml CI workflow is disabled for APEX tasks. The orchestrator creates issues; I implement them
  directly.
- Realistic ceiling: 0.90-0.92 AUC, 95% top-3. The 0.95 AUC target is aspirational but motivating.

### 2026-04-12 00:15 CEST: Track GIANTS launched — three-layer biological model informed by Moriniere and GenoPHI

#### Executive summary

Launched Track GIANTS after literature review of Moriniere 2026 (receptor specificity from genomes) and Noonan 2025
(GenoPHI strain-level prediction). Rather than benchmarking against these papers, we integrate their key findings into a
three-layer biological prediction model: depolymerase→capsule compatibility (Gate 1), RBP→OMP receptor compatibility
(Gate 2), and host defense survival (Gate 3). Baseline: AUTORESEARCH all-pairs 0.810 AUC.

#### Strategic decision

APEX Phase 1-2 exhausted the feature-engineering approach: PLM embeddings, physicochemical descriptors, and undirected
cross-terms all proved neutral or harmful on the all-pairs architecture. The literature reveals why:

1. **Receptor specificity is discrete and localized** — short motifs at specific genomic loci, not global protein
   properties. k=5 amino acid k-mers predict receptor class at AUROC 0.99 (Moriniere 2026).
2. **ML pipeline matters more than features** — algorithm + training strategy explains 5-18x more variance than genomic
   representation (Noonan 2025, 13.2M training runs).
3. **Our hosts are diverse clinical isolates** with 99 capsule features and 12 OMP receptor variants — the signal is
   available if features are directed correctly.

The three-layer architecture matches the actual infection mechanism instead of treating all features as flat inputs.
Defense features (Gate 3) are unrestricted (all known systems, not just CRISPR/RM/Abi) and gated on adsorption success.

#### Baseline policy

AUTORESEARCH all-pairs (0.810 AUC, 90.8% top-3, 0.167 Brier) is the single canonical baseline. TL18 (0.823) has
feature integrity issues (DefenseFinder version drift, soft-leaky pairwise features). Per-phage blend (0.830) is not
deployable and not a valid comparison target.

### 2026-04-12 23:05 CEST: GIANTS → SPANDEX transition

#### Executive summary

Track GIANTS hit a feature-engineering ceiling (7 independent null results, GT04-GT08). The biggest remaining gains
came from label quality (+3.1pp top-3 from excluding ambiguous 'n' pairs, GT09) and evaluation methodology (top-3
discards ranking and potency information). SPANDEX shifts focus from features to data quality, evaluation rigor, and
panel expansion. GT09 cancelled and folded into SX02/SX03. Details in track_SPANDEX.md.

### 2026-04-12 23:05 CEST: Track SPANDEX — evaluation overhaul and panel expansion

#### Executive summary

Track GIANTS established the 0.823 AUC ceiling through 7 independent null results (GT04-GT08) and identified label
quality as the binding constraint (+3.1pp top-3 from excluding ambiguous 'n' pairs). SPANDEX addresses the three
remaining levers: evaluation methodology, label quality, and panel expansion. GT09 (BASEL panel expansion) is cancelled
and folded into SPANDEX.

#### Strategic decisions

1. **Top-3 metric retired.** Top-3 collapses the full phage ranking into a binary "lucky in the first 3 slots" signal.
   It cannot distinguish a model that ranks a true positive 4th from one that ranks it 148th, and it cannot handle
   partial ground truth from multi-source panels. Replaced by nDCG (graded relevance) + mAP (binary retrieval quality).

2. **Graded evaluation using MLC 0-4 dilution potency.** The interaction matrix encodes 5-level lysis strength, not
   just binary. MLC=4 (lysis at 10,000x dilution) is clinically superior to MLC=1 (barely lyses at neat concentration).
   nDCG with graded relevance rewards ranking potent phages higher. This was always available in our data but discarded
   by the binary any_lysis label policy.

3. **k-fold CV replaces fixed ST03 holdout.** 10-fold bacteria-stratified CV gives robust performance estimates not
   dependent on which 65 bacteria happen to be in one holdout split. ST03 kept as a single-fold comparison column for
   backwards compatibility.

4. **BASEL integration done properly.** Track K TK02 was invalidated — zero BASEL rows joined training because features
   were never computed for BASEL phages. SPANDEX runs the full annotation pipeline (Pharokka + DepoScope + feature
   slots) on 52 BASEL genomes before integrating. BASEL data is binary only (single-concentration spot test at >10^9
   pfu/ml), mapping to relevance=1 for positives. No graded upstream data exists.

5. **Pre-flight gates on every ticket.** Each SX ticket has a variance/separability check that runs before the main
   work. SX01 checks whether existing predictions separate MLC grades (if not, graded nDCG is cosmetic and SX04 is
   dead). SX02 checks BASEL phage diversity and depolymerase coverage. SX03 checks feature-space overlap between BASEL
   and Guelin phages.

#### Why SPANDEX, not more GIANTS tickets

GIANTS proved that the feature-engineering approach has diminishing returns within the current evaluation framework.
The next gains come from: (a) measuring more carefully (graded metrics, k-fold CV), (b) cleaning the training signal
(excluding ambiguous labels), and (c) expanding the phage panel (BASEL). These are evaluation and data changes, not
feature changes — a different track with a different philosophy.

#### Future: Two-tower / set-aware phage encoder track

LightGBM requires fixed-width feature vectors, forcing mean-pooling across variable-length sequences (residues within a
protein, proteins within a phage). This dilutes position-specific binding signal (RBP tip motifs) and functional
identity (depolymerase vs tail fiber). The cluster-membership depolymerase features have the same generalization
flaw — they're essentially memorization of training phage sequences. A set-aware architecture (Set Transformer, Deep
Sets, two-tower with cross-attention over phage proteins) would handle variable-length phage inputs natively, with
learned pooling conditioned on the host features. PHIStruct (Kuchi et al. 2024) validates this direction for
phage-host prediction.

Revisit when SPANDEX completes and a ceiling is confirmed that's attributable to pooling/cluster-membership artifacts
(e.g., SX05 closes the zero-filled BASEL gap but generalization AUC still lags within-panel by >5pp, or SX06 shows
that continuous-but-pooled depo features don't rescue k-fold performance). At that point, scope a dedicated track for
a two-tower neural encoder feeding learned phage embeddings into LightGBM (hybrid) or a full end-to-end neural
predictor.

### 2026-04-13 16:00 CEST: SPANDEX plan reshape — MLC=4 pipeline fix elevated to SX05

#### Executive summary

Re-reading the paper Methods revealed that our pipeline's `DILUTION_WEIGHT_MAP` assigns MLC=4 from a dilution the paper
explicitly excluded (`log_dilution=-4`, 5×10⁴ pfu/ml, unreplicated). The MLC=4 pipeline fix was promoted to the first
unimplemented SPANDEX ticket (SX05); the old SX09 empirical MLC sensitivity study was deleted as superseded. Later
tickets shifted by one. Full technical rationale in track_SPANDEX.md (2026-04-13 15:55 CEST entry).

#### Strategic points

1. **MLC=1 is not suspect.** Earlier drafts had flagged MLC=1 for LFW/Abi concerns. Re-reading the paper Methods shows
   the LFW caveat is a general biological note about lawn-clearing-without-plaques, not a label quality concern
   specific to MLC=1. MLC=1 stays as a standard low-potency interaction.

2. **MLC=4 is the real anomaly.** The paper's MLC=4 is a morphological distinction at 5×10⁶ (entire lawn lysis vs
   countable plaques at the same concentration — biologically suspicious because at MOI ~0.1 you would expect plaques,
   not full clearing). Our binary 0/1/n raw data cannot represent this morphological split, and our repurposed MLC=4
   (derived from the unreplicated 5×10⁴ observation) is a different, noisier label.

3. **Executive decision over sensitivity study.** The pipeline fix is principled (paper protocol), not empirical
   (a comparative study of label cohorts). An empirical study would have been cosmetic — it cannot make the
   morphological split detectable, and the label noise is structural.

4. **Revised trigger references for the Two-tower future note above:** the "SX05" in that revisit trigger (written
   before this reshape) refers to the TL17 BASEL projection — now SX06. The "SX06" refers to continuous depolymerase
   bitscore — now SX08. Future `/replan` runs should resolve trigger conditions against the current plan.yml IDs.

#### Ticket mapping

| Old ID | New ID | Task |
|--------|--------|------|
| — | SX05 | MLC mapping fix (new — elevated this entry) |
| SX05 | SX06 | BASEL TL17 phage_projection |
| SX06 | SX07 | BASEL PLM embeddings (restoration hint pinned to PR 393 / SHA d9717ab) |
| SX07 | SX08 | Continuous depolymerase bitscore |
| SX08 | SX09 | Per-functional-class PLM blocks |
| SX09 | — | Deleted (MLC sensitivity study superseded by SX05) |
| SX10 | SX10 | Final consolidation (deps updated, MLC axis removed) |

### 2026-04-14 01:37 CEST: SPANDEX wave-2 plan (SX11–SX15)

#### Executive summary

Five new SPANDEX tickets planned for the second wave. Scope was explicitly narrowed down from ~15 candidates to 5
after walking through the twelve post-SPANDEX shortcomings identified during the closing review: the three
highest-EV structural attacks plus one consolidation ticket and one evaluation-framework ticket. Vision-based
label re-reading was cut because prior spot checks already showed poor image quality (see
`label-vision-reading-spot-checked-dead` knowledge unit). Six other candidate feature families (CRISPR spacers,
phage anti-defense genes, Kaptive K-typing, phylogroup residualization, within-family HVS detection,
structure-based tip features) were deferred to keep wave-2 focused.

#### Strategic framing

Wave 1 (SX05–SX10) fixed label semantics (SX05 MLC 0–3) and reclaimed the BASEL zero-fill gap (SX06 real TL17
features). It left three structural shortcomings untouched:

1. **Binary training target** can't learn potency — the single biggest structural lever I had flagged as
   unexploited at wave-1 close.
2. **OMP homogeneity** collapses every receptor × host cross-term we've tried.
3. **Whole-protein phage features** average over the 5-mer motifs where receptor specificity actually lives.

Wave 2 attacks all three with one ticket each (SX11, SX13, SX12 respectively), plus one for stratified
evaluation honesty (SX14) and one for unified-evaluation-framework simplification (SX15). Everything else
recognized as important is tracked in knowledge / open questions but deferred out of wave-2 scope.

#### Ticket lineup (see `track_SPANDEX.md` for full details)

| ID | Topic | Key hypothesis |
|----|-------|----------------|
| SX11 | Potency loss-function ablation (hurdle / LambdaRank / ordinal all-threshold) | SX04's null was loss-choice-specific, not a verdict on the concept |
| SX12 | Moriniere 5-mer phage features (direct) | GT06's intermediate-classifier path discarded signal; direct-feature use is the honest test |
| SX13 | OMP k-mer host features + SX12 × SX13 cross-term | Mirroring Moriniere onto the host side recovers allelic variation the HMM score averages away |
| SX14 | Wave-2 consolidation + stratified evaluation | Honest reporting by within-family / cross-family / narrow-host / phylogroup-orphan subsets |
| SX15 | Unified Guelin+BASEL k-fold framework (bacteria + phage axes) | Replaces the three-arm SX01/SX03 maze with one honest number; tests deployability on genuinely unseen phages |

#### Cuts and why

- **Vision label re-read (was SX12 in an earlier draft):** prior spot checks showed the plaque images are too
  noisy for a vision model to confidently disambiguate the `n` scores. Dead-end recorded in knowledge.
- **CRISPR spacer pairwise match, phage anti-defense, Kaptive K-typing, phylogroup residualization:** all
  scientifically sound, all deferred. If wave-2 leaves the 10.9 pp cross-panel AUC gap intact, these are the
  first candidates for wave 3.
- **Two-tower / GNN architectures:** out of scope per wave-1 "keep LightGBM" scope discipline.
- **Within-family HVS detection, structure-based tip features, reference-DB OMP allele typing:** higher
  engineering cost for uncertain marginal gain relative to the k-mer approach (SX12 + SX13).

#### Execution plan

- Commit the plan, tick the orchestrator, then implement SX11–SX15 sequentially as PRs with local self-review
  subagents, same workflow as SX05–SX10.
- At wave close (after SX15), produce a fundamentals-style science/biology review of results with case-by-case
  explanations, mirroring the post-SPANDEX wave-1 review format.

### 2026-04-15 22:17 CEST: SPANDEX wave-2 closing review — fundamentals

#### Executive summary

Wave-2 attacked three structural shortcomings identified at SPANDEX wave-1 close: (1) binary training target
not learning potency, (2) omp-score-homogeneity collapsing receptor × host cross-terms, (3) whole-protein
phage features averaging over binding-interface motifs. All three attacks **failed the +2 pp aggregate
acceptance gate**. But SX14's stratified evaluation layer changed the reading: **SX11 ordinal/LambdaRank
losses deliver +2.7–3.5 pp within-family nDCG (LambdaRank CI disjoint from baseline)**, a real stratum-
specific win hidden in aggregate because cross-family pairs (69% of panel) dilute to zero. SX15 then gave
the first honest deployability estimate for unseen phages (AUC 0.8988, nDCG 0.7229) and confirmed BASEL
phages generalize as well as Guelin on phage-axis. The canonical production baseline (`spandex-final-baseline`
= SX10) is unchanged; the wave's real deliverables are three new knowledge artifacts (`spandex-wave-2-baseline`,
`stratified-eval-framework`, `spandex-unified-kfold-baseline`) and the `case-by-case` skill.

#### The wave in plain terms

We set out to fix three things we thought were wrong with SX10:

1. **Potency blindness.** SX10 treats MLC=1 and MLC=3 as identical positive labels. A phage that kills at 5×10⁸
   pfu/ml (weak) looks the same as one that kills at 5×10⁶ (strong). Surely explicit potency supervision helps?

2. **OMP homogeneity.** Our 369 clinical *E. coli* all express the 12 core OMPs at near-identical HMM bit
   scores (CV 0.01-0.17). Every cross-term like `phage_is_OmpC × host_OmpC_score` collapses to phage-only.
   But OmpC has ~50 allelic variants across clinical strains — that variation must matter somewhere?

3. **Whole-protein phage features.** The TL17 phage_projection averages over whole proteomes. Moriniere 2026
   showed receptor class is determined by 5-mer motifs at RBP tip regions. Surely feeding the model those
   specific motifs helps over the averaged features?

For each, we built the test and ran it honestly. All three came back null on aggregate. But stratified
reporting showed that the first hypothesis wasn't fully wrong — it just helps where the data supports it
(within-family).

#### Ticket-by-ticket fundamentals review

**SX11 — Potency loss-function ablation**

**The hypothesis:** Four loss formulations might recover the potency signal the binary label discards.
SX04 had tested vanilla MLC regression on 79%-zero-inflated data and it failed. We tested three losses
that explicitly handle zero-inflation + ordinal structure:

- *Hurdle (two-stage):* `P(y>=1) × E[MLC | y>=1]`. Decouples lysis/no-lysis from potency grade. Canonical
  statistical fix for zero-inflated ordinals.
- *LambdaRank:* directly optimizes nDCG at training time, query-grouped by bacterium. If ranking is what
  we want, train for ranking.
- *Ordinal all-threshold:* three binary classifiers for y>=1, y>=2, y>=3 combined via cumulative
  probabilities. Respects the ladder explicitly.

**The aggregate verdict:** validated null. Best arm (ordinal) was +1.33 pp nDCG over SX11's own binary
baseline, which was still below the +2 pp gate. But the aggregate hid what was actually happening.

**What SX14 stratified eval revealed:** **Within-family (≥3 training-positive pairs for this bacterium's
cv_group), all three alternative losses deliver +2.7–3.5 pp nDCG over SX10.** LambdaRank's within-family
CI [0.838, 0.874] is **disjoint** from SX10's [0.800, 0.840] — the strongest stratum-level signal in the
entire wave. Cross-family (69% of pairs) shows no improvement for any arm, which dilutes to zero
aggregate.

**Why this is biologically sensible:** When the model has seen this phage family hitting hosts of this
cv_group before, it has learned the family's host-range signature. An explicit ordinal loss lets it
reshape the per-bacterium potency gradient among positives — "rank MLC=3 above MLC=1 within this
bacterium's lysers." This costs ~1 pp mAP (some lysis/no-lysis boundary softening) but gains 2.7-3.5 pp
nDCG within-family. In cross-family cold-start, the model has no family-specific signal to refine, so
loss choice is irrelevant — and it shows.

**What it means for production:** SX11 alternative losses are not yet adoption-ready because SX11's arms
ran *without* per-phage blending (the SX10 architectural gain). Integrating ordinal/LambdaRank with
per-phage blending is the natural wave-3 follow-up; if the within-family gain persists, stratum-aware
inference routing becomes worth building.

**Case-by-case notes from SX11 analysis:** NILS20 (3 positives, MLC 3/1/1) is the clearest positive.
Baseline ranks LF82_P8 (MLC=3) at rank 4 behind 3 false positives; ordinal promotes to rank 1. But
H1-006-0003-S-L (2 positives, both MLC=1) is the clearest negative — binary ranks them #1 and #2, ordinal
demotes to #4 and #5 because MLC=0 phages get boosted by their high P(y>=2)/P(y>=3). The trade-off is
mechanical: equal-weighted threshold combination compresses the MLC=1 vs MLC=0 gap.

**SX12 — Moriniere phage 5-mer features (direct)**

**The hypothesis:** Moriniere 2026's 815 receptor-predictive k-mers achieve AUROC 0.99 for receptor class
prediction on K-12 derivatives. GT06 tested them as intermediate-classifier features and saw zero lift.
But maybe feeding them directly as phage features lets LightGBM find the signal GT06's intermediate
classification discarded.

**The verdict:** null at every granularity. Aggregate +0.23 pp AUC (CIs overlap massively). SX14
stratified showed no stratum-specific lift either (all within ±0.2 pp). RFE keeps 95 of 815 k-mers
(11.7%) but they contribute ~5% of total feature importance and **zero appear in the top-20 features**.

**Why it's biologically null in our panel:**

1. **The k-mers are information-redundant with `phage_projection`.** Both encode phage sequence
   similarity — k-mers via shared 5-mer presence, TL17 BLAST via shared protein hits. Two phages sharing
   80% of Moriniere k-mers almost always share TL17 BLAST hits. RFE keeps both but the trees preferentially
   split on the already-engineered `phage_projection` features because they have higher information gain
   per split (top feature `host_serotype` imp=2076 vs top k-mer `NVSVG` imp=11, a 190× gap).

2. **Moriniere selected k-mers for the wrong prediction task.** They discriminate receptor class on
   K-12/BL21 (which lack capsule and O-antigen). In those hosts, receptor identity IS the dominant
   factor. In our clinical *E. coli*, **the bottleneck sits upstream of receptor recognition** —
   polysaccharide barriers gate physical access to OMP receptors. That's what `depo × capsule` features
   (Gate 1 in `three-layer-hypothesis`) capture, and it's the validated +1.2 pp mechanism in GT03.
   Knowing the phage's target receptor doesn't help when the surface is blocked.

3. **Per `same-receptor-uncorrelated-hosts`:** phages sharing a predicted receptor have **weakly
   correlated host ranges** (Tsx phages Jaccard 0.091, below the random-pair baseline of 0.17). Receptor
   identity tells you which OMP a phage *could* bind; it doesn't tell you which strains will actually
   support lysis. That extra information is what we don't have.

**Case-by-case:** ~10 bacteria win meaningfully (IAI78 +0.12 nDCG, ECOR-25 +0.17, NILS31 +0.16), exactly
offset by ~10 losers (ECOR-19 -0.36, EDL933 -0.24, NILS38 -0.13). RFE's effective feature budget is
fixed, so adding 95 correlated k-mers knocks out other useful features in some folds. NILS53 — the
canonical narrow-host failure — gains only +0.011 nDCG. K-mers do not break narrow-host prior collapse.

**SX13 — Host OMP 5-mer features (four-arm ablation)**

**The hypothesis:** Invert the Moriniere trick onto the host side. OmpC has 50 allelic variants across
our clinical panel; surely the extracellular loop variations encoded in 5-mer presence-absence predict
which strains a given phage can actually infect.

**The verdict:** null across all four arms (marginal, cross-term, path1 cluster, baseline). All deltas
within ±0.4 pp of SX10. Permutation test on cross-term aggregate: 73% of random prediction swaps produce
deltas this extreme or larger. Signal indistinguishable from noise.

**What's biologically surprising:** host OMP variation IS substantial. Our 369 clinical hosts span 28
MMseqs2 clusters on BTUB at 99% identity, 49 on OMPC (matching the knowledge model's 50 OmpC variants
claim), 16 on FHUA, 32 on NFRA. The 5% – 95% frequency band retains **5,546 informative k-mers** across
the 12 core OMPs (BTUB alone contributes 1,495). The variation is real and recoverable.

**But variation ≠ prediction.** None of the 4 arms could extract lysis-predictive signal from that
variation. The `omp-score-homogeneity` unit was pointing at the right *biological* bottleneck but at the
wrong *representational* level — HMM scores ARE homogeneous, but escalating to finer representations
doesn't rescue prediction because **host-range variance lives downstream of OMP recognition** in clinical
strains. Polysaccharide access (Gate 1, captured by `depo × capsule`) + post-adsorption factors
(intracellular defenses, injection efficiency, co-evolutionary dynamics) do the real work. OMP-level
features are looking where the light is, not where the key was dropped.

**Case-by-case notes:**

- NILS53 improves +2.59 pp nDCG under cross_term, BUT 51 bacteria at the same lysis rate (12% ±3pp) had
  mean Δ=-0.001. NILS53 is at the 76th percentile of its peer group — **one outlier draw, not a
  reproducible narrow-host rescue**. The `case-by-case` skill's peer-percentile check caught this.
- Small directional signal in marginal arm for moderate-narrow deciles (lysis 5-11%, mean +1-2 pp nDCG).
  Biologically consistent with OMP variation mattering more for specialists, but not significant by sign
  test (75/137 positive, p=0.31).

**SX14 — Wave-2 consolidation + stratified evaluation**

This ticket was the structural deliverable that transformed the wave-2 reading. Per spec, arms that
failed their acceptance gate are excluded from the consolidated model; since SX11/12/13 all failed, the
"consolidated wave-2 final" is literally SX10 unchanged. The stratified eval layer is what made the wave
interpretable.

**The insight:** aggregate metrics on a 33k-pair evaluation can hide real stratum-specific effects. SX11
LambdaRank's +3.48 pp within-family is real; it's invisible in aggregate because cross-family pairs
(69% weight) pull the weighted mean to zero. This generalizes: **every future ticket should report
stratified metrics alongside aggregate**. That's now codified in `stratified-eval-framework` knowledge
and in the `case-by-case` skill.

**SX15 — Unified Guelin+BASEL k-fold**

This ticket gave the first honest deployability estimate on unseen phages. Bacteria-axis (SX10
configuration with BASEL training pairs added to shared ECOR hosts) is essentially identical to SX10 —
adding BASEL's 1240 pairs buys no aggregate lift (consistent with `external-data-neutral`). Phage-axis
(unseen phages, all-pairs model only) gives AUC 0.8988 — *higher* than within-panel — but nDCG 0.7229,
7 pp lower than bacteria-axis.

**The AUC-vs-nDCG divergence is structural:** AUC is threshold-independent and ranks pairs globally; the
all-pairs model's features preserve lysis/no-lysis discrimination cleanly. nDCG and mAP are per-bacterium
ranking metrics; when 15 of each bacterium's 148 predictions are for held-out phages with no per-phage
model, their exact ordering is noisier even though the threshold separation is clean. **~7 pp nDCG is
the per-phage blending tax — the cost of the cold-start-phage scenario.**

**The reassuring finding:** BASEL phages generalize as well as Guelin phages on phage-axis — AUCs
essentially identical (holdout_phage_basel 0.8965 vs holdout_phage_guelin 0.8986). Feature engineering
(TL17 phage_projection from SX06) transfers cleanly across phage collections, not just numerically but
across the full stratified panel.

#### The narrow-host specialist problem remains unresolved

The most clinically relevant failure mode — narrow-host phages with <15% panel-wide lysis rate — stays
unresolved through the wave. All 10 SX14-evaluated arms cluster at nDCG 0.67–0.71 on narrow-host pairs.
None breaks past ~0.71. This is the `narrow-host-prior-collapse` knowledge unit: when broad-host
Straboviridae with 60-70% panel-wide lysis rates dominate model rankings, narrow-host specialists
(Dhillonvirus, Kagunavirus, Autographiviridae at 10-25% lysis) get buried below the top-k cutoff. No
feature family tested in wave 2 reaches it. Future tracks must either (a) rethink the ranking architecture
(per-phage-family priors instead of global all-pairs priors) or (b) expand the panel with more narrow-host
positives to teach the model that broad doesn't mean probable.

#### What we actually learned (beyond nulls)

1. **The binary model's implicit potency signal is stronger than expected.** Spearman 0.25 between
   P(lysis) and MLC among positives — not because we asked, but because higher-MLC pairs have more
   consistent positive training rows across dilutions. Explicit ordinal supervision adds ~3 pp within-family
   nDCG but costs calibration (+2.6 pp Brier for LambdaRank). The binary target captures most of the
   easily-extractable potency signal.

2. **OMP homogeneity was only half-right.** HMM scores are homogeneous. Loop-level k-mer / cluster
   variation is substantial (5546 informative 5-mers, 49 OMPC clusters). But variation doesn't predict
   lysis — the bottleneck is upstream (polysaccharide access) or downstream (post-adsorption), not at
   the receptor level.

3. **Cross-source generalization works.** BASEL phages behave like unseen Guelin phages in the phage-axis
   split. The TL17 family projection features are source-agnostic. This is encouraging for expanding the
   phage panel with other collections (genoPHI, VHRdb, KlebPhaCol as non-*E. coli* benchmarks).

4. **Per-phage blending is load-bearing.** SX10's +2.0 pp aggregate AUC over the all-pairs baseline comes
   from per-phage blending (AX02). It's also 7 pp of nDCG on the phage-axis eval — exactly the gap
   between 0.80 (seen phage) and 0.72 (unseen phage). Per-phage blending works when the phage is in
   training; it cannot be applied to new phages at inference time. That's the deployability bound.

5. **Aggregate-only reporting is dangerous.** If we'd stopped at SX11's aggregate null, we'd have
   wrongly concluded loss-function alternatives offer nothing. Stratified eval changed that. Future
   tracks inherit this as a default.

#### Production baseline: unchanged, but better characterized

The canonical production configuration remains **`spandex-final-baseline` (SX10)** — GT03 all_gates_rfe +
AX02 per-phage blending on SX05-corrected MLC labels with SX06 real TL17 BASEL features. Within-panel:
nDCG 0.7958, AUC 0.8699. Wave-2 did not displace it. But wave-2 added four things the baseline didn't
have before:

1. **Stratified metric decomposition** (within_family / cross_family / narrow_host_phage /
   phylogroup_orphan) as the default reporting format.
2. **Unified Guelin+BASEL k-fold** on two axes, with phage-axis providing the honest unseen-phage
   deployability number (`spandex-unified-kfold-baseline`).
3. **Quantified per-phage blending tax** (~7 pp nDCG when phage is unseen). Important for production
   cost/benefit thinking.
4. **`case-by-case` skill** — automated per-bacterium audit with permutation test, peer-percentile
   outlier detection, and Decision headline classification.

#### What wave-2 didn't address (open questions)

- **Narrow-host ranking.** Persistent failure mode; no arm broke the ~0.71 nDCG ceiling on narrow-host
  phages. `narrow-host-prior-collapse` remains active.
- **~10% ambiguous `n` labels.** Vision re-read was cut; label noise floor remains.
- **Cross-panel cross-family cold-start.** Phage-axis cross-family nDCG 0.67, AUC 0.79 — the honest
  floor for "new phage + new family on new host subpopulation". Panel expansion is the direct path.
- **Deferred feature families:** CRISPR spacer pairwise match, phage anti-defense, Kaptive K-typing,
  phylogroup residualization. Still scientifically sound, explicitly deferred to keep wave-2 scope
  tight. First candidates for wave-3 if within-family SX11 losses don't pan out in per-phage integration.

#### Strategic read on wave cost

Wave-2 produced zero new production configurations and three validated nulls. That's a "success" in
science terms — we learned what doesn't work, tightened our knowledge model, and now have sharper tools
(stratified eval, case-by-case skill, unified k-fold). But it's a "failure" in product terms — four
weeks of engineering for zero prediction-quality lift in production. The trade-off was implicit in the
wave-2 plan's "highest-EV structural attacks" framing; the specific attacks we chose were biologically
motivated but turned out to be dominated by more fundamental bottlenecks (polysaccharide access, panel
size, narrow-host prior collapse) that feature engineering on 96-148 phages can't overcome.

The prescription for future tracks: **lead with stratified-eval-first design.** Don't run aggregate-only
experiments. Treat the +2 pp gate as one signal among several; stratum-specific significance (CI
disjoint) is more meaningful than aggregate point-estimate lift. And accept that the 0.82 AUC ceiling is
probably panel-bound, not feature-bound — closing it means expanding the panel, not adding more features.

### 2026-04-19 07:30 CEST: CHISEL framework pivot — MLC retired, AUC+Brier scorecard, concentration-aware training

#### Executive summary

Track CHISEL replaces SPANDEX as the active modeling track. Three coordinated changes: (1) the
atomic training unit switches from pair-level `any_lysis` rollup to the raw
`(bacterium, phage, concentration, replicate) → {0, 1}` row, with concentration entering the model
as a plain numeric feature; (2) MLC is retired as a training and evaluation label — it is a derived
rule-based rollup on top of the raw observations and forced a metric change at SX05 that we no
longer need; (3) the scorecard collapses to AUC and Brier only. Ranking metrics (top-3, nDCG, mAP)
are dropped because ranking is a product-layer concern rather than a property of a pairwise
biological model. This is a framework pivot on top of the earlier top-3 → nDCG pivot; both were
driven by separation-of-concerns reasoning but only CHISEL commits to it fully.

#### Why MLC is retired as a training/evaluation label

MLC is a derived rollup with two problems. First, its definition depends on the measurement
panel — the Guelin protocol records three replicated dilutions so MLC runs 0–3 (post-SX05 cap)
while BASEL has a single high-titer spot that maps to "≥1 or 0", so MLC arithmetic is panel-
specific. Second, every integration decision (BASEL+→MLC=1/2/3 in SX15) became a label-convention
knob rather than a data-structure choice. Training on the raw observations sidesteps both: each
source contributes the rows it actually produced. MLC survives only as a reporting shorthand in
glossaries and lab notebooks, not as a field in training data.

#### Why AUC+Brier replace ranking metrics

The biological question the model can answer is "given this bacterium, this phage, and this
applied concentration, what is P(lysis)?" That's a discrimination question (AUC) plus a
calibration question (Brier). Everything else — how many phages to administer as a cocktail, what
threshold to apply, how to bias toward narrow-host specialists — is a downstream product-layer
decision that depends on operational constraints the biological model has no view of. Ranking
metrics bake product-layer policy into the scientific scorecard and reward the model for optimising
a decision it doesn't own. Dropping them also resolves the BASEL-vs-Guelin nDCG comparison artifact
(fewer nDCG rungs on binary BASEL labels → cosmetically different nDCG numbers for equal
discrimination).

#### Why (pair, concentration, binary) is the training unit

Three reasons. One, it's the data's own atomic unit — `raw_interactions.csv` is keyed on
`(bacteria, phage, image, replicate, plate, log_dilution, X, Y, score)`; collapsing it imposes a
convention the raw data does not. Two, it integrates BASEL natively (single concentration = one
row) without a label-convention knob. Three, it gives each Guelin pair up to nine training rows
(3 replicates at 5×10⁸ pfu/ml, 2 at 5×10⁷, 3 at 5×10⁶, 1 at 5×10⁴) instead of one rolled-up
label, which should help calibration. Concentration enters as `log_dilution ∈
{0, -1, -2, -4}` — a plain numeric feature, not a dose-response functional form. Rows with
`score='n'` are dropped from training (not silently treated as negatives). Evaluation is scored at
each pair's highest observed concentration so within-panel and cross-source comparisons share a
comparable operating point.

#### Sequencing (CH02 through CH07)

CH01 is documentation only — this entry, the SPANDEX closure note, glossary deprecations, and
knowledge-unit rewrites. The actual work is seven tickets:

- **CH02**: cv_group leakage fix. `assign_bacteria_folds` currently hashes on bacterium name, not
  `cv_group`; 52 of 301 cv_groups split across folds. Fix and revalidate the SX10 canonical
  configuration under the corrected fold design — expected small shift, not a headline change.
  Establishes the baseline the rest of CHISEL builds on.
- **CH03**: row-expanded training matrix with `any_lysis` semantics preserved — a regression
  safety-net. Should give |ΔAUC| < 0.005 vs CH02 (same labels, same features, same model;
  different row layout).
- **CH04**: the behavioral flip. Per-row binary labels, concentration as a numeric feature,
  AUC+Brier scorecard. Establishes `chisel-baseline` and supersedes `spandex-final-baseline`.
- **CH05**: phage-axis + cross-source under CHISEL (SX15 successor). Establishes
  `chisel-unified-kfold-baseline`.
- **CH06**: both-axis 10×10 double-CV holdout — addresses the critic framing that SPANDEX
  closed prematurely on the panel-size-ceiling claim without actually evaluating unseen bacteria
  × unseen phages simultaneously. If AUC crashes to marginal-riding range, that's a new dead-end
  knowledge unit.
- **CH07**: SX12/SX13 feature-family null re-audit under the CHISEL frame. SX11 retires entirely
  because MLC-graded losses no longer apply.

#### The load-bearing argument: separation of concerns

This pivot is not about picking better metrics. It is about drawing a clean line between the
biological model (what the physics of phage infection determines) and the product layer (what to
do with those predictions operationally). SPANDEX blurred them — the scorecard rewarded the model
for making cocktail-composition-aware decisions that belong to a downstream system. CHISEL
commits to predicting calibrated pairwise lysis probabilities and delegating policy to a separate
component. Both the previous top-3 → nDCG pivot and this pivot trace back to that same argument,
and it was implicit in the CHISEL design from the start: training on raw observations, scoring on
AUC+Brier, and leaving ranking to the product layer all express the same decomposition.
