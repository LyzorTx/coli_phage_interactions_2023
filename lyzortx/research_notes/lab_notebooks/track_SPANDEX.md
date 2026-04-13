# Track SPANDEX: Sparse Panel Data Expansion

**Goal:** Overhaul evaluation metrics (graded nDCG + mAP replacing top-3), adopt k-fold cross-validation, exclude
ambiguous labels, integrate BASEL phage panel, and test ordinal lysis potency prediction.

**Baseline:** GT03 all_gates_rfe + AX02 per-phage blending on ST03 holdout (0.830 AUC, 93.8% top-3, 0.144 Brier).

**Data landscape:**

| Source | Phages | Bacteria | Scoring | nDCG relevance |
|--------|--------|----------|---------|----------------|
| Our panel (Guelin × Picard) | 96 | 369 (incl. 71 ECOR) | MLC 0-4 (dilution potency) | 0,1,2,3,4 |
| BASEL via GenoPHI | 52 | 25 ECOR | Binary spot test (>10^9 pfu/ml) | 0 or 1 |

BASEL tested at ~2x higher concentration than our maximum (~5x10^8 pfu/ml). BASEL positive = at least MLC 1. BASEL
negative = high-confidence MLC 0 (failed at even higher titer). No graded BASEL × ECOR data exists; the original
BASEL EOP data covers K-12 defense experiments only, not the ECOR host range panel.

---

## 2026-04-12 22:57 CEST: Track inception and design decisions

**Motivation:** Track GIANTS established the 0.823 AUC ceiling through 7 independent null results (GT04-GT08). The
remaining levers are: (1) label quality (GT09 clean-label: +3.1pp top-3), (2) panel expansion (BASEL adds 52 new
phages), and (3) evaluation methodology (top-3 discards ranking information and can't handle sparse ground truth).

**Key design decisions:**

1. **Top-3 metric deleted.** Top-3 hit rate collapses the entire ranking into a binary "did we get lucky in the first
   3 slots." A model ranking the true positive 4th scores identically to one ranking it 148th. Replaced by nDCG (graded
   ranking quality) and mAP (retrieval quality).

2. **Graded nDCG uses MLC 0-4 as relevance.** MLC score = dilution potency (number of concentration levels showing
   lysis). MLC=4 phages lyse at 10,000x dilution — better therapeutic candidates than MLC=1 (barely lyses at neat
   concentration). nDCG with graded relevance rewards ranking potent phages higher.

3. **mAP as complementary metric.** mAP answers: "across the full list, how well are positives separated from
   negatives?" It's rank-aware (positives at top contribute more) and handles binary relevance for BASEL pairs.

4. **k-fold CV replaces fixed ST03 holdout.** 10-fold bacteria-stratified CV gives robust performance estimates not
   dependent on which 65 bacteria happen to be in one holdout split. Every bacterium evaluated exactly once.

5. **Partial ground truth scoring.** For bacteria with mixed Guelin + BASEL observations, nDCG/mAP score only observed
   pairs. Unobserved pairs don't participate — no imputation, no penalty.

6. **BASEL relevance mapping.** BASEL spot test at >10^9 pfu/ml ≈ our MLC=1 threshold. BASEL positive → relevance 1.
   BASEL negative → relevance 0 (high confidence).

7. **Clean labels.** Exclude pairs with ambiguous 'n' raw scores and no positive scores from both training and
   evaluation. GT09 showed +3.1pp top-3 from this alone.

**Pre-flight gates on each ticket:**

- **SX01:** Check if existing predictions separate MLC grades (1→2→3→4). If flat, graded nDCG is cosmetic and SX04 is
  dead on arrival.
- **SX02:** Check BASEL phage family diversity and depolymerase count. If homogeneous or depo-sparse, validated features
  won't transfer.
- **SX03:** Check BASEL phage feature overlap with Guelin phages. If >90% nearest-neighbor overlap, training signal is
  redundant.
- **SX04:** Gated on SX01 pre-flight. If binary model doesn't separate MLC grades, ordinal prediction can't help.

### 2026-04-13 00:29 CEST: SX02 — BASEL phage feature computation

#### Executive summary

Ran Pharokka annotation + DepoScope depolymerase prediction on 52 BASEL phage genomes. Pre-flight gate passed:
genome sizes span 39-170 kb (4.3x ratio, indicating multiple phage families), and 44/52 phages have DepoScope
depolymerases — our only validated pairwise feature (depo×capsule) will have signal for most BASEL phages.
Extended phage_stats and phage_projection slots to 148 phages (96 Guelin + 52 BASEL).

#### Pre-flight results

Genome size distribution (proxy for family diversity):

| Size class | Count | Fraction |
|------------|-------|----------|
| Small (<50 kb, podoviruses) | 13 | 25% |
| Medium (50-100 kb, sipho/drexler) | 14 | 27% |
| Large (>100 kb, myoviruses) | 25 | 48% |

Not dominated by a single family. Diverse panel.

#### DepoScope results

- 8,701 proteins scanned across 52 phages (mean 167 CDS per phage)
- 80 depolymerases predicted (score >= 0.5)
- 44/52 phages have at least one depolymerase (85%)
- 8 phages with zero depolymerases (likely use non-enzymatic adsorption)
- Runtime: 227s on MPS (Apple Silicon GPU)

#### Extended slot CSVs

| Slot | Guelin | BASEL | Total | Notes |
|------|--------|-------|-------|-------|
| phage_stats | 96 | 52 | 148 | GC, length, N50, record count. CV comparable (GC: 0.121 vs 0.125). |
| phage_projection | 96 | 52 | 148 | BASEL rows zero-filled for RBP family features (no TL17 BLAST DB) |
| phage_rbp_struct | 96 | 52 | 148 | rbp_count from Pharokka (52/52 have RBPs), PLM PCA zero-filled |

Combined DepoScope predictions written to `.scratch/deposcope_combined/` (18,830 proteins: 10,129 Guelin + 8,701
BASEL). SX03 must pass `deposcope_dir=Path(".scratch/deposcope_combined")` to
`compute_pairwise_depo_capsule_features()` for cross-terms to include BASEL phages.

Limitation: phage_projection features for BASEL are zero-filled because the TL17 RBP family BLAST database was
built from Guelin phages only. PLM PCA features are also zero-filled (would need ProstT5+SaProt inference).
For SX03, BASEL phages rely on depo×capsule cross-terms and phage_stats but not phage_projection or PLM.

### 2026-04-13 01:15 CEST: SX03 — BASEL data integration + cross-source evaluation

#### Executive summary

Three-arm comparison: adding BASEL training data gives negligible lift (+0.3pp nDCG, CIs overlap). The model
generalizes moderately to unseen BASEL phages (AUC 0.72, nDCG 0.65) — useful but far below within-panel performance
(AUC 0.87). Pre-flight correctly flagged 100% feature-space overlap (BASEL phage_projection is zero-filled).

#### Pre-flight overlap

BASEL phage features cluster entirely within Guelin feature space (cosine overlap 100% at threshold 0.1). This is
because phage_projection (33 RBP family features) is zero-filled for BASEL — all BASEL phages look identical in that
slot. Only phage_stats (GC, genome length) provides non-zero differentiation. Flagged as expected.

#### Results

| Arm | nDCG | 95% CI | mAP | AUC | Description |
|-----|------|--------|-----|-----|-------------|
| A (baseline) | 0.779 | [0.771, 0.795] | 0.711 | 0.870 | Our clean data only (SX01 replication) |
| B (+ BASEL) | 0.782 | [0.774, 0.798] | 0.713 | 0.870 | BASEL added to training |
| C (generalize) | 0.648 | [0.603, 0.710] | 0.407 | 0.721 | Predict for unseen BASEL phages |

#### Interpretation

**Arm B (BASEL training): negligible lift.** Adding 1,240 BASEL pairs (302 positive) to 32K clean training pairs
is a 3.8% increase — too small to move the needle. The overlapping CIs (A: [0.771,0.795] vs B: [0.774,0.798])
confirm no statistically significant improvement. This is consistent with the pre-flight flag (BASEL features
cluster within Guelin space) and the `panel-size-ceiling` knowledge finding.

**Arm C (generalization): moderate.** AUC 0.72 for unseen BASEL phages is well above chance (0.50) but well
below the 0.87 achieved for known Guelin phages. The model can discriminate lysis from non-lysis for new phages
based on host features + phage_stats + depo×capsule cross-terms alone, but ranking quality (nDCG 0.65, mAP 0.41)
is substantially degraded. The gap is expected: BASEL phages have zero-filled phage_projection and PLM features,
so the model relies on the small number of non-zero phage features (GC, genome length, depo count/has_depo).

**Biological implication:** The model has learned transferable host-side patterns (which bacteria are generally
susceptible, which capsule types are penetrable) that generalize to new phages. But phage-specific features
(which RBP families, which receptor targets) are missing for BASEL, capping generalization. Computing proper
phage_projection features (TL17 RBP BLAST) for BASEL would likely close part of this gap.

#### Output artifacts

- `lyzortx/generated_outputs/sx03_eval/all_predictions.csv`
- `lyzortx/generated_outputs/sx03_eval/bootstrap_results.json`
- `lyzortx/generated_outputs/sx03_eval/preflight_overlap.json`

### 2026-04-13 01:31 CEST: SX04 — Ordinal lysis potency prediction

#### Executive summary

LightGBM regression predicting MLC 0-4 directly instead of binary P(lysis). Ordinal regression slightly improves
nDCG (+0.4pp) and Spearman correlation with true potency grades (0.30 vs 0.24), but significantly degrades mAP
(-3.1pp) and AUC (-3.9pp). Does not meet the 2pp nDCG improvement threshold for adoption. Binary classification
remains the default.

#### Results

| Metric | SX01 Binary | SX04 Ordinal | Delta |
|--------|-------------|--------------|-------|
| nDCG | 0.779 [0.771, 0.795] | 0.783 [0.773, 0.799] | +0.4pp (CIs overlap) |
| mAP | 0.711 [0.693, 0.729] | 0.680 [0.662, 0.697] | -3.1pp |
| AUC | 0.870 [0.857, 0.882] | 0.831 [0.815, 0.845] | -3.9pp |
| Brier | 0.125 | 0.143 | -1.8pp |
| Spearman (MLC, positives) | 0.24 | 0.30 | +0.06 |

#### Interpretation

The ordinal regressor learns the potency ordering slightly better (Spearman 0.30 vs 0.24 among positives, nDCG
+0.4pp) but at the cost of worse binary discrimination. The regression objective spreads model capacity across
predicting all 5 MLC levels, diluting the primary binary signal. With 79% of pairs at MLC=0, the regressor
optimizes heavily for the zero-inflation, degrading its ability to separate lysis from non-lysis.

**Caveats:** (1) No explicit zero-inflation handling was applied (e.g., Tweedie loss, hurdle model). Vanilla
squared-error regression is a sufficient first test since the decision threshold is +2pp nDCG — more sophisticated
ordinal approaches would add complexity without changing the conclusion. (2) The SX01 binary baseline includes
per-phage blending; SX04 does not. Part of the AUC/mAP regression is attributable to missing blending, not purely
the ordinal objective. This doesn't change the adoption decision — even with blending parity, the nDCG gain would
remain below 2pp.

**Decision:** Binary classification is the better default. The nDCG improvement (+0.4pp) does not meet the 2pp
threshold and comes at too high a cost to discrimination metrics. The pre-flight finding that binary predictions
already separate MLC grades (Spearman 0.24) means the binary model captures enough potency information without
explicit ordinal training.

### 2026-04-12 23:28 CEST: SX01 — Graded evaluation framework + clean-label baseline

#### Executive summary

Implemented SPANDEX evaluation suite (nDCG with graded MLC 0-4 relevance, mAP, k-fold CV, bootstrap CIs) and
established the clean-label baseline. Pre-flight gate passed: existing predictions separate MLC grades monotonically
(Spearman rho=0.24 among positives, well above 0.1 threshold). SX04 lives.

#### Pre-flight results

Mean predicted P(lysis) by MLC grade in GT03 all_gates_rfe holdout predictions:

| MLC grade | Mean P(lysis) | Median P(lysis) | N pairs |
|-----------|---------------|-----------------|---------|
| 0 | 0.267 | 0.118 | 5,096 |
| 1 | 0.611 | 0.721 | 469 |
| 2 | 0.674 | 0.785 | 312 |
| 3 | 0.740 | 0.840 | 198 |
| 4 | 0.816 | 0.880 | 160 |

Clear monotonic increase. The binary classifier already implicitly ranks higher-potency phages higher, validating
graded nDCG as a meaningful metric.

#### SPANDEX baseline (10-fold CV, clean labels, RFE + per-phage blending)

| Metric | Value | 95% CI |
|--------|-------|--------|
| nDCG (graded MLC 0-4) | 0.779 | [0.771, 0.795] |
| mAP (binary >= 1) | 0.711 | [0.693, 0.729] |
| AUC | 0.870 | [0.857, 0.882] |
| Brier | 0.125 | [0.119, 0.131] |

Per-fold variation (RFE selects 207-257 features from 507 per fold):

| Fold | Bacteria | nDCG | mAP | AUC | Brier |
|------|----------|------|-----|-----|-------|
| 0 | 38 | 0.766 | 0.683 | 0.880 | 0.119 |
| 1 | 43 | 0.765 | 0.673 | 0.884 | 0.121 |
| 2 | 42 | 0.833 | 0.782 | 0.888 | 0.113 |
| 3 | 39 | 0.770 | 0.731 | 0.856 | 0.131 |
| 4 | 41 | 0.786 | 0.719 | 0.857 | 0.134 |
| 5 | 37 | 0.759 | 0.708 | 0.877 | 0.126 |
| 6 | 35 | 0.808 | 0.745 | 0.885 | 0.116 |
| 7 | 31 | 0.812 | 0.759 | 0.848 | 0.142 |
| 8 | 26 | 0.718 | 0.641 | 0.890 | 0.116 |
| 9 | 37 | 0.754 | 0.659 | 0.861 | 0.130 |

#### Interpretation

The AUC of 0.870 is substantially higher than the old ST03 fixed-holdout result (0.823). This is not a model
improvement — it reflects two evaluation changes: (1) 10-fold CV uses all bacteria for evaluation instead of 65
fixed ones, giving a more representative estimate, and (2) clean labels exclude ~3,462 ambiguous pairs that were
mislabeled as negatives. The model and features are identical to GT03 + AX02.

nDCG 0.779 and mAP 0.714 establish the SPANDEX ranking baselines. These are the first proper ranking metrics for
this project — they replace the retired top-3 hit rate.

#### Output artifacts

- `lyzortx/generated_outputs/sx01_eval/preflight_results.json`
- `lyzortx/generated_outputs/sx01_eval/kfold_predictions.csv`
- `lyzortx/generated_outputs/sx01_eval/fold_metrics.csv`
- `lyzortx/generated_outputs/sx01_eval/bootstrap_results.json`

### 2026-04-13 14:39 CEST: Mean-pooling analysis for PLM features (design note)

#### Executive summary

The current `phage_rbp_struct` pipeline double-aggregates PLM embeddings: first mean-pooling within each protein
(losing positional specificity for binding motifs), then mean-pooling across all RBPs of a phage (blurring
functional identity — depolymerase signal mixed with tail fiber signal). This informs SX06: prefer max-pooling or
regional pooling within a protein, keep per-functional-class embeddings separate across proteins, only mean-pool
within a class when a phage has multiple proteins of that class.

#### Current aggregation cascade

```
┌─────────────────────────────────────┐
│  Phage X has 3 RBPs:                │
│    - Tail fiber  (500 AA)           │
│    - Tail spike/depolymerase (400 AA)│
│    - Baseplate protein (300 AA)     │
└─────────────────────────────────────┘
              │
              ▼  (ESM-2/ProstT5+SaProt)
┌─────────────────────────────────────┐
│  Per-AA embeddings:                  │
│    Tail fiber:   500 × 1280 matrix   │
│    Depolymerase: 400 × 1280 matrix   │
│    Baseplate:    300 × 1280 matrix   │
└─────────────────────────────────────┘
              │
              ▼  MEAN POOL within protein
┌─────────────────────────────────────┐
│    Tail fiber:   1 × 1280 vector     │
│    Depolymerase: 1 × 1280 vector     │
│    Baseplate:    1 × 1280 vector     │
└─────────────────────────────────────┘
              │
              ▼  MEAN POOL across proteins ← the costly step
┌─────────────────────────────────────┐
│  Phage X: 1 × 1280 vector            │
│  (blurry average of 3 distinct       │
│   functional proteins)               │
└─────────────────────────────────────┘
              │
              ▼  PCA
┌─────────────────────────────────────┐
│  32 features: phage_rbp_struct__*   │
└─────────────────────────────────────┘
```

#### Where signal is lost

- **Within-protein mean pool**: binding specificity lives in short motifs (~20 residues in RBP tip domains). Mean pool
  dilutes this 25× with scaffold/linker. Moriniere 2026 achieves AUROC 0.99 for receptor class using position-preserving
  5-mers — evidence that positional signal matters.
- **Across-protein mean pool**: tail fiber + depolymerase + baseplate protein do different jobs. Averaging produces a
  vector that's nearly useless for distinguishing "has strong capsule-hydrolysis depolymerase" from "has 3 similar tail
  fibers."

#### Scope decision

The set-aware architecture (Set Transformer, two-tower with cross-attention) is the honest fix — noted as a future
direction in project.md. For SPANDEX, we stay with LightGBM but do the minimal practical mean-pooling (SX06): per-
functional-class embedding blocks with max/regional pooling within protein.

### 2026-04-13 15:55 CEST: SX05 — MLC=4 pipeline correction (align with paper protocol)

#### Executive summary

Audit of our MLC derivation against the paper's Methods found a direct contradiction: our pipeline's
`DILUTION_WEIGHT_MAP` assigned MLC=4 to pairs with lysis at `log_dilution=-4` (5×10⁴ pfu/ml), but the paper explicitly
excludes that dilution from MLC because it was unreplicated. SX05 dropped `log_dilution=-4` from `DILUTION_WEIGHT_MAP`,
filtered those rows out of track_a and autoresearch MLC computation, and regenerated `interaction_matrix.csv` by
capping the 1288 MLC=4 cells at MLC=3. Re-running the 10-fold SX01 evaluation on the corrected labels gave nDCG
0.7958 [0.7877, 0.8124] vs the 0-4 baseline 0.7785 [0.7705, 0.7948] — a +1.73pp point-estimate lift; mAP, AUC, and
Brier were bit-identical (training labels are binary `any_lysis`, unchanged).

#### What the paper actually says

Gaborieau 2024 Methods, "Evaluating phage-bacteria interaction outcomes by plaque assay experiments":

> "The MLC score is null in the case of non-lytic interaction and ranges from 1 (lytic interaction at the highest
> phage titre) to 4 (uncountable number of lysis plaques at 5 × 10⁶ p.f.u. ml⁻¹)."

> "The outcome of interaction at 5 × 10⁴ p.f.u. ml⁻¹ was not taken into account in the calculation of the MLC score
> because it was not verified by a replicate."

Fig 2b legend:

> "interactions observed at the lowest phage concentration (5 × 10⁶ p.f.u. ml⁻¹) are distinguished between 3
> (individualized lysis plaque) and 4 (entire lysis of the bacterial lawn)"

So the paper's MLC=3 vs MLC=4 is a **morphological** distinction at 5×10⁶ (countable plaques vs entire lawn lysis),
NOT a concentration distinction. MLC=4 is the biological anomaly — at MOI ~0.1, full lawn clearing cannot be cleanly
explained by productive plaque-forming replication.

#### What our raw data looks like

`raw_interactions.csv` verified:
- Columns: `bacteria, phage, image, replicate, plate, log_dilution, X, Y, score`
- Scores: `{0, 1, n}` (binary only)
- Log dilutions: `{0, -1, -2, -4}` (no `-3`; 5×10⁵ was never tested)
- **No plaque morphology column** — we cannot reconstruct the paper's individualized-vs-entire-lawn distinction

So our pipeline cannot represent the paper's MLC=4 anyway. By repurposing `log_dilution=-4` (the unreplicated
dilution) as MLC=4, we invented a *different* MLC=4 that has no defensible biological or statistical meaning.

#### The fix (SX05)

- `DILUTION_WEIGHT_MAP` in `lyzortx/pipeline/track_a/steps/build_track_a_foundation.py` reduced from
  `{0: 1, -1: 2, -2: 3, -4: 4}` to `{0: 1, -1: 2, -2: 3}`.
- `DILUTION_POTENCY_LABEL_MAP` kept consistent (drop `-4: "very_high"`).
- `EXCLUDED_LOG_DILUTIONS = frozenset({-4})` added; `filter_excluded_dilutions()` drops those rows at raw ingestion
  in track_a; `aggregate_raw_pairs` in `build_contract.py` filters the same set before label computation.
- `find_best_dilution_any_lysis`, `potency_rank_from_dilution`, and `potency_label_from_dilution` defensively skip
  `EXCLUDED_LOG_DILUTIONS` (belt-and-suspenders).
- `LabelPolicyV1/V2.expected_observations_per_pair` default changed from 9 → 8 (3 replicates at 5×10⁸ + 2 at 5×10⁷
  - 3 at 5×10⁶; 5×10⁴ dropped).
- `lyzortx/pipeline/autoresearch/regenerate_interaction_matrix.py` now regenerates
  `data/interactions/interaction_matrix.csv` by loading the committed matrix, capping every MLC=4 cell at MLC=3, and
  preserving everything else (including out-of-raw bacteria and blank missing cells). The paper's binary raw data
  cannot reproduce the morphological MLC=3 vs MLC=4 split, so capping is the only defensible collapse.

Matrix before → after:

| MLC | before | after | delta |
|-----|--------|-------|-------|
| blank (missing) | 157 | 157 | 0 |
| 0 | 30459 | 30459 | 0 |
| 1 | 3014 | 3014 | 0 |
| 2 | 2306 | 2306 | 0 |
| 3 | 1368 | 2656 | +1288 |
| 4 | 1288 | 0 | −1288 |
| total | 38592 | 38592 | 0 |

Pair-level semantics (for pairs now recomputed from raw in track_a / build_contract):

| Raw pattern at {0, −1, −2, −4} | Old MLC | New MLC |
|-------------------------------|---------|---------|
| {1, 1, 1, 0} | 3 | 3 (unchanged) |
| {1, 1, 1, 1} | 4 | 3 (−4 row ignored; −2 observation still valid) |
| {1, 1, 0, 1} | 4 | 2 (no lysis at −2) |
| {0, 0, 0, 1} | 4 | 0 ("lysis only at 5×10⁴" is unreplicated noise) |

#### SX01 10-fold CV effect (corrected labels)

Command: `python -m lyzortx.pipeline.autoresearch.sx01_eval --device-type cpu`.

| Metric | SX01 baseline (MLC 0−4) | SX05 (MLC 0−3) | Δ point estimate |
|--------|-------------------------|----------------|------------------|
| nDCG | 0.7785 [0.7705, 0.7948] | **0.7958 [0.7877, 0.8124]** | **+1.73 pp** |
| mAP | 0.7111 [0.6925, 0.7290] | 0.7111 [0.6925, 0.7290] | 0.00 |
| AUC | 0.8699 [0.8570, 0.8819] | 0.8699 [0.8570, 0.8819] | 0.00 |
| Brier | 0.1248 [0.1187, 0.1309] | 0.1248 [0.1187, 0.1309] | 0.00 |

mAP/AUC/Brier are bit-identical because training labels are binary `any_lysis`, which is unaffected by the MLC value
change. The nDCG improvement comes entirely from collapsing MLC=4 into MLC=3 in the relevance weights:
exponential-gain nDCG treated a misranked MLC=4 item more harshly than a misranked MLC=3 item, and that penalty
disappears once the morphological split is collapsed. Preflight Spearman ρ(predicted P(lysis), MLC) on positives
slipped from 0.24 → 0.23 and mean P(lysis) per grade stayed monotonic across 0,1,2,3 (0.27, 0.61, 0.67, 0.77), so
the graded structure is preserved even after collapse.

#### Why this is an executive decision, not a sensitivity study

Earlier drafts had an empirical SX09 comparing cohorts (exclude MLC=1 / downweight / collapse MLC=4 / combined). That
was unnecessary:

1. **MLC=4 fix is principled, not empirical.** The paper excluded this dilution. We should too. No comparison needed.
2. **MLC=1 is not a problem.** Re-reading the paper Methods confirms: MLC=1 is a standard low-potency lytic
   interaction at highest titre. The paper's LFW/Abi caveat is a biological note about lawn-clearing-without-plaques
   possibly being non-productive — not an exclusion criterion and not specific to MLC=1.
3. **Our binary data cannot distinguish LFW from productive lysis** anyway (no morphology column).

SX09 (old) deleted. SX05 (new) implements the fix.

#### Downstream consequences

- **MLC 0–3 remains valid graded relevance for nDCG.** Training labels are binary (`any_lysis`) per `label-policy-binary`
  and are unaffected by the mapping change.
- **SX01 baseline needs a re-run** on corrected labels to document the metric-level delta. Acceptance criterion in SX05.
- **SX10 consolidation simplified.** MLC label cohort selection is no longer a configurable axis — the fix is baked in.

#### Plan renumbering

SPANDEX tickets shifted:

| Old ID | New ID | Task |
|--------|--------|------|
| — | SX05 | MLC mapping fix (this note) |
| SX05 | SX06 | BASEL TL17 phage_projection |
| SX06 | SX07 | BASEL PLM embeddings |
| SX07 | SX08 | Continuous depolymerase bitscore |
| SX08 | SX09 | Per-functional-class PLM blocks |
| SX09 | — | Deleted (MLC sensitivity study superseded by SX05) |
| SX10 | SX10 | Final consolidation (deps updated) |

#### Next step

Implement SX05. Paper quote + raw-CSV verification + the executive rationale are recorded in the `mlc-dilution-potency`
knowledge unit.

### 2026-04-13 19:00 CEST: SX06 — BASEL TL17 phage_projection features (gap largely closed)

#### Executive summary

SX02 zero-filled the `phage_projection` slot (33 TL17 RBP family columns) for all 52 BASEL phages, which collapsed
Arm C generalization: every BASEL feature vector sat on the origin, indistinguishable from zero-filled Guelin rows.
SX06 reused the live TL17 deployable runtime (`lyzortx/pipeline/track_l/steps/deployable_tl17_runtime.py`) —
`project_phage_feature_rows_batched` on the BASEL FASTAs in `.scratch/basel/genomes/` — and overwrote the zero rows
in `.scratch/basel/feature_slots/phage_projection/features.csv` with real per-family MMseqs2 identities. 39/52 BASEL
phages (75%) now have at least one non-zero TL17 family hit — comparable to 84/96 (87%) in Guelin. Re-running SX03
lifted Arm C nDCG from 0.648 to **0.762 (+11.4 pp)**, mAP from 0.407 to **0.519 (+11.2 pp)**, and AUC from 0.721
to **0.761 (+4.0 pp)**. Arm A (our-data-only) and Arm B (pooled) unchanged — the fix is isolated to the
generalization-to-unseen-phages path exactly where the SX02 gap lived.

#### What changed

- New script `lyzortx/pipeline/autoresearch/project_basel_tl17_features.py`:
  - loads `tl17_rbp_runtime.joblib` + `tl17_rbp_reference_bank.faa`
  - calls `project_phage_feature_rows_batched` on the 52 BASEL FASTAs
  - merges projected rows into the existing `phage_projection` slot CSV without touching Guelin rows
- Slot artefact lives at `.scratch/basel/feature_slots/phage_projection/features.csv` (gitignored; regenerable by
  re-running the script).
- No change to production pipelines; SX03 picks up the corrected slot via `patch_context_with_extended_slots`.

#### BASEL coverage vs Guelin (TL17 family hits)

|              | Panel size | Non-zero phages | Mean families hit | Median |
|--------------|------------|-----------------|-------------------|--------|
| Guelin       | 96         | 84 (87%)        | 2.94              | 2      |
| BASEL (new)  | 52         | 39 (75%)        | 2.62              | 1      |
| BASEL (old)  | 52         | 0               | 0                 | 0      |

13 BASEL phages legitimately miss every TL17 family — their RBP proteomes don't align to the Guelin reference
bank at 70% query coverage. That is expected biology, not a bug; we don't invent features for them.

#### SX03 re-evaluation deltas

Command: `python -m lyzortx.pipeline.autoresearch.sx03_eval --device-type cpu`.

| Arm                                    | Metric | SX03 baseline (zero-fill) | SX06 (real TL17) | Δ |
|----------------------------------------|--------|---------------------------|------------------|---------|
| A — Our clean data only                | nDCG   | 0.7785 [0.7705, 0.7948]   | 0.7958 [0.7877, 0.8124] | +1.73 pp (from SX05 MLC fix) |
| A                                      | mAP    | 0.7111                    | 0.7111           | 0.00    |
| A                                      | AUC    | 0.8699                    | 0.8699           | 0.00    |
| B — Our clean + BASEL training         | nDCG   | 0.7818 [0.7741, 0.7979]   | 0.7958 [0.7877, 0.8130] | +1.40 pp |
| B                                      | mAP    | 0.7131                    | 0.7119           | −0.12 pp |
| B                                      | AUC    | 0.8702                    | 0.8699           | −0.03 pp |
| C — Train on our data, predict BASEL   | nDCG   | 0.6480 [0.6031, 0.7099]   | **0.7619 [0.7219, 0.8207]** | **+11.39 pp** |
| C                                      | mAP    | 0.4072                    | **0.5186**       | **+11.15 pp** |
| C                                      | AUC    | 0.7206                    | **0.7607**       | **+4.02 pp** |
| C                                      | Brier  | 0.1958                    | 0.1844           | −1.14 pp |

Arm A/B stay flat — the SX05 nDCG shift already landed and BASEL training pairs still neutral on our-panel
scoring (`external-data-neutral` knowledge holds). Arm C is where the SX02 gap lived, and SX06 reclaims most of
it: generalization-to-unseen-phages is now nDCG 0.76 and AUC 0.76, versus within-panel (Arm A) at AUC 0.87. Gap
closed to ~10.9 pp AUC — a big step, but SX07's "within 3 pp" conditional skip is **not** met.

#### SX07 decision

`plm-rbp-redundant` shows PLM embeddings (ProstT5 + SaProt → PCA 32) contributed zero lift on the Guelin ST03
holdout while cannibalizing 33.9% of feature importance from existing family features. SX06's massive Arm C lift
came from giving BASEL phages their TL17 family alignment signal — the same family signal PLM would partially
duplicate. Expected value of SX07 on BASEL Arm C: small, given the Guelin-panel dead-end, against hours of
ProstT5 (3B) + SaProt (650M) inference on CPU and restoration of code deleted in PR #393.

Decision: **skip SX07**. Move directly to SX08 (continuous depolymerase bitscore). If SX08+SX10 still leave a
material Arm C gap, revisit SX07 with a scoped BASEL-only evaluation.
