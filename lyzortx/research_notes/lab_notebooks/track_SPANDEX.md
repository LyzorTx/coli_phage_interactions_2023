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
| phage_stats | 96 | 52 | 148 | GC, length, N50, record count |
| phage_projection | 96 | 52 | 148 | BASEL rows zero-filled for RBP family features (no TL17 BLAST DB) |

Limitation: phage_projection features for BASEL are zero-filled because the TL17 RBP family BLAST database was
built from Guelin phages only. BASEL phages would need BLAST against this DB to get non-zero RBP family
memberships. For SX03 integration, this means BASEL phages will rely on depo×capsule cross-terms and phage_stats
but not phage_projection RBP family features.

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
