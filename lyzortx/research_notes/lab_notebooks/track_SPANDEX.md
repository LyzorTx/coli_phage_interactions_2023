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
