### 2026-03-22: TF01 ST03 v1 benchmark protocol and bootstrap CIs

#### Executive summary

Locked the ST0.3 grouped host split as the canonical v1 benchmark protocol by surfacing the protocol ID in the split
artifacts and carrying it through the TG02 benchmark summary. Added 1000-strain bootstrap confidence intervals for the
v1 benchmark's dual-slice metrics on the isotonic LightGBM outputs: ROC-AUC, Brier score, ECE, and top-3 hit rate.
The benchmark summary now records both full-label and strict-confidence results against the same ST0.3 contract.

#### Interpretation

- The split contract is now explicit rather than implied by the salt string alone, which makes downstream benchmark
  reporting easier to audit and compare across model revisions.
- Bootstrapping at the holdout-strain level keeps the confidence intervals aligned with the evaluation denominator used
  by the recommendation metric, instead of treating pairs as independent samples.
- The strict-confidence slice remains a materially smaller and harder evaluation set, so the dual-slice reporting is
  still necessary to avoid overstating benchmark performance.

### 2026-03-22: TF02 v0 vs v1 holdout comparison and error buckets

#### Executive summary

Compared the v0 metadata logreg baseline from ST0.5/0.7 against the locked v1 genomic GBM winner from TG05
(`defense + OMP + phage-genomic`). v1 improves every tracked holdout metric and reduces the holdout miss count from
10 to 8, but it does not eliminate the hard strains.

#### Headline v0 -> v1 comparison

| Metric | v0 (metadata logreg) | v1 (genomic GBM) | Delta |
|--------|----------------------|------------------|-------|
| ROC-AUC | 0.826948 | 0.910766 | +0.083818 |
| Top-3 hit rate (all strains) | 84.615% | 87.692% | +3.077% |
| Top-3 hit rate (susceptible only) | 87.302% | 90.476% | +3.174% |
| Brier score | 0.171223 | 0.109543 | -0.061680 |

#### Error bucket analysis

- **Fixed by v1:** `ECOR-14`, `NILS41`, `NILS70`, `ROAR205`. These were v0 ranking misses where the best true
  positive was already close to the cutoff; v1 lifts the correct phage to rank 1 in all four cases.
- **Still unpredictable:** `ECOR-06`, `ECOR-69`, `FN-B4`, `H1-002-0060-C-T`, `H1-003-0088-B-J`, `IAI67`, `NILS24`,
  `NILS53`.
- **Why they remain hard:**
  - `FN-B4` and `NILS24` are abstention failures: zero positives means any top-3 recommendation is wrong.
  - `ECOR-06` and `H1-002-0060-C-T` are narrow-recall misses: each has one positive, but it stays outside the top 3.
  - `ECOR-69` and `NILS53` are family-collapse cases: v1 still pushes the wrong phage-family block to the top.
  - `H1-003-0088-B-J` and `IAI67` are new v1 misses where the best true positive only reaches ranks 4 and 6.

#### Interpretation

- The v1 feature blocks do fix a real subset of the v0 failures, especially the same-family ordering misses.
- The remaining misses are not a single bug class; they split cleanly into abstention, narrow-recall, and
  cross-family-collapse modes.
- The final model is better overall, but the honest deployment message is still strain-specific: some strains are now
  solvable, and some remain effectively unpredictable with this feature set.
