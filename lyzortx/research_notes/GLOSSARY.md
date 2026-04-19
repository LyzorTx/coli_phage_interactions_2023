# Glossary

Plain-language definitions of terms that recur across lab notebooks, knowledge units, and PR discussions.
Entries are ordered alphabetically. Keep definitions short, concrete, and grounded in how we use the term in
*this* project — not textbook generality. Link to the source notebook / knowledge unit where the term
originated.

---

## AUC (Area Under the ROC Curve)

A threshold-independent measure of how well a binary classifier separates positives from negatives.
Range [0, 1], where 0.5 = random guessing, 1.0 = perfect separation, <0.5 = worse than random
(flip the predictions).

**Mechanism**: take every possible (positive, negative) pair in your evaluation set. AUC is the
fraction of pairs where the model assigns a higher score to the positive than the negative. An AUC
of 0.87 means: if you pick a random lysing pair and a random non-lysing pair, 87% of the time the
model scores the lysing one higher.

**Why we use it**: threshold-free (don't have to pick a P(lysis) cutoff), insensitive to class
imbalance (our data is 79% MLC=0). Dominant discrimination metric across the project.

**What AUC misses**: it's a **global** pairwise metric — it pools all pairs across all bacteria.
A model can have great AUC but still bury the one real positive for a specific bacterium below
broad-host generalists. That's why we also track per-bacterium nDCG / mAP, which reward ranking
quality **within a query**. SX15 showed the two can diverge: phage-axis AUC 0.8988 with nDCG
0.7229 (7 pp gap = the per-phage blending tax).

**Rough interpretation scale**: 0.5 random, 0.7 weak-but-usable, 0.8 good, 0.9 very strong. Our
within-panel baseline is 0.87, cross-panel 0.76.

## Bacteria-axis evaluation

A cross-validation split where **bacteria are held out** but all phages remain in training. For
each fold, a subset of bacteria goes to test; the model is trained on all pairs involving the
remaining bacteria, then asked to predict lysis on the held-out bacteria × all phages.

**What it tests**: generalization to unseen hosts. "If I sequence a new clinical *E. coli* isolate,
how well can I predict which of my catalogued phages will lyse it?"

**What stays easy**: every phage was in training — the model has a learned per-phage representation,
can use per-phage blending (AX02), and has seen how each phage behaves across other hosts.

**Our canonical setup**: 10-fold CV over the 369 unified panel bacteria (SX15), stratified by
(source, phylogroup). Result: aggregate nDCG 0.7965, AUC 0.8685 — essentially identical to SX10.
See `spandex-unified-kfold-baseline`. Contrast with **phage-axis evaluation**.

## Brier score

A probability **calibration** metric for binary predictions. Mean squared error between predicted
probability and the binary outcome: `Brier = mean((P_pred - y_true)²)`.

Range [0, 1]. **Lower is better.** Perfect: 0 (every prediction is 0 or 1 and matches the outcome).
A constant predictor at the base rate (0.21 for us) gets ~0.16. Always-0.5 gets 0.25.

**Why it's different from AUC**: AUC only cares about the *ranking* of predictions. A model that
outputs 0.501 for every positive and 0.499 for every negative has AUC = 1.0 (perfect ranking) but
terrible calibration — the probabilities don't reflect true event rates. Brier punishes that; it
rewards predictions where P_pred = 0.7 actually corresponds to a 70% lysis rate in reality.

**Why we track both**: AUC tells us if we rank right, Brier tells us if our probability outputs
are trustworthy for downstream decisions (cocktail composition, threshold tuning, expected-value
calculations). A well-calibrated 0.3 prediction means "30% chance of lysis" — actionable. An
uncalibrated 0.3 from LambdaRank means only "ranked between a 0.2 and a 0.4 within this query" —
not actionable across bacteria.

**What happens when it moves**: SX11 LambdaRank's Brier rose 0.1279 → 0.1508 (+2.3 pp) because
per-query ranking normalization destroys global probability meaning. SX10 baseline sits at 0.1248;
CatBoost (GT05) improved Brier to 0.152 vs LightGBM's 0.161 on GT03 features — useful for
calibration even when AUC was flat.

**Rough interpretation scale** (for our class imbalance, base rate ~0.21): 0.16 = predicting the
base rate uninformatively, 0.12–0.13 = our working range (informative but with room), 0.08 = would
be excellent, <0.05 = typical of much easier problems. Our SX10 Brier of 0.125 means we're
materially better-calibrated than guessing the base rate.

## Bootstrap confidence interval

A non-parametric way to put a confidence interval on a statistic (like AUC or Brier) when there's
no clean closed-form formula for its sampling distribution. The trick: treat your observed sample
as if it were the full population, then simulate "another draw" by **sampling with replacement**
from your own data. Compute the statistic on each resample. The 2.5th and 97.5th percentiles of
the resulting 1,000 numbers are your 95% CI.

**Why the bootstrap at all, instead of a classical formula?** AUC is a U-statistic — it's the
rank of positives among negatives (see AUC entry), not a mean. Every (positive, negative) pair
contributes interactively, so classical variance formulas (Hanley-McNeil, DeLong) require
distributional assumptions that break under class imbalance, correlated predictions, or small
samples. Bootstrap sidesteps all of it: you never have to write down a formula for the sampling
distribution — you simulate it.

**Why "with replacement"**: without replacement, every resample is just a permutation of the
original — the statistic would never change. Replacement lets some observations appear multiple
times and others not at all, which is what genuine sampling variability looks like. On any given
resample, roughly 37% (1/*e*) of the original items are absent.

**Why bacterium-level, not pair-level, for us**: pairs within a single bacterium share the host
genome and host features. If the model's host encoding is wrong for bacterium X, all ~96 of X's
pair predictions are miscalibrated together — they're correlated, not independent. Bootstrap
theory requires independent resampling units. A pair-level bootstrap would hide this clustering
and produce **falsely narrow CIs**, because it'd almost never leave all of X's pairs out at once.
Effective sample size drops from pair-count (≈ 35,000) to bacterium-count (≈ 369) when you
account for the correlation — orders of magnitude. We resample 369 bacteria with replacement,
pull in all their pair predictions, compute AUC+Brier, repeat 1,000 times.

**Which axis to resample depends on what's held out**: bacteria-axis CV (SX10/CH04) holds
bacteria out, so bacteria are the independent unit — resample bacteria. Phage-axis CV (CH05)
holds phages out — resample phages. Both-axis (CH06) held-out region — potentially resample
both. The bootstrap unit tracks the generalization question the metric answers, not just the
observation count. See `bootstrap_auc_brier_by_bacterium` (CH04),
`bootstrap_spandex_cis` (SPANDEX era).

See also: `error-buckets` and `st03-canonical-benchmark` for how small holdouts (65 bacteria)
make CI width load-bearing for interpretation.

## Concentration-aware training unit

The atomic training row from Track CHISEL onward: one observation = one `(bacterium, phage,
log_dilution, replicate) → {0, 1}` tuple from `raw_interactions.csv`. Each row is a physically
observed plaque assay outcome; concentration enters the model as a plain numeric feature
(`log_dilution ∈ {0, -1, -2, -4}`), letting LightGBM learn its own interaction with host and phage
features rather than baking in a dose-response shape.

**Contrast with the SPANDEX frame**: SPANDEX rolled up all concentrations and replicates for a
(bacterium, phage) pair into a single `any_lysis` label and then sometimes used a derived MLC
grade (0–3) on top. CHISEL drops the rollup. The consequences:

- rows with `score='n'` are excluded from training instead of silently counted as negatives
- BASEL (single-concentration) observations fit the training frame natively — one row per pair
- each Guelin pair now contributes up to 9 rows (3 replicates at 5×10⁸ pfu/ml, 2 at 5×10⁷,
  3 at 5×10⁶, 1 at 5×10⁴), giving the model more signal per pair than a single rollup label

**Evaluation** is scored at each pair's highest observed concentration so within-panel and
cross-source results share a comparable operating point. See `label-policy-binary` for the
authoritative definition and CH02–CH04 in the plan for the migration sequence.

## Hurdle two-stage loss

A two-model composition for zero-inflated ordinal data (our MLC is 79% zeros).

1. A binary classifier for **does it lyse at all?** (y ≥ 1)
2. A regressor, trained **only on positives**, for **given it lyses, how potent?** (E[MLC | y≥1])

Combined score: `P(lyses) × E[MLC | lyses]`. Canonical statistical fix for zero-inflated ordinals.
Biologically decouples adsorption (the lysis/no-lysis question) from replication efficiency (the
potency question), which can have different determinants.

Tested in SX11 — +0.94 pp nDCG, sub-threshold. See `ordinal-regression-not-better`.

## LambdaRank

A LightGBM training objective (`objective="lambdarank"`) that directly optimizes a ranking metric
(nDCG) during training instead of predicting a label.

Mechanism: training pairs are grouped by query (bacterium). For each pair of phages (i, j) within
that bacterium where i has higher true MLC than j, the loss penalizes the model whenever its
predicted score for j exceeds i's. The penalty is **weighted by how much nDCG would improve** if
that swap were fixed — mistakes at the top of the ranking cost more than mistakes at the bottom.
That weight is the "lambda."

**Why we tried it:** our eval metric is per-bacterium nDCG; train for the metric you want to win.

**What we learned (SX11):** aggregate null on AUC (-3.4 pp regression, because per-query scores
aren't cross-calibrated). The stratum-specific within-family nDCG gain (+3.48 pp) does not survive
the CHISEL metric change — ranking metrics are retired, so this is no longer a live arm. See
`ordinal-regression-not-better`, `ranking-metrics-retired`.

## MLC (minimum lytic concentration)

**DEPRECATED as training/evaluation label from Track CHISEL onward.** CHISEL switches to the
(bacterium, phage, concentration, replicate) → {0, 1} atomic unit; MLC is retained only as
shorthand when summarising raw observations by eye. See `label-policy-binary` and
`Concentration-aware training unit` below.

Historical definition (SPANDEX and earlier): integer 0–3 potency label where 0 = no lysis observed,
3 = lysis observed at the lowest phage concentration tested (5×10⁶ pfu/ml = most potent). Higher
MLC = more potent phage–host interaction. Derived from binary spot-test scores at three replicated
dilutions (5×10⁸, 5×10⁷, 5×10⁶); the unreplicated 5×10⁴ dilution is excluded. Post-SX05, MLC=4 from
the paper's matrix is capped to MLC=3 because our binary data cannot reconstruct the paper's
plaque-morphology split at 5×10⁶. See `mlc-dilution-potency`.

MLC is **ordinal**, not metric: 3 > 2 > 1 > 0 is meaningful, but "distance from 1 to 2" is not
guaranteed to equal "distance from 2 to 3."

## nDCG (Normalized Discounted Cumulative Gain)

**DEPRECATED from Track CHISEL onward.** Ranking metrics (nDCG, mAP, top-3) are retired from the
scorecard — the CHISEL biological model predicts P(lysis | bacterium, phage, concentration); any
ranking or cocktail-selection policy lives in a downstream product layer. Scored on AUC + Brier
only. See `ranking-metrics-retired`.

Historical definition: a ranking metric for graded relevance. For one bacterium's ranked list of
phages:

- **Gain**: the true MLC of each phage (0, 1, 2, or 3).
- **Discounted**: gain at rank k is divided by log₂(k+1) — hits at the top count more than hits
  deep in the list.
- **Cumulative**: sum the discounted gains across the list.
- **Normalized**: divide by the DCG of the ideal ranking (all true positives sorted by MLC
  descending). Result ∈ [0, 1], where 1 = perfect ranking.

## Ordinal all-threshold loss

A loss that respects MLC's ordering without treating it as a real number. Instead of one model,
train **three binary classifiers** in parallel, each answering a threshold question:

- Model A: is y ≥ 1? (any lysis)
- Model B: is y ≥ 2?
- Model C: is y ≥ 3?

The three classifiers are implicitly coupled by the ladder (if y≥2, then y≥1). At prediction time
combine the probabilities. In SX11 we used equal-weight averaging:
`score = (P(y≥1) + P(y≥2) + P(y≥3)) / 3`.

**How it's "ordinal":** respects that 3 > 2 > 1 > 0 without assuming the gaps are metric — sidesteps
the mistake plain regression makes on ordered categorical labels.

**What we learned (SX11):** best alternative arm at +1.33 pp nDCG over SX11's binary baseline but
below the +2 pp gate. Within-positive potency ranking improves (Kendall tau 0.208→0.290, 62% of
bacteria improve). But equal-weight averaging **demotes 55% of MLC=1 pairs** because a borderline
broad-host MLC=0 phage can earn higher P(y≥2)/P(y≥3) from training on other bacteria. Net: cleaner
2-vs-3 separation, weaker lysis/no-lysis boundary. See `ordinal-regression-not-better`.

## Per-phage blending (and the "blending tax")

An architectural trick from AX02 that adds a second predictor on top of the all-pairs model.

**How it works**: the production model has two components.

1. **All-pairs LightGBM**: one model trained on all (bacterium, phage) pairs using features that
   describe both sides. Works for any phage with features computed, including unseen ones.
2. **Per-phage LightGBM**: for each phage in training, a small dedicated 32-tree LightGBM fit on
   **host-only features** — answering "does host X support lysis by *this specific phage*?" At
   inference, the per-phage score is blended 50/50 with the all-pairs score.

Per-phage blending contributes the dominant architectural gain over the all-pairs baseline (+2 pp
AUC on ST03 holdout). See `per-phage-blending-dominant`.

**Why it's not universally deployable**: per-phage sub-models require training pairs for that
phage. For a newly-isolated phage with zero training rows, there's nothing to fit. At deployment
on unseen phages we can only use the all-pairs component.

**The "blending tax"**: the performance gap between bacterium-axis CV (every phage was in training
— can use blending) and phage-axis CV (held-out phages have no sub-model — can't). Under SX10's
canonical nDCG scoring the gap was ~7 pp; under AUC it's ~3 pp (0.870 bacteria-axis vs 0.899
phage-axis). Either way it's the cost paid per unseen phage and the dominant deployability bound
for cocktail-against-new-phage use cases. See `per-phage-not-deployable`, `deployment-goal`,
`spandex-unified-kfold-baseline`.

Note that "tax" is our shorthand for "structural cost from a missing component," not a modeling
failure — there is no fix without per-phage training data.

## Phage-axis evaluation

The symmetric cross-validation split to bacteria-axis: **phages are held out** but all bacteria
remain in training. For each fold, a subset of phages goes to test; the model is trained on all
pairs involving the remaining phages, then asked to predict lysis on all bacteria × held-out phages.

**What it tests**: generalization to unseen phages — the actual deployment scenario. "I'm
considering adding a newly-isolated phage to my catalogue. Before running any wet-lab screens,
which hosts will it likely lyse?"

**What becomes hard**: per-phage blending (AX02) is **impossible** for held-out phages — there are
zero training pairs for them, so no per-phage sub-model can be trained. The model must fall back
to its all-pairs component, relying entirely on phage-side features (TL17 family projection, depo
cross-terms, phage_stats) to represent the unseen phage. That's why phage-axis is the honest
deployability floor. See `per-phage-not-deployable`, `deployment-goal`.

**Our canonical setup**: 10-fold CV over the 148 unified panel phages (CH05), stratified by ICTV
family + "other" (<10 phages) + "UNKNOWN" (no family) catch-all. Result under CH05: aggregate AUC
0.8850, Brier 0.1348. Cross-source AUC parity holds (Guelin 0.8863 vs BASEL 0.8818) **but**
cross-source calibration diverges (Brier 0.1329 vs 0.1884) — the earlier SX15 claim "BASEL phages
generalize as well as Guelin" was a discrimination-only finding, not deployment readiness. See
`chisel-unified-kfold-baseline` for the three separate findings (discrimination parity, calibration
divergence, BASEL bacteria-axis deficit).

**Why AUC goes up but nDCG goes down vs bacteria-axis**: AUC pools all pairs globally and the
all-pairs features preserve lysis/no-lysis discrimination cleanly. nDCG is per-bacterium; each
bacterium now has ~15 held-out phages with no per-phage model, so within-list ordering is noisier.
~7 pp nDCG is the **per-phage blending tax** — the dominant deployability constraint.

**The third axis** (both bacterium and phage held out) is the hardest scenario and isn't in our
current evaluation. Not blocked by panel size — a 10×10 double-CV on 148 × 369 still yields ~330
observed pairs per cell, enough for pair-level bootstrap. Not prioritized because (a) 10× compute
cost vs single-axis, (b) the deployment scenarios of interest are single-axis (new phage against
catalogued hosts, or new host against catalogued phages), and (c) per-bacterium ranking gets
tight (~3 positives per held-out bacterium in each cell).

## Phage projection (TL17 BLAST features, "zero-vector projection")

`phage_projection` is a slot of ~33 phage-side features derived by BLASTing each phage's protein
complement against TL17, a Guelin-derived reference bank of protein families. Each feature reports
whether the phage has a protein hit in a specific TL17 family (thresholded bit-score or
presence-absence). It is the dominant phage-side feature family — roughly the phage analogue of
the host_surface slot.

**Zero-vector projection** (or "zero projection") describes a phage whose full `phage_projection`
feature vector is identically zero: none of its proteins matched any TL17 family above threshold.
The BLAST ran and produced no hits — this is a finding about the phage's genomic novelty relative
to the Guelin reference bank, not a missing-data bug.

**Why it's load-bearing**: under the phage-axis split, zero-vector phages calibrate *correctly*
because the model has no phage signal to misuse and falls back to the host-side prior. Non-zero
phages whose projection vectors map into Guelin-calibrated TL17 neighborhoods can be actively
*worse* calibrated than zero-vector phages when their actual host ranges don't match their
Guelin-neighbor priors. CH05 documented this for BASEL: 13 zero-vector BASEL phages (Brier 0.12)
vs 39 non-zero BASEL phages (Brier 0.31) on bacteria-axis. See the CH05 entry under
`chisel-unified-kfold-baseline` and `plm-rbp-redundant` (the same mechanism at cross-family scale).

**Practical consequence for future work**: more phages in TL17-underpopulated neighborhoods would
help more than engineering new phage features — the deficiency is reference-bank coverage, not
feature representation.

## Spearman correlation (Spearman's ρ)

A measure of how monotonically two variables move together, ranging from −1 to +1.

- +1 = perfectly monotonic increasing (X up ⇒ Y always up)
- 0 = no monotonic relationship
- −1 = perfectly monotonic decreasing

**The trick**: instead of correlating raw values, correlate their **ranks**. Sort X smallest to
largest and label 1, 2, 3…; do the same for Y; compute Pearson correlation on the rank labels.
Robust to non-linear-but-monotonic relationships and to outliers.

**Why we use it here**: MLC is ordinal, so Pearson-correlating MLC with P(lysis) is dishonest
(assumes equal distances between grades). Spearman asks the right question: "when MLC goes up,
does P(lysis) tend to go up too?"

**Rough interpretation scale**: 0.1 trivial, 0.25 weak-but-real, 0.5 strong, 0.7+ very strong. Our
SX11 finding of Spearman 0.246 between binary P(lysis) and MLC among positives means the binary
model is a **rough** potency ranker — getting it directionally right with real noise.

## Strata (SX14 stratified evaluation framework)

Categorical buckets for decomposing aggregate metrics. Introduced SX14 after realizing aggregate
nDCG can hide real stratum-specific signals. Four strata, each a per-pair boolean:

- **within_family**: phage's family has ≥3 training-positive pairs on this bacterium's cv_group.
  The model has family-specific priors to refine. (~28% of pairs)
- **cross_family**: zero training-positive pairs for that family × cv_group. Cold-start on the
  phage side. (~69% of pairs)
- **narrow_host_phage**: phage has <30% panel-wide lysis rate. Specialist phages. (~72% of pairs —
  overlaps with both above)
- **phylogroup_orphan**: holdout bacterium has ≤2 phylogroup siblings in its CV fold. Rare; often
  too small for reliable bootstrap (~1.6%).

**Why we kept them**: under the CHISEL AUC+Brier scorecard the strata are rescoped to a narrow-host
diagnostic rather than mandatory reporting — aggregate and stratified metrics are much more aligned
once ranking metrics are retired. The four-stratum decomposition still helps diagnose where a
feature is moving predictions and is especially useful for investigating narrow-host-phage rescue.
See `stratified-eval-framework`.

Not to be confused with **lysis-rate deciles** used by the `case-by-case` skill (10 buckets of
phage panel-wide lysis rate, 0-10%, 10-20%, …, 90-100%) — those are a separate tool for
per-bacterium audits.

## Within-positive potency ranking

Among the phages that DO lyse a given bacterium, do we rank the more potent (higher MLC) ones
higher than the less potent (lower MLC) ones? Measured by Kendall tau or Spearman ρ on (predicted
score vs true MLC), restricted to pairs with MLC > 0.

Distinct from the **lysis/no-lysis boundary** (separating MLC=0 from MLC≥1), which is what binary
AUC and the top-k ranking mostly depend on. SX11 showed the two goals trade off under ordinal
losses: +3.5 pp within-family potency ranking, -1 pp mAP at the lysis boundary.
