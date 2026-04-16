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

**What we learned (SX11):** aggregate null (-1 pp mAP, AUC -3.4 pp because per-query scores aren't
cross-calibrated), but **within-family nDCG +3.48 pp with CI disjoint from baseline** — the
strongest stratum-specific signal in wave-2. The loss is sound; aggregate dilution by 69%
cross-family pairs hides it. Integrating with per-phage blending is a wave-3 candidate. See
`spandex-wave-2-baseline`, `stratified-eval-framework`.

## MLC (minimum lytic concentration)

Our potency label: integer 0–3 where 0 = no lysis observed, 3 = lysis observed at the lowest phage
concentration tested (5×10⁶ pfu/ml = most potent). Higher MLC = more potent phage–host interaction.
Derived from binary spot-test scores at three replicated dilutions (5×10⁸, 5×10⁷, 5×10⁶); the
unreplicated 5×10⁴ dilution is excluded. Post-SX05, MLC=4 from the paper's matrix is capped to
MLC=3 because our binary data cannot reconstruct the paper's plaque-morphology split at 5×10⁶. See
`mlc-dilution-potency`.

MLC is **ordinal**, not metric: 3 > 2 > 1 > 0 is meaningful, but "distance from 1 to 2" is not
guaranteed to equal "distance from 2 to 3."

## nDCG (Normalized Discounted Cumulative Gain)

A ranking metric for graded relevance. For one bacterium's ranked list of phages:

- **Gain**: the true MLC of each phage (0, 1, 2, or 3).
- **Discounted**: gain at rank k is divided by log₂(k+1) — hits at the top count more than hits
  deep in the list.
- **Cumulative**: sum the discounted gains across the list.
- **Normalized**: divide by the DCG of the ideal ranking (all true positives sorted by MLC
  descending). Result ∈ [0, 1], where 1 = perfect ranking.

Why we use it: naturally handles **graded** relevance (MLC=3 is "more relevant" than MLC=1), unlike
top-3 hit rate which treats all positives identically. See `top3-metric-retired`.

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

**Our canonical setup**: 10-fold CV over the 148 unified panel phages (SX15), stratified by ICTV
family so each family appears in both train and test. Result: aggregate AUC 0.8988 (higher than
bacteria-axis!) but nDCG 0.7229 (7 pp lower). Reassuring sub-finding: BASEL phages generalize as
well as Guelin phages (AUC 0.896 vs 0.899) — feature engineering transfers across phage collections.

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

**Why we care**: SX11 LambdaRank's +3.48 pp within-family nDCG was invisible in aggregate because
cross-family pairs diluted it to zero. Every ticket from SX14 onward reports stratified metrics
alongside aggregate. See `stratified-eval-framework`, `spandex-wave-2-baseline`.

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
