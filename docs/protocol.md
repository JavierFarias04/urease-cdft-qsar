# Analysis protocol

Fixed before the analysis was run, and reproduced here unchanged. Numbers in
the text come from `results/tables/`; nothing in this document is retyped from
the manuscript.

## 1. Compound set

Urease inhibition records with a reported IC₅₀ in µM, plus records for which no
inhibition was detected. The latter are treated as right-censored at 1 mM. This
is a modeling decision, not a measurement: the underlying reports state only
that no inhibition was observed at the concentrations tested.

Structures are standardized with RDKit — fragment parent, cleanup, and charge
neutralisation — with stereochemistry retained, then grouped by full InChIKey.
Within a group of identical structures carrying different IC₅₀ values the rule
is: a measured value takes precedence over a censored one; two measured values
within one log unit are replaced by their geometric mean; a wider disagreement
drops the group. No group was dropped on this dataset.

## 2. Descriptor blocks

Three blocks, all computed on the **same** standardized structures, so
descriptor content is the only factor that differs between them.

- **Physicochemical and topological** — CDK descriptors from the 2D structure.
- **CDFT global reactivity** — global reactivity indices derived from the
  B3LYP-D3/6-31G frontier orbital eigenvalues through the Koopmans
  approximation, in hartree.
- **Combined** — the union of the two.

Each block is run twice: with the chemical-family attribute supplied as an
explicit categorical feature, and with it withheld.

## 3. Activity definitions

Six, all applied to the same structures:

- single cutoff at 5, 10 and 25 µM — inhibitor if IC₅₀ ≤ cutoff;
- two cutoff at 5, 10 and 25 µM against 50 µM — inhibitor if IC₅₀ ≤ the lower
  cutoff, non-inhibitor if IC₅₀ > 50 µM, compounds between the two thresholds
  excluded.

Under the two-cutoff definitions the non-inhibitor class is defined
invariantly, so its size does not change with the lower threshold.

## 4. Partitioning

Three protocols, applied to the same model population:

1. **600 records, random 80/20** — reproduces the protocol of the earlier
   literature, duplicates included.
2. **387 structures, random 80/20** — isolates the effect of duplicate records.
3. **387 structures, Bemis–Murcko scaffold split** — no scaffold appears in
   both the training and the test partition.

The third is the one to read. The set contains 86 scaffolds for 600 records and
83% of structures have an ECFP4 Tanimoto neighbor above 0.7, so a random split
places close analogues on both sides of the partition.

## 5. Preprocessing

Fitted on the training partition only, in this order: drop descriptors missing
in at least 80% of training rows; median-impute the remainder; remove
near-zero-variance descriptors by the caret criterion (frequency ratio 19,
unique-value fraction 10%); standardize to zero mean and unit variance; apply
SMOTE to the training rows only.

Nothing is fitted on the test partition, and SMOTE never touches it.

## 6. Feature selection

Three regimes: none; Boruta (20 iterations, 200-tree random forest, one-sided
binomial test with Bonferroni correction at α = 0.05); and retention of the
descriptors carrying the first 95% of cumulative XGBoost gain.

## 7. Classifiers

Seven: random forest, XGBoost, decision tree (CART), support vector machine
with an RBF kernel, k-nearest neighbors, naïve Bayes, logistic regression.

**Hyperparameters are fixed and identical across the three descriptor blocks.**
This is deliberate: unequal tuning effort between blocks would confound a
comparison whose only intended factor is descriptor content. The consequence is
that absolute performance is lower than a tuned model would reach, and the
figures should be read as block contrasts, not as a performance ceiling.

CART substitutes for the C5.0 implementation used in the earlier R workflow.

## 8. Replication and statistics

Ten repeats per configuration, each with its own seed derived from the master
seed. The reported metric is the Matthews correlation coefficient; balanced
accuracy, sensitivity, specificity and ROC AUC are also recorded.

Blocks are compared **paired on matched configurations** — the same tuple of
(family condition, activity definition, selection regime, classifier, repeat) —
with a Friedman test across the three blocks, Nemenyi post-hoc comparisons, and
Wilcoxon signed-rank tests for each pairwise contrast. Comparing maxima over
model populations instead of paired configurations inflates differences between
blocks and is not used here.

## 9. Interpretation

SHAP attribution (TreeExplainer) on gradient-boosted models fitted to the
combined block under scaffold splitting, aggregated by descriptor origin.
Feature-selection stability is the fraction of runs in which each descriptor
was retained.

Attribution identifies what a fitted model used. It is not evidence of a
chemical mechanism, and no causal claim is drawn from it.
