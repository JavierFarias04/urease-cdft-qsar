# Quantum-chemical reactivity descriptors add no detectable predictive value for *Helicobacter pylori* urease inhibition

Data and code to reproduce every number, table and figure in the accompanying
article.

## What the study found

Conceptual-DFT global reactivity descriptors, computed through the Koopmans
approximation from B3LYP-D3/6-31G frontier orbital energies, were added to a
block of 287 physicochemical and topological descriptors and the two blocks
were compared on **identical compounds** under matched model configurations.

Adding the reactivity block changes the Matthews correlation coefficient by a
median of **0.000** under all three evaluation protocols. Under the strictest
protocol — a Bemis–Murcko scaffold split — the mean change is **−0.012 ±
0.005** over 2,520 matched configurations, and the combined block is the better
of the two in only 34.5% of them. The reactivity block on its own is
informative but clearly weaker than the physicochemical block.

Three further results bound how the original dataset should be used:

- The 600 records describe **387 unique structures** and **86 Bemis–Murcko
  scaffolds**; 83% of structures have an ECFP4 Tanimoto neighbor above 0.7. A
  random split therefore places close analogues on both sides of the partition,
  and scaffold-based splitting costs 27–39% of the MCC depending on the block.
- The apparent benefit of supplying the chemical family as a feature
  (+0.25 MCC on the 600 records) falls to **0.000** once duplicate records are
  removed. It was a duplication artifact.
- The reactivity block is nominally 19-dimensional but effectively
  3-dimensional: three principal components carry 95% of its variance, and
  feature selection retains only the frontier-orbital gap and its algebraic
  relatives. SHAP attribution assigns it 5.1 ± 0.6% of the total.

The conclusion is about **redundancy, not irrelevance**. The gap is the
electronic quantity that carries predictive signal; the remaining indices are
non-linear reexpressions of the same two eigenvalues and add nothing to
descriptors that already encode molecular size, shape and polarity.

## Layout

```
data/raw/          three input workbooks + six original R benchmark results, all checksummed
data/processed/    standardized compounds, duplication audit, similarity matrix
scripts/           01…10 in run order, fig1…fig6, run_all.sh
src/ureasecdft/    config (paths, seeds, thresholds), input checks, figure style
results/tables/    every result table as CSV
results/figures/   the six publication figures, PNG at 400 dpi and vector PDF
                   (no fitted models are serialised: hyperparameters are fixed,
                    so any model is cheaper to refit than to store)
docs/              analysis protocol, data dictionary, main-text tables as markdown
```

`scripts/` holds the pipeline order; `src/` holds the shared logic and every
frozen parameter. No script contains a hardcoded path or a retyped result.

## Setup

```bash
conda env create -f environment.yml && conda activate ureasecdft
# or: python -m venv .venv && source .venv/bin/activate && pip install -r requirements.txt
```

RDKit is the one dependency that is markedly easier to install through conda.

## Running

```bash
bash scripts/run_all.sh            # everything, in order
N_JOBS=4 bash scripts/run_all.sh   # fewer cores
```

Step 01 verifies the SHA-256 of all nine raw inputs — the three input workbooks
and the six original R benchmark workbooks — against `src/ureasecdft/config.py`
and stops if any differs, so a run cannot silently proceed on a different data
snapshot.

Expect about **75 minutes** in total on 10 cores and under 4 GB of memory.
Step 03 dominates: it fits 22,680 models (3 protocols × 3 blocks × 2 family
conditions × 6 activity definitions × 3 selection regimes × 7 classifiers × 10
repeats). Every other step runs in seconds to a few minutes.

## What each step does

1. `01_build_datasets.py` — verify input checksums; standardize structures with
   RDKit; group by InChIKey and resolve IC₅₀ conflicts; audit duplication at
   four levels of structural identity; assemble the three descriptor blocks.
2. `02_cdft_audit_eda.py` — re-derive all fifteen reactivity descriptors from
   the three orbital eigenvalues and check them against the tabulated values;
   PCA effective dimensionality; correlation structure; ECFP4 fingerprints and
   the nearest-neighbor similarity distribution; scaffold populations;
   frontier-orbital physicality diagnostics.
3. `03_benchmark.py` — the controlled factorial benchmark. Preprocessing is
   fitted on the training partition only and SMOTE never touches the test rows.
4. `04_reconcile_R_results.py` — tidy the six original R workbooks into one
   table of 756 models and apply the same paired statistics, so the central
   result does not depend on this reimplementation.
5. `05_statistics.py` — Friedman, Nemenyi and Wilcoxon on matched
   configurations.
6. `06_interpretation.py` — SHAP attribution by descriptor origin;
   feature-selection stability.
7. `07_roc.py` — pooled ROC curves with percentile bands.
8. `08_tables.py` — build every result table, with arithmetic self-checks.
9. `09_md_tables.py` — render the main-text tables as markdown.
10. `10_data_dictionary.py` — build `docs/data_dictionary.md` from the files on
    disk, so it cannot drift from the data it describes.

`fig1`…`fig6` build the publication figures from `results/tables/` and
`data/processed/`.

## Reproduction

Every script in this repository was executed end to end on the inputs in
`data/raw/` before release. Values reproduced exactly:

| Quantity | Published | This repository |
|---|---|---|
| records analyzed | 600 | 600 |
| unique canonical SMILES | 475 | 475 |
| unique InChIKeys | 387 | 387 |
| unique connectivity skeletons | 331 | 331 |
| Bemis–Murcko scaffolds (600 records) | 86 | 86 |
| median nearest-neighbor Tanimoto | 0.820 | 0.820 |
| structures with a neighbor > 0.7 | 82.7% | 82.7% |
| records with ε_HOMO > 0 | 68 | 68 |
| records with ε_LUMO > 0 | 271 | 271 |
| physicochemical block (computed / constant / informative) | 287 / 61 / 226 | 287 / 61 / 226 |
| combined block | 306 / 62 / 244 | 306 / 62 / 244 |
| reactivity block | 19 / 1 / 18 | 19 / 1 / 18 |

The benchmark was re-executed in full for this release (22,680 fits, 70.9 min on
10 cores) and the paired contrasts recomputed from it:

| Quantity | Published | This repository |
|---|---|---|
| median ΔMCC, combined − physicochemical (all three protocols) | 0.000 | 0.000 |
| mean ΔMCC, scaffold split | −0.012 ± 0.005 | -0.012 ± 0.005 |
| configurations favoring the combined block, scaffold split | 34.5% | 34.5% |
| Wilcoxon p, scaffold split | 2.1 × 10⁻⁷ | 2.1 × 10⁻⁷ |
| mean ΔMCC, 600 records random split | −0.009 | -0.009 |
| mean ΔMCC, 387 structures random split | +0.002 | +0.002 |
| matched configurations per protocol | 2,520 | 2,520 |
| SHAP share, reactivity block | 5.1 ± 0.6% | 5.1 ± 0.6% |
| SHAP share, physicochemical block | 94.6 ± 0.6% | 94.6 ± 0.6% |
| SHAP share, chemical-family indicators | 0.27 ± 0.09% | 0.27 ± 0.09% |
| scaffold-split AUC, reactivity block | 0.676 ± 0.087 | 0.676 ± 0.087 |
| scaffold-split AUC, physicochemical block | 0.925 ± 0.035 | 0.925 ± 0.035 |
| scaffold-split AUC, combined block | 0.928 ± 0.037 | 0.928 ± 0.037 |
| model fits, step 03 | 22,680 (0 failures) | 22,680 (0 failures) |

Every value above was read back from `results/tables/` after the run, not
copied from the manuscript. No published number failed to reproduce.


## Methodological caveats a reuser must know

**Hyperparameters are fixed, not tuned.** They are identical across the three
descriptor blocks on purpose: unequal tuning effort would confound a comparison
whose only intended factor is descriptor content. Absolute performance is
therefore lower than a tuned model would reach. Read the results as block
contrasts, not as a performance ceiling.

**The reactivity descriptors are partly unphysical.** 11% of records have a
positive ε_HOMO and 45% a positive ε_LUMO, so the Koopmans ionization potential
and electron affinity are negative for those records. This follows from a
gas-phase basis set without diffuse functions applied to anionic microspecies
at pH 7.2. Whether a basis set with diffuse functions and implicit solvation
would change the central result is **untested** and is the most important open
question about this work. Do not reuse these columns as physical IP or EA
values.

**Censoring is an assumption.** Records with no reported IC₅₀ are treated as
right-censored at 1 mM. The sources state only that no inhibition was observed
at the concentrations tested.

**Attribution is not mechanism.** SHAP values identify what a fitted model
used. Nothing here supports a claim about how a compound inhibits urease.

**The earlier 677-compound set could not be reconstructed.** The 2024 study
reports 677 compounds and 207 variables; the files available here yield 600
records and 387 unique structures, and this repository reports the latter. The
discrepancy is not resolved.

**Row-order coupling.** `data/processed/tanimoto_similarity_dedup.npy` is a
square matrix whose row order follows `data/processed/compounds_dedup.csv`.
Join through that file; do not assume its order matches any other table.

## Known differences on re-running

- **Figure dimensions differ from the published PNGs by design.** The published
  figures were saved with matplotlib's `savefig.bbox="tight"`, which crops each
  canvas to its own content, so the six figures came out between 6.67 and
  7.26 in wide despite a declared width of 7.20 in. This repository sets
  `savefig.bbox=None` (`src/ureasecdft/style.py`), so all six are exactly
  7.20 in wide. The plotted content is unchanged; only the surrounding margin
  is. No reported number depends on this.
- **Stochastic steps.** The train/test splits, SMOTE resampling, random-forest
  and XGBoost fitting, and Boruta's shadow features are all seeded from
  `config.SEED`, so a rerun on the same library versions reproduces the tables
  exactly. Across different scikit-learn, XGBoost or RDKit versions the
  per-configuration MCC values may shift in the third decimal. The central
  result is a median difference of 0.000 over thousands of matched
  configurations and a sign test over them, so it does not depend on any single
  fit.
- **C5.0 is not available in Python.** Step 03 uses CART where the original R
  workflow used C5.0. Step 04 reconciles against the original R results, so
  both implementations are represented.

## Citation

See `CITATION.cff`. Author names, affiliations, ORCIDs, the article DOI and the
copyright holder are marked `[PENDING]` and must be filled in before the
repository is made public.

## License

MIT, see `LICENSE`. The compound data were compiled from the primary
literature; please cite the original reports for individual measurements.
