# Raw inputs

Immutable. Nothing in this directory is written by the pipeline; every derived
dataset goes to `data/processed/`.

## Files

| File | Contents | SHA-256 |
|---|---|---|
| `all_descriptor_smile.xlsx` | CDK physicochemical and topological descriptors with the SMILES of each compound | `4caded9c92e246b3…` |
| `datos.xlsx` | the merged compound table the original R workflow consumed, carrying IC₅₀ and the chemical-family attribute | `9ff8d91d43a646d6…` |
| `00koopmans_new.xlsx` | conceptual-DFT global reactivity descriptors from the B3LYP-D3/6-31G frontier orbital energies | `baf7ece81cc98c1a…` |
| `r_benchmark/*.xlsx` | the six result workbooks written by the original R benchmark, reconciled in step 04 | one digest per workbook, keyed `r_benchmark/<name>.xlsx` |

All nine digests are in `src/ureasecdft/config.py` (`INPUT_SHA256`, keyed by
path relative to this directory), and `scripts/01_build_datasets.py` verifies
every one of them before reading anything. A mismatch stops the run rather
than silently analyzing a different snapshot.

## Provenance

Urease inhibition data were compiled from the primary literature by the
authors; the compilation underlying `datos.xlsx` is the one described in the
2024 study of physicochemical filters for urease inhibitors
(*Int. J. Mol. Sci.* **25**, 4303, https://doi.org/10.3390/ijms25084303).
IC₅₀ values are in µM. Records with no reported IC₅₀ denote compounds for
which no inhibition was detected; they are treated as right-censored at 1 mM
(`CENSORED_IC50_UM` in `config.py`), which is a modeling decision of this
work, not a measured value.

The reactivity descriptors were computed by the authors: structures were drawn
in 2D, protonated at pH 7.2 with Epik, geometry-optimized at B3LYP-D3/6-31G in
Gaussian 16, and the global reactivity indices derived from the resulting
frontier orbital eigenvalues through the Koopmans approximation. Energies are
in hartree. The Gaussian input and output files are not redistributed here;
`00koopmans_new.xlsx` is the tabulated result, and step 02 re-derives every
descriptor from the three eigenvalues so the table can be checked without
them.

## Two cautions for reusers

**The reactivity descriptors are not uniformly physical.** 68 of the 600
records have ε_HOMO > 0 and 271 have ε_LUMO > 0, so the Koopmans ionization
potential is negative for 11% of the set and the electron affinity for 45%.
This follows from applying a gas-phase basis set without diffuse functions to
anionic microspecies. Step 02 writes the per-record diagnostic to
`results/tables/SI_T5_frontier_orbital_diagnostics.csv`. Do not reuse these
descriptors as physical ionization potentials or electron affinities.

**The 600 records are not 600 distinct compounds.** They collapse to 387
unique standardized InChIKeys and 86 Bemis–Murcko scaffolds. Any split that
ignores this leaks close analogues between training and test partitions. Step
01 writes the full duplication audit to
`data/processed/duplication_audit.csv`.

## descriptor_definitions.csv

Definitions for all 306 descriptors, joined onto Table S6 by `scripts/08_tables.py`.
Not an input to the analysis: nothing in steps 01-07 reads it, and it carries no data,
only documentation of what each descriptor column means.

Provenance of the 287 physicochemical entries. Harvested from the Chemistry Development
Kit, which is the implementation `rcdk` v3.6.0 calls:

* `descriptor-algorithms.owl`, the CDK descriptor dictionary, for class-level definitions.
* The javadoc of the CDK descriptor classes for the per-component meaning of each indexed
  name (for example `WTPT-4` is the sum of path lengths starting from oxygens, and
  `BCUTc-1h` is the highest Burden eigenvalue weighted by partial charge).
* `KierHallSmartsDescriptor` for the 79 `khs.*` entries, each reported with the SMARTS
  pattern whose matches it counts.

Two points a reader should know, both read off the CDK source rather than inferred:

* The autocorrelation descriptors are indexed from topological distance zero. `ATSc1`,
  `ATSm1` and `ATSp1` are the distance-0 terms; `ATSc5`, `ATSm5` and `ATSp5` are the
  distance-4 terms.
* The CDK javadoc for `CarbonTypesDescriptor` contains the typo "triply hound" for
  "triply bound", and the CPSA table prints `FFSA-2` for `FPSA-2`. Both are corrected
  here; no meaning was changed.

WHIM symbols are defined only to the level CDK documents them - 11 directional components
(lambda 1-3, nu 1-2, gamma 1-3, eta 1-3) and 6 non-directional (T, A, V, K, G, D), with
the symbol conventions attributed to Todeschini et al. All entries in this release use
unit atom weighting (`.unity`).

Provenance of the 19 CDFT entries: the equations audited in Table S7 and reproduced in
`scripts/02_cdft_audit_eda.py`, including the implied reference
epsilon_HOMO(tetracyanoethylene) = -0.344290 hartree used by the nucleophilicity index.
