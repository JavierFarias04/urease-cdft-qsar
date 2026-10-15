# Data dictionary

Built from the files on disk by `scripts/10_data_dictionary.py`. Wide descriptor tables list their first columns only; the remainder keep the names used by the descriptor software.

Row order in `data/processed/tanimoto_similarity_dedup.npy` follows `data/processed/compounds_dedup.csv`. Join through that file rather than assuming row order matches any other table.

## Raw inputs

### `data/raw/00koopmans_new.xlsx`

649 rows x 21 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `homo_energy` | float64 | hartree |  |
| `homo_.1_energy` | float64 | hartree |  |
| `lumo_energy` | float64 | hartree |  |
| `gap` | float64 | hartree |  |
| `potential_ionization` | float64 | hartree |  |
| `affinity_electronic` | float64 | hartree |  |
| `chemical_potential` | float64 | hartree |  |
| `hardness` | float64 | hartree |  |
| `softness` | float64 | hartree⁻¹ |  |
| `electronegativity` | float64 | hartree |  |
| `electrophilicity` | float64 | hartree |  |
| `nucleophilicity` | float64 | hartree |  |
| … | | | 9 further descriptor columns, named as in the source workbook |

### `data/raw/all_descriptor_smile.xlsx`

600 rows x 309 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `tipo_familia` | str |  | chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros) |
| `homo_energy` | float64 | hartree |  |
| `homo_.1_energy` | float64 | hartree |  |
| `lumo_energy` | float64 | hartree |  |
| `gap` | float64 | hartree |  |
| `potential_ionization` | float64 | hartree |  |
| `affinity_electronic` | float64 | hartree |  |
| `chemical_potential` | float64 | hartree |  |
| `hardness` | float64 | hartree |  |
| `softness` | float64 | hartree⁻¹ |  |
| `electronegativity` | float64 | hartree |  |
| `electrophilicity` | float64 | hartree |  |
| … | | | 297 further descriptor columns, named as in the source workbook |

### `data/raw/datos.xlsx`

628 rows x 290 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `molecula` | str |  |  |
| `Nombre real` | str |  |  |
| `IC50` | float64 |  |  |
| `Fsp3` | float64 |  |  |
| `nSmallRings` | int64 |  |  |
| `nAromRings` | int64 |  |  |
| `nRingBlocks` | int64 |  |  |
| `nAromBlocks` | int64 |  |  |
| `nRings3` | int64 |  |  |
| `nRings4` | int64 |  |  |
| `nRings5` | int64 |  |  |
| `nRings6` | int64 |  |  |
| … | | | 278 further descriptor columns, named as in the source workbook |

## Processed datasets

### `data/processed/benchmark_raw_results.csv`

22680 rows x 17 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `repeat` | int64 |  |  |
| `selection` | str |  |  |
| `n_selected` | int64 |  |  |
| `algorithm` | str |  |  |
| `n_train` | int64 |  |  |
| `n_test` | int64 |  |  |
| `accuracy` | float64 |  |  |
| `balanced_accuracy` | float64 |  |  |
| `mcc` | float64 | dimensionless |  |
| `f1` | float64 |  |  |
| … | | | 3 further descriptor columns, named as in the source workbook |

### `data/processed/block_composition.csv`

3 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `n_descriptors` | int64 |  |  |
| `n_constant` | int64 |  |  |
| `n_informative` | int64 |  |  |

### `data/processed/compounds_dedup.csv`

387 rows x 320 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `tipo_familia` | str |  | chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros) |
| `homo_energy` | float64 | hartree |  |
| `homo_.1_energy` | float64 | hartree |  |
| `lumo_energy` | float64 | hartree |  |
| `gap` | float64 | hartree |  |
| `potential_ionization` | float64 | hartree |  |
| `affinity_electronic` | float64 | hartree |  |
| `chemical_potential` | float64 | hartree |  |
| `hardness` | float64 | hartree |  |
| `softness` | float64 | hartree⁻¹ |  |
| `electronegativity` | float64 | hartree |  |
| `electrophilicity` | float64 | hartree |  |
| `nucleophilicity` | float64 | hartree |  |
| `acceptor_power` | float64 | hartree |  |
| … | | | 306 further descriptor columns, named as in the source workbook |

### `data/processed/compounds_raw.csv`

600 rows x 318 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `tipo_familia` | str |  | chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros) |
| `homo_energy` | float64 | hartree |  |
| `homo_.1_energy` | float64 | hartree |  |
| `lumo_energy` | float64 | hartree |  |
| `gap` | float64 | hartree |  |
| `potential_ionization` | float64 | hartree |  |
| `affinity_electronic` | float64 | hartree |  |
| `chemical_potential` | float64 | hartree |  |
| `hardness` | float64 | hartree |  |
| `softness` | float64 | hartree⁻¹ |  |
| `electronegativity` | float64 | hartree |  |
| `electrophilicity` | float64 | hartree |  |
| `nucleophilicity` | float64 | hartree |  |
| `acceptor_power` | float64 | hartree |  |
| … | | | 304 further descriptor columns, named as in the source workbook |

### `data/processed/duplication_audit.csv`

4 rows x 9 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `level` | str |  |  |
| `n_unique` | int64 |  |  |
| `n_duplicate_groups` | int64 |  |  |
| `n_rows_in_duplicate_groups` | int64 |  |  |
| `n_groups_discordant_IC50` | int64 |  |  |
| `n_groups_mixing_censored_and_measured` | int64 |  |  |
| `n_groups_discordant_class_5uM` | int64 |  |  |
| `n_groups_discordant_class_10uM` | int64 |  |  |
| `n_groups_discordant_class_25uM` | int64 |  |  |

### `data/processed/feature_selection_records.csv`

24857 rows x 7 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `method` | str |  |  |
| `scheme` | str |  |  |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `repeat` | int64 |  |  |
| `feature` | str |  |  |
| `n_candidates` | int64 |  |  |

### `data/processed/original_R_paired_mcc.csv`

252 rows x 7 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `cdft` | float64 |  |  |
| `physchem` | float64 |  |  |
| `combined` | float64 |  |  |

### `data/processed/original_R_results_tidy.csv`

756 rows x 17 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `algorithm_raw` | str |  |  |
| `accuracy` | float64 |  |  |
| `sensitivity` | float64 |  |  |
| `specificity` | float64 |  |  |
| `precision` | float64 |  |  |
| `f1` | float64 |  |  |
| `balanced_accuracy` | float64 |  |  |
| `mcc` | float64 | dimensionless |  |
| `strategy_raw` | str |  |  |
| `selection_raw` | str |  |  |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `source_file` | str |  |  |
| `n_compounds` | int64 |  |  |
| … | | | 3 further descriptor columns, named as in the source workbook |

### `data/processed/paired_mcc_dedup387_random.csv`

2520 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `repeat` | int64 |  |  |
| `cdft` | float64 |  |  |
| `physchem` | float64 |  |  |
| `combined` | float64 |  |  |

### `data/processed/paired_mcc_dedup387_scaffold.csv`

2520 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `repeat` | int64 |  |  |
| `cdft` | float64 |  |  |
| `physchem` | float64 |  |  |
| `combined` | float64 |  |  |

### `data/processed/paired_mcc_raw600_random.csv`

2520 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `strategy` | str |  |  |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `repeat` | int64 |  |  |
| `cdft` | float64 |  |  |
| `physchem` | float64 |  |  |
| `combined` | float64 |  |  |

## Result tables

### `results/tables/SI_T1_cdft_equation_audit.csv`

16 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `descriptor` | str |  |  |
| `equation_reproduced` | bool |  |  |
| `max_abs_deviation` | float64 |  |  |
| `pearson_r` | float64 |  |  |

### `results/tables/SI_T2_cdft_correlation_matrix.csv`

18 rows x 19 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Unnamed: 0` | str |  |  |
| `homo_energy` | float64 | hartree |  |
| `homo_.1_energy` | float64 | hartree |  |
| `lumo_energy` | float64 | hartree |  |
| `gap` | float64 | hartree |  |
| `potential_ionization` | float64 | hartree |  |
| `affinity_electronic` | float64 | hartree |  |
| `chemical_potential` | float64 | hartree |  |
| `hardness` | float64 | hartree |  |
| `softness` | float64 | hartree⁻¹ |  |
| `electronegativity` | float64 | hartree |  |
| `electrophilicity` | float64 | hartree |  |
| `nucleophilicity` | float64 | hartree |  |
| `acceptor_power` | float64 | hartree |  |
| `donor_power` | float64 | hartree |  |
| `net_electrophilicity` | float64 | hartree |  |
| `hyperhardness` | float64 | hartree |  |
| `chemical_potential_pos` | float64 | hartree |  |
| `chemical_potential_neg` | float64 | hartree |  |

### `results/tables/SI_T3_cdft_variance_from_frontier_orbitals.csv`

18 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `descriptor` | str |  |  |
| `R2_on_HOMO_LUMO_linear` | float64 |  |  |

### `results/tables/SI_T4_scaffold_population.csv`

85 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `murcko_scaffold` | str |  | Bemis-Murcko scaffold SMILES; the grouping variable for the scaffold split |
| `n_structures` | int64 |  |  |

### `results/tables/SI_T5_frontier_orbital_diagnostics.csv`

4 rows x 3 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `diagnostic` | str |  |  |
| `n` | int64 |  |  |
| `pct` | float64 |  |  |

### `results/tables/SI_T6_benchmark_summary_by_configuration.csv`

378 rows x 16 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `n_repeats` | int64 |  |  |
| `mcc_mean` | float64 |  |  |
| `mcc_ci95` | float64 |  |  |
| `bacc_mean` | float64 |  |  |
| `bacc_ci95` | float64 |  |  |
| `auc_mean` | float64 |  |  |
| `auc_ci95` | float64 |  |  |
| `acc_mean` | float64 |  |  |
| `sens_mean` | float64 |  |  |
| `spec_mean` | float64 |  |  |
| `n_selected` | float64 |  |  |

### `results/tables/SI_T7_shap_all_family_excluded.csv`

226 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Unnamed: 0` | str |  |  |
| `0` | float64 |  |  |

### `results/tables/SI_T7_shap_all_family_included.csv`

230 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Unnamed: 0` | str |  |  |
| `0` | float64 |  |  |

### `results/tables/SI_T8_feature_selection_stability.csv`

312 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `method` | str |  |  |
| `feature` | str |  |  |
| `origin` | str |  |  |
| `selection_frequency` | float64 |  |  |
| `n_runs` | int64 |  |  |

### `results/tables/SI_cdft_equation_audit.csv`

16 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `descriptor` | str |  |  |
| `equation_reproduced` | bool |  |  |
| `max_abs_deviation` | float64 |  |  |
| `pearson_r` | float64 |  |  |

### `results/tables/SI_cdft_from_frontier_orbitals.csv`

18 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `descriptor` | str |  |  |
| `R2_on_HOMO_LUMO_linear` | float64 |  |  |

### `results/tables/SI_descriptor_list.csv`

306 rows x 2 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `descriptor` | str |  |  |
| `block` | str |  | descriptor block: physchem, cdft or combined |

### `results/tables/SI_effective_dimensionality.csv`

3 rows x 9 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Descriptor block` | str |  |  |
| `n_variables` | float64 |  |  |
| `n_PC_50` | float64 |  |  |
| `n_PC_90` | float64 |  |  |
| `n_PC_95` | float64 |  |  |
| `n_PC_99` | float64 |  |  |
| `var_PC1` | float64 |  |  |
| `var_PC1_PC2` | float64 |  |  |
| `participation_ratio` | float64 |  |  |

### `results/tables/SI_family_composition.csv`

4 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `tipo_familia` | str |  | chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros) |
| `n` | float64 |  |  |
| `n_censored` | float64 |  |  |
| `median_IC50` | float64 |  |  |
| `pct_below_10uM` | float64 |  |  |

### `results/tables/SI_family_vs_clusters.csv`

7 rows x 3 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Clusters` | float64 |  |  |
| `Adjusted Rand index` | float64 |  |  |
| `Adjusted mutual information` | float64 |  |  |

### `results/tables/SI_feature_stability.csv`

312 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `method` | str |  |  |
| `feature` | str |  |  |
| `origin` | str |  |  |
| `selection_frequency` | float64 |  |  |
| `n_runs` | float64 |  |  |

### `results/tables/SI_orbital_diagnostics.csv`

4 rows x 3 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `diagnostic` | str |  |  |
| `n` | float64 |  |  |
| `pct` | float64 |  |  |

### `results/tables/SI_original_R_contrasts.csv`

3 rows x 9 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `comparison` | str |  |  |
| `median_delta` | float64 |  |  |
| `mean_delta` | float64 |  |  |
| `q1` | float64 |  |  |
| `q3` | float64 |  |  |
| `n_pairs` | float64 |  |  |
| `n_favoring_second` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |
| `nemenyi_p` | float64 |  |  |

### `results/tables/SI_shap_attribution.csv`

5 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `origin` | str |  |  |
| `n_models` | float64 |  |  |
| `mean_share_pct` | float64 |  |  |
| `ci95_pct` | float64 |  |  |

### `results/tables/T10_roc_auc.csv`

6 rows x 6 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Split scheme` | str |  |  |
| `Descriptor block` | str |  |  |
| `Mean AUC (dimensionless)` | float64 |  |  |
| `95% CI half-width (dimensionless)` | float64 |  |  |
| `Minimum over repeats (dimensionless)` | float64 |  |  |
| `Maximum over repeats (dimensionless)` | float64 |  |  |

### `results/tables/T1_block_composition.csv`

3 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Descriptor block` | str |  |  |
| `Descriptors computed (n)` | int64 |  |  |
| `Constant or near-constant (n)` | int64 |  |  |
| `Informative on this set (n)` | int64 |  |  |
| `Source` | str |  |  |

### `results/tables/T2_structural_redundancy.csv`

4 rows x 7 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Level of structural identity` | str |  |  |
| `Distinct entities (n)` | int64 |  |  |
| `Groups with >1 record (n)` | int64 |  |  |
| `Records in such groups (n)` | int64 |  |  |
| `Groups with discordant IC50 (n)` | int64 |  |  |
| `Groups mixing censored and measured (n)` | int64 |  |  |
| `Groups with discordant class at 10 uM (n)` | int64 |  |  |

### `results/tables/T3_class_composition.csv`

12 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Compound set` | str |  |  |
| `Activity definition` | str |  |  |
| `Type` | str |  |  |
| `Inhibitors, UI (n)` | int64 |  |  |
| `Non-inhibitors, nUI (n)` | int64 |  |  |
| `Total modeled (n)` | int64 |  |  |
| `Excluded, intermediate zone (n)` | int64 |  |  |
| `Minority class (%)` | float64 |  |  |

### `results/tables/T4_partition_sizes.csv`

12 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Compound set` | str |  |  |
| `Activity definition` | str |  |  |
| `Modeled compounds (n)` | int64 |  |  |
| `Training, 80% (n)` | int64 |  |  |
| `Test, 20% (n)` | int64 |  |  |

### `results/tables/T5_classifier_configuration.csv`

7 rows x 3 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Classifier` | str |  |  |
| `Implementation` | str |  |  |
| `Setting (identical for every descriptor block)` | str |  |  |

### `results/tables/T6_performance_by_arm.csv`

18 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Evaluation protocol` | str |  |  |
| `Descriptor block` | str |  |  |
| `Chemical family` | str |  |  |
| `Models (n)` | int64 |  |  |
| `MCC median (dimensionless)` | float64 |  |  |
| `MCC Q1 (dimensionless)` | float64 |  |  |
| `MCC Q3 (dimensionless)` | float64 |  |  |
| `MCC max (dimensionless)` | float64 |  |  |
| `Balanced accuracy median (dimensionless)` | float64 |  |  |
| `AUC median (dimensionless)` | float64 |  |  |

### `results/tables/T7a_friedman.csv`

3 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Evaluation protocol` | str |  |  |
| `Matched configurations (n)` | int64 |  |  |
| `Friedman chi-square (dimensionless)` | float64 |  |  |
| `p (dimensionless)` | float64 |  |  |
| `Mean rank CDFT` | float64 |  |  |
| `Mean rank physicochemical` | float64 |  |  |
| `Mean rank combined` | float64 |  |  |
| `Median MCC CDFT (dimensionless)` | float64 |  |  |
| `Median MCC physicochemical (dimensionless)` | float64 |  |  |
| `Median MCC combined (dimensionless)` | float64 |  |  |

### `results/tables/T7b_block_contrasts.csv`

9 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Evaluation protocol` | str |  |  |
| `Contrast` | str |  |  |
| `Pairs (n)` | int64 |  |  |
| `Median dMCC (dimensionless)` | float64 |  |  |
| `Mean dMCC (dimensionless)` | float64 |  |  |
| `95% CI half-width (dimensionless)` | float64 |  |  |
| `% favoring the first term` | float64 |  |  |
| `Rank-biserial r (dimensionless)` | float64 |  |  |
| `Wilcoxon p (dimensionless)` | float64 |  |  |
| `Nemenyi p (dimensionless)` | float64 |  |  |

### `results/tables/T7c_block_contrasts_family_excluded.csv`

9 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Evaluation protocol` | str |  |  |
| `Contrast` | str |  |  |
| `Pairs (n)` | int64 |  |  |
| `Median dMCC (dimensionless)` | float64 |  |  |
| `Mean dMCC (dimensionless)` | float64 |  |  |
| `95% CI half-width (dimensionless)` | float64 |  |  |
| `% favoring the first term` | float64 |  |  |
| `Rank-biserial r (dimensionless)` | float64 |  |  |
| `Wilcoxon p (dimensionless)` | float64 |  |  |
| `Nemenyi p (dimensionless)` | float64 |  |  |

### `results/tables/T8_protocol_effect.csv`

12 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Descriptor block` | str |  |  |
| `Chemical family` | str |  |  |
| `Configurations (n)` | int64 |  |  |
| `Mean MCC, first protocol (dimensionless)` | float64 |  |  |
| `Mean MCC, second protocol (dimensionless)` | float64 |  |  |
| `Mean change (dimensionless)` | float64 |  |  |
| `95% CI half-width (dimensionless)` | float64 |  |  |
| `Relative change (%)` | float64 |  |  |
| `Wilcoxon p (dimensionless)` | float64 |  |  |
| `Contrast` | str |  |  |

### `results/tables/T8b_family_attribute.csv`

9 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Evaluation protocol` | str |  |  |
| `Descriptor block` | str |  |  |
| `Pairs (n)` | int64 |  |  |
| `Median dMCC (dimensionless)` | float64 |  |  |
| `Mean dMCC (dimensionless)` | float64 |  |  |
| `95% CI half-width (dimensionless)` | float64 |  |  |
| `% favoring inclusion` | float64 |  |  |
| `Wilcoxon p (dimensionless)` | float64 |  |  |

### `results/tables/T9_best_configurations.csv`

15 rows x 11 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Descriptor block` | str |  |  |
| `Chemical family` | str |  |  |
| `Feature selection` | str |  |  |
| `Classifier` | str |  |  |
| `MCC (dimensionless)` | float64 |  |  |
| `95% CI (dimensionless)` | float64 |  |  |
| `Balanced accuracy (dimensionless)` | float64 |  |  |
| `AUC (dimensionless)` | float64 |  |  |
| `Sensitivity (dimensionless)` | float64 |  |  |
| `Specificity (dimensionless)` | float64 |  |  |
| `Descriptors retained (n)` | int64 |  |  |

### `results/tables/T_benchmark_by_arm.csv`

18 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `n_models` | int64 |  |  |
| `mcc_median` | float64 |  |  |
| `mcc_q1` | float64 |  |  |
| `mcc_q3` | float64 |  |  |
| `mcc_max` | float64 |  |  |
| `bacc_median` | float64 |  |  |
| `auc_median` | float64 |  |  |

### `results/tables/T_block_pairwise_by_regime.csv`

9 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `comparison` | str |  |  |
| `n_pairs` | int64 |  |  |
| `median_delta` | float64 |  |  |
| `mean_delta` | float64 |  |  |
| `ci95_delta` | float64 |  |  |
| `pct_favoring_second` | float64 |  |  |
| `rank_biserial` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |
| `nemenyi_p` | float64 |  |  |

### `results/tables/T_block_pairwise_family_excluded.csv`

9 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `comparison` | str |  |  |
| `n_pairs` | int64 |  |  |
| `median_delta` | float64 |  |  |
| `mean_delta` | float64 |  |  |
| `ci95_delta` | float64 |  |  |
| `pct_favoring_second` | float64 |  |  |
| `rank_biserial` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |
| `nemenyi_p` | float64 |  |  |

### `results/tables/T_class_composition.csv`

12 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `dataset` | str |  |  |
| `strategy` | str |  |  |
| `kind` | str |  |  |
| `n_UI` | int64 |  |  |
| `n_nUI` | int64 |  |  |
| `n_total` | int64 |  |  |
| `n_excluded` | int64 |  |  |
| `pct_minority` | float64 |  |  |

### `results/tables/T_effective_dimensionality.csv`

3 rows x 9 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `n_variables` | int64 |  |  |
| `n_PC_50` | int64 |  |  |
| `n_PC_90` | int64 |  |  |
| `n_PC_95` | int64 |  |  |
| `n_PC_99` | int64 |  |  |
| `var_PC1` | float64 |  |  |
| `var_PC1_PC2` | float64 |  |  |
| `participation_ratio` | float64 |  |  |

### `results/tables/T_family_attribute_effect.csv`

9 rows x 8 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `n_pairs` | int64 |  |  |
| `median_delta_mcc` | float64 |  |  |
| `mean_delta_mcc` | float64 |  |  |
| `ci95` | float64 |  |  |
| `pct_favoring_inclusion` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |

### `results/tables/T_family_composition.csv`

4 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `tipo_familia` | str |  | chemical family as recorded in the source compilation (Coumarin, Triazol, Tiazol, Otros) |
| `n` | int64 |  |  |
| `n_censored` | int64 |  |  |
| `median_IC50` | float64 |  |  |
| `pct_below_10uM` | float64 |  |  |

### `results/tables/T_family_vs_fingerprint_clusters.csv`

7 rows x 3 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `n_clusters` | int64 |  |  |
| `adjusted_rand_index` | float64 |  |  |
| `adjusted_mutual_information` | float64 |  |  |

### `results/tables/T_friedman_by_regime.csv`

3 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `regime` | str |  | compound set and split scheme, one of raw600/random, dedup387/random, dedup387/scaffold |
| `n_matched_configurations` | int64 |  |  |
| `friedman_chi2` | float64 |  |  |
| `friedman_p` | float64 |  |  |
| `mean_rank_cdft` | float64 |  |  |
| `mean_rank_physchem` | float64 |  |  |
| `mean_rank_combined` | float64 |  |  |
| `median_mcc_cdft` | float64 |  |  |
| `median_mcc_physchem` | float64 |  |  |
| `median_mcc_combined` | float64 |  |  |

### `results/tables/T_original_R_block_comparison.csv`

3 rows x 9 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `comparison` | str |  |  |
| `median_delta` | float64 |  |  |
| `mean_delta` | float64 |  |  |
| `q1` | float64 |  |  |
| `q3` | float64 |  |  |
| `n_pairs` | int64 |  |  |
| `n_favoring_second` | int64 |  |  |
| `wilcoxon_p` | float64 |  |  |
| `nemenyi_p` | float64 |  |  |

### `results/tables/T_original_R_family_effect.csv`

3 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `n_pairs` | int64 |  |  |
| `median_delta_mcc` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |

### `results/tables/T_original_R_maxima_vs_median.csv`

6 rows x 6 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `n` | int64 |  |  |
| `median` | float64 |  |  |
| `q90` | float64 |  |  |
| `max` | float64 |  |  |

### `results/tables/T_roc_auc.csv`

6 rows x 6 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `scheme` | str |  |  |
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `auc_mean` | float64 |  |  |
| `auc_ci95` | float64 |  |  |
| `auc_min` | float64 |  |  |
| `auc_max` | float64 |  |  |

### `results/tables/T_shap_attribution_by_origin.csv`

5 rows x 5 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `origin` | str |  |  |
| `n_models` | int64 |  |  |
| `mean_share_pct` | float64 |  |  |
| `ci95_pct` | float64 |  |  |

### `results/tables/T_shap_top25_family_excluded.csv`

25 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Unnamed: 0` | str |  |  |
| `mean_shap_share` | float64 |  |  |
| `origin` | str |  |  |
| `rank` | int64 |  |  |

### `results/tables/T_shap_top25_family_included.csv`

25 rows x 4 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `Unnamed: 0` | str |  |  |
| `mean_shap_share` | float64 |  |  |
| `origin` | str |  |  |
| `rank` | int64 |  |  |

### `results/tables/T_split_scheme_effect.csv`

12 rows x 10 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `n_configurations` | int64 |  |  |
| `mcc_random` | float64 |  |  |
| `mcc_scaffold` | float64 |  |  |
| `mean_delta` | float64 |  |  |
| `ci95` | float64 |  |  |
| `pct_relative_drop` | float64 |  |  |
| `wilcoxon_p` | float64 |  |  |
| `contrast` | str |  |  |

### `results/tables/T_top_configurations_scaffold.csv`

15 rows x 11 columns.

| column | dtype | unit | notes |
|---|---|---|---|
| `block` | str |  | descriptor block: physchem, cdft or combined |
| `family` | str |  | whether the chemical-family attribute was supplied to the model (included) or withheld (excluded) |
| `selection` | str |  |  |
| `algorithm` | str |  |  |
| `mcc_mean` | float64 |  |  |
| `mcc_ci95` | float64 |  |  |
| `bacc_mean` | float64 |  |  |
| `auc_mean` | float64 |  |  |
| `sens_mean` | float64 |  |  |
| `spec_mean` | float64 |  |  |
| `n_selected` | float64 |  |  |
