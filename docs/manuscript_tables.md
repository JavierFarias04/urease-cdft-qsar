
## Tables

**Table 1.** Composition of the three descriptor blocks. The blocks were built on the same 387 standardized structures, so descriptor content is the only factor that differs between them. Constant or near-constant descriptors were removed on the training partition by the criterion in Section 2.7; the CDFT column removed is a dual reactivity index that takes the value 1 for every compound.

| Descriptor block | Descriptors computed (n) | Constant or near-constant (n) | Informative on this set (n) | Source |
|---|---|---|---|---|
| Physicochemical and topological (CDK) | 287 | 61 | 226 | rcdk / CDK, 2D structure |
| CDFT global reactivity (Koopmans) | 19 | 1 | 18 | B3LYP-D3/6-31G frontier orbital energies |
| Combined | 306 | 62 | 244 | union of the two blocks |


**Table 2.** Structural redundancy of the 600 analyzed records at four levels of identity. Discordance is evaluated within groups sharing the same entity: 'mixing censored and measured' counts groups in which one record reported no detectable inhibition while another reported a measured IC50, so that identical structures carry opposite activity labels.

| Level of structural identity | Distinct entities (n) | Groups with >1 record (n) | Records in such groups (n) | Groups with discordant IC50 (n) | Groups mixing censored and measured (n) | Groups with discordant class at 10 uM (n) |
|---|---|---|---|---|---|---|
| canonical SMILES (standardized) | 475 | 114 | 239 | 23 | 20 | 20 |
| InChIKey (full) | 387 | 157 | 370 | 16 | 12 | 11 |
| InChIKey connectivity block | 331 | 155 | 424 | 16 | 12 | 11 |
| Bemis-Murcko scaffold | 86 | 65 | 579 | 45 | 15 | 12 |


**Table 3.** Class composition of the six activity definitions on the 387 standardized structures. Under the two-cutoff definitions the non-inhibitor class is defined invariantly as IC50 > 50 uM, so its size does not vary with the lower threshold; compounds between the two thresholds are excluded. The corresponding composition of the 600-record set is given in the Supporting Information. UI, urease inhibitor; nUI, non-inhibitor.

| Activity definition | Type | Inhibitors, UI (n) | Non-inhibitors, nUI (n) | Total modeled (n) | Excluded, intermediate zone (n) | Minority class (%) |
|---|---|---|---|---|---|---|
| 5 uM | single-cutoff | 87 | 300 | 387 | 0 | 22.5 |
| 5 and 50 uM | two-cutoff | 87 | 184 | 271 | 116 | 32.1 |
| 10 uM | single-cutoff | 97 | 290 | 387 | 0 | 25.1 |
| 10 and 50 uM | two-cutoff | 97 | 184 | 281 | 106 | 34.5 |
| 25 uM | single-cutoff | 129 | 258 | 387 | 0 | 33.3 |
| 25 and 50 uM | two-cutoff | 129 | 184 | 313 | 74 | 41.2 |


**Table 4.** Classifier configuration. Hyperparameters were fixed and identical for every descriptor block, so that unequal tuning effort could not confound the block comparison (Section 2.11). p denotes the number of predictors after selection.

| Classifier | Implementation | Setting (identical for every descriptor block) |
|---|---|---|
| Random forest | scikit-learn RandomForestClassifier | 500 trees, sqrt(p) features per split, unlimited depth |
| SVM (RBF kernel) | scikit-learn SVC | C = 1.0, gamma = 1/(p var(X)), probability estimates enabled |
| Decision tree (CART) | scikit-learn DecisionTreeClassifier | Gini impurity, unlimited depth |
| XGBoost | xgboost XGBClassifier | 200 rounds, max_depth 4, eta 0.1, subsample 0.9, colsample 0.8, lambda 1 |
| k-nearest neighbors | scikit-learn KNeighborsClassifier | k = 5, Euclidean distance on standardized inputs |
| Gaussian naive Bayes | scikit-learn GaussianNB | default variance smoothing 1e-9 |
| Logistic regression | scikit-learn LogisticRegression | L2 penalty, C = 1.0, lbfgs, 3000 iterations |


**Table 5.** Performance of every arm under the three evaluation protocols, over the full population of matched configurations (7 classifiers x 6 activity definitions x 3 selection regimes x 10 repeats). Maxima are shown for comparison with the literature but are not used as summary statistics (Section 2.12).

| Evaluation protocol | Descriptor block | Chemical family | Models | MCC median | MCC Q1 | MCC Q3 | MCC max | Bal. acc. median | AUC median |
|---|---|---|---|---|---|---|---|---|---|
| 387 structures, random split | CDFT reactivity | excluded | 1260 | 0.508 | 0.411 | 0.597 | 0.845 | 0.753 | 0.809 |
| 387 structures, random split | CDFT reactivity | included | 1260 | 0.571 | 0.471 | 0.669 | 0.886 | 0.786 | 0.853 |
| 387 structures, random split | Combined | excluded | 1260 | 0.738 | 0.628 | 0.825 | 1.000 | 0.866 | 0.927 |
| 387 structures, random split | Combined | included | 1260 | 0.738 | 0.634 | 0.830 | 1.000 | 0.867 | 0.928 |
| 387 structures, random split | Physicochemical | excluded | 1260 | 0.736 | 0.635 | 0.824 | 1.000 | 0.866 | 0.927 |
| 387 structures, random split | Physicochemical | included | 1260 | 0.736 | 0.632 | 0.823 | 1.000 | 0.866 | 0.927 |
| 387 structures, scaffold split | CDFT reactivity | excluded | 1260 | 0.292 | 0.106 | 0.478 | 0.868 | 0.648 | 0.713 |
| 387 structures, scaffold split | CDFT reactivity | included | 1260 | 0.357 | 0.175 | 0.527 | 1.000 | 0.676 | 0.747 |
| 387 structures, scaffold split | Combined | excluded | 1260 | 0.537 | 0.354 | 0.702 | 1.000 | 0.779 | 0.864 |
| 387 structures, scaffold split | Combined | included | 1260 | 0.535 | 0.347 | 0.700 | 1.000 | 0.774 | 0.860 |
| 387 structures, scaffold split | Physicochemical | excluded | 1260 | 0.554 | 0.357 | 0.719 | 1.000 | 0.782 | 0.864 |
| 387 structures, scaffold split | Physicochemical | included | 1260 | 0.560 | 0.355 | 0.717 | 1.000 | 0.784 | 0.863 |
| 600 records, random split | CDFT reactivity | excluded | 1260 | 0.332 | 0.254 | 0.421 | 0.760 | 0.677 | 0.751 |
| 600 records, random split | CDFT reactivity | included | 1260 | 0.588 | 0.512 | 0.663 | 0.862 | 0.811 | 0.876 |
| 600 records, random split | Combined | excluded | 1260 | 0.634 | 0.557 | 0.698 | 0.848 | 0.825 | 0.905 |
| 600 records, random split | Combined | included | 1260 | 0.782 | 0.664 | 0.866 | 1.000 | 0.895 | 0.948 |
| 600 records, random split | Physicochemical | excluded | 1260 | 0.640 | 0.565 | 0.706 | 0.862 | 0.827 | 0.910 |
| 600 records, random split | Physicochemical | included | 1260 | 0.798 | 0.685 | 0.871 | 1.000 | 0.901 | 0.952 |


**Table 6.** Friedman test across the three descriptor blocks, paired on the tuple (chemical-family condition, activity definition, feature-selection regime, classifier, repeat). A mean rank of 1 is best.

| Evaluation protocol | Matched configurations | Friedman chi-square | p | Mean rank CDFT | Mean rank physicochemical | Mean rank combined | Median MCC CDFT | Median MCC physicochemical | Median MCC combined |
|---|---|---|---|---|---|---|---|---|---|
| 387 structures, random split | 2520 | 2796.2 | < 1e-300 | 2.830 | 1.600 | 1.570 | 0.537 | 0.736 | 0.738 |
| 387 structures, scaffold split | 2520 | 1059.6 | 8.1 × 10^-231^ | 2.510 | 1.700 | 1.790 | 0.329 | 0.556 | 0.536 |
| 600 records, random split | 2520 | 3258.7 | < 1e-300 | 2.900 | 1.480 | 1.620 | 0.472 | 0.696 | 0.684 |


**Table 7.** Pairwise contrasts between descriptor blocks in MCC, paired as in Table 6. A positive value favors the first term of the contrast. The combined-minus-physicochemical contrast is the test of the central hypothesis: its median is zero under every protocol and its mean is negative under scaffold splitting. Nemenyi p values are from the post hoc test following Friedman.

| Evaluation protocol | Contrast | Pairs | Median dMCC | Mean dMCC | 95% CI | % favoring first | Rank-biserial r | Wilcoxon p | Nemenyi p |
|---|---|---|---|---|---|---|---|---|---|
| 387 structures, random split | Physicochemical - CDFT reactivity | 2520 | 0.183 | 0.189 | 0.006 | 91.3 | 0.833 | < 1e-300 | < 1e-300 |
| 387 structures, random split | Combined - CDFT reactivity | 2520 | 0.188 | 0.191 | 0.006 | 91.3 | 0.831 | < 1e-300 | < 1e-300 |
| 387 structures, random split | Combined - Physicochemical | 2520 | 0.000 | 0.002 | 0.002 | 39.2 | 0.028 | 0.055 | 0.604 |
| 387 structures, scaffold split | Physicochemical - CDFT reactivity | 2520 | 0.206 | 0.212 | 0.013 | 75.6 | 0.519 | 7.8 × 10^-173^ | < 1e-300 |
| 387 structures, scaffold split | Combined - CDFT reactivity | 2520 | 0.192 | 0.200 | 0.013 | 74.7 | 0.502 | 3.9 × 10^-165^ | < 1e-300 |
| 387 structures, scaffold split | Combined - Physicochemical | 2520 | 0.000 | -0.012 | 0.005 | 34.5 | -0.075 | 2.1 × 10^-7^ | 0.009 |
| 600 records, random split | Physicochemical - CDFT reactivity | 2520 | 0.234 | 0.235 | 0.005 | 94.6 | 0.896 | < 1e-300 | < 1e-300 |
| 600 records, random split | Combined - CDFT reactivity | 2520 | 0.224 | 0.227 | 0.005 | 94.7 | 0.898 | < 1e-300 | < 1e-300 |
| 600 records, random split | Combined - Physicochemical | 2520 | 0.000 | -0.009 | 0.002 | 31.5 | -0.145 | 6.9 × 10^-21^ | 9.8 × 10^-7^ |


**Table 8.** Effect of the evaluation protocol on apparent performance. Each row pairs the same (activity definition, selection regime, classifier) configuration between two protocols. A positive relative change denotes a loss of performance in the second protocol. Note that removing duplicate records raises performance when the chemical family is withheld, because it eliminates records carrying contradictory activity labels.

| Contrast | Descriptor block | Chemical family | Configurations | Mean MCC, first | Mean MCC, second | Mean change | 95% CI | Relative change (%) | Wilcoxon p |
|---|---|---|---|---|---|---|---|---|---|
| scaffold vs random split (387 structures) | CDFT reactivity | included | 126 | 0.563 | 0.342 | -0.221 | 0.020 | 39.3 | 5.7 × 10^-22^ |
| scaffold vs random split (387 structures) | CDFT reactivity | excluded | 126 | 0.492 | 0.284 | -0.208 | 0.017 | 42.2 | 2.5 × 10^-22^ |
| scaffold vs random split (387 structures) | Physicochemical | included | 126 | 0.716 | 0.525 | -0.192 | 0.018 | 26.7 | 5.5 × 10^-22^ |
| scaffold vs random split (387 structures) | Physicochemical | excluded | 126 | 0.717 | 0.526 | -0.190 | 0.018 | 26.6 | 6.4 × 10^-22^ |
| scaffold vs random split (387 structures) | Combined | included | 126 | 0.719 | 0.510 | -0.209 | 0.018 | 29.0 | 4.5 × 10^-22^ |
| scaffold vs random split (387 structures) | Combined | excluded | 126 | 0.718 | 0.517 | -0.201 | 0.018 | 28.0 | 5.4 × 10^-22^ |
| duplicate removal (600 -> 387, both random) | CDFT reactivity | included | 126 | 0.585 | 0.563 | -0.022 | 0.013 | 3.8 | 0.003 |
| duplicate removal (600 -> 387, both random) | CDFT reactivity | excluded | 126 | 0.334 | 0.492 | 0.158 | 0.016 | -47.3 | 4.0 × 10^-22^ |
| duplicate removal (600 -> 387, both random) | Physicochemical | included | 126 | 0.761 | 0.716 | -0.045 | 0.010 | 5.9 | 6.8 × 10^-13^ |
| duplicate removal (600 -> 387, both random) | Physicochemical | excluded | 126 | 0.629 | 0.717 | 0.088 | 0.013 | -14.0 | 2.4 × 10^-21^ |
| duplicate removal (600 -> 387, both random) | Combined | included | 126 | 0.752 | 0.719 | -0.033 | 0.009 | 4.4 | 1.6 × 10^-9^ |
| duplicate removal (600 -> 387, both random) | Combined | excluded | 126 | 0.620 | 0.718 | 0.098 | 0.012 | -15.8 | 2.6 × 10^-22^ |


**Table 9.** The ten best-performing configurations under scaffold-based partitioning, ranked by mean MCC over 10 repeats and six activity definitions, with the 95% confidence interval of the mean. The leading configurations from the combined and physicochemical blocks are statistically indistinguishable. Descriptors retained is the mean number surviving feature selection.

| Descriptor block | Chemical family | Feature selection | Classifier | MCC | 95% CI | Balanced accuracy | AUC | Sensitivity | Specificity | Descriptors retained |
|---|---|---|---|---|---|---|---|---|---|---|
| Combined | excluded | XGBoost gain | XGBoost | 0.599 | 0.050 | 0.804 | 0.879 | 0.788 | 0.820 | 76 |
| Physicochemical | included | Boruta | k-nearest neighbors | 0.591 | 0.064 | 0.803 | 0.863 | 0.807 | 0.799 | 26 |
| Combined | included | XGBoost gain | XGBoost | 0.591 | 0.053 | 0.799 | 0.879 | 0.779 | 0.820 | 76 |
| Physicochemical | included | Boruta | Random forest | 0.590 | 0.061 | 0.797 | 0.878 | 0.778 | 0.816 | 26 |
| Physicochemical | excluded | Boruta | k-nearest neighbors | 0.587 | 0.065 | 0.800 | 0.857 | 0.808 | 0.792 | 26 |
| Physicochemical | excluded | Boruta | Random forest | 0.585 | 0.063 | 0.793 | 0.876 | 0.780 | 0.806 | 26 |
| Physicochemical | included | Boruta | XGBoost | 0.584 | 0.064 | 0.796 | 0.858 | 0.777 | 0.814 | 26 |
| Physicochemical | included | XGBoost gain | XGBoost | 0.579 | 0.052 | 0.796 | 0.871 | 0.772 | 0.821 | 72 |
| Physicochemical | excluded | XGBoost gain | XGBoost | 0.578 | 0.055 | 0.795 | 0.874 | 0.770 | 0.819 | 72 |
| Physicochemical | excluded | none | Logistic regression | 0.578 | 0.054 | 0.796 | 0.865 | 0.785 | 0.808 | 205 |


**Table 10.** Effect of supplying the chemical family as an explicit attribute, paired on (activity definition, selection regime, classifier, repeat). The substantial benefit seen on the 600 records falls to zero for the physicochemical and combined blocks once structures are deduplicated.

| Evaluation protocol | Descriptor block | Pairs | Median dMCC | Mean dMCC | 95% CI | % favoring inclusion | Wilcoxon p |
|---|---|---|---|---|---|---|---|
| 387 structures, random split | CDFT reactivity | 1260 | 0.064 | 0.071 | 0.006 | 67.6 | 1.9 × 10^-92^ |
| 387 structures, random split | Combined | 1260 | 0.000 | 0.001 | 0.002 | 21.9 | 0.932 |
| 387 structures, random split | Physicochemical | 1260 | 0.000 | -0.000 | 0.002 | 20.3 | 0.306 |
| 387 structures, scaffold split | CDFT reactivity | 1260 | 0.031 | 0.057 | 0.013 | 54.0 | 1.5 × 10^-21^ |
| 387 structures, scaffold split | Combined | 1260 | 0.000 | -0.007 | 0.005 | 23.6 | 0.228 |
| 387 structures, scaffold split | Physicochemical | 1260 | 0.000 | -0.002 | 0.006 | 24.0 | 0.748 |
| 600 records, random split | CDFT reactivity | 1260 | 0.249 | 0.251 | 0.007 | 97.1 | 1.4 × 10^-204^ |
| 600 records, random split | Combined | 1260 | 0.118 | 0.132 | 0.006 | 92.1 | 3.2 × 10^-192^ |
| 600 records, random split | Physicochemical | 1260 | 0.121 | 0.133 | 0.006 | 91.6 | 1.1 × 10^-191^ |
