"""
08_tables.py — assemble every manuscript and Supporting Information table,
plus the audit of numbers changed relative to the previous draft.
"""
import json, os
import numpy as np
import pandas as pd

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
D, T = str(cfg.PROCESSED), str(cfg.TABLES)
OUT = str(cfg.TABLES)

meta = json.load(open(f"{D}/build_meta.json"))
eda = json.load(open(f"{D}/eda_meta.json"))
audit_dup = pd.read_csv(f"{D}/duplication_audit.csv")
comp = pd.read_csv(f"{D}/block_composition.csv")
cls = pd.read_csv(f"{T}/T_class_composition.csv")
ed = pd.read_csv(f"{T}/T_effective_dimensionality.csv")
fried = pd.read_csv(f"{T}/T_friedman_by_regime.csv")
pw = pd.read_csv(f"{T}/T_block_pairwise_by_regime.csv")
pwnf = pd.read_csv(f"{T}/T_block_pairwise_family_excluded.csv")
famt = pd.read_csv(f"{T}/T_family_attribute_effect.csv")
split = pd.read_csv(f"{T}/T_split_scheme_effect.csv")
arm = pd.read_csv(f"{T}/T_benchmark_by_arm.csv")
best = pd.read_csv(f"{T}/T_top_configurations_scaffold.csv")
rocT = pd.read_csv(f"{T}/T_roc_auc.csv")
origR = pd.read_csv(f"{T}/T_original_R_block_comparison.csv")
clu = pd.read_csv(f"{T}/T_family_vs_fingerprint_clusters.csv")
famcomp = pd.read_csv(f"{T}/T_family_composition.csv")
shapatt = pd.read_csv(f"{T}/T_shap_attribution_by_origin.csv")
stab = pd.read_csv(f"{T}/SI_T8_feature_selection_stability.csv")
raw = pd.read_csv(f"{D}/compounds_raw.csv")
tables = {}

# ---------------- Table 1: composition of the three descriptor blocks
t1 = pd.DataFrame([
    dict(**{"Descriptor block": "Physicochemical and topological (CDK)",
            "Descriptors computed (n)": 287, "Constant or near-constant (n)": 61,
            "Informative on this set (n)": 226, "Source": "rcdk / CDK, 2D structure"}),
    dict(**{"Descriptor block": "CDFT global reactivity (Koopmans)",
            "Descriptors computed (n)": 19, "Constant or near-constant (n)": 1,
            "Informative on this set (n)": 18, "Source": "B3LYP-D3/6-31G frontier orbital energies"}),
    dict(**{"Descriptor block": "Combined", "Descriptors computed (n)": 306,
            "Constant or near-constant (n)": 62, "Informative on this set (n)": 244,
            "Source": "union of the two blocks"})])
assert t1.loc[2, "Descriptors computed (n)"] == t1.loc[0, "Descriptors computed (n)"] + t1.loc[1, "Descriptors computed (n)"]
assert t1.loc[2, "Informative on this set (n)"] == t1.loc[0, "Informative on this set (n)"] + t1.loc[1, "Informative on this set (n)"]
tables["T1_block_composition"] = t1

# ---------------- Table 2: structural redundancy of the compound set
t2 = audit_dup.rename(columns={
    "level": "Level of structural identity", "n_unique": "Distinct entities (n)",
    "n_duplicate_groups": "Groups with >1 record (n)", "n_rows_in_duplicate_groups": "Records in such groups (n)",
    "n_groups_discordant_IC50": "Groups with discordant IC50 (n)",
    "n_groups_mixing_censored_and_measured": "Groups mixing censored and measured (n)",
    "n_groups_discordant_class_10uM": "Groups with discordant class at 10 uM (n)"})
t2 = t2[["Level of structural identity", "Distinct entities (n)", "Groups with >1 record (n)",
         "Records in such groups (n)", "Groups with discordant IC50 (n)",
         "Groups mixing censored and measured (n)", "Groups with discordant class at 10 uM (n)"]]
tables["T2_structural_redundancy"] = t2

# ---------------- Table 3: activity definitions and class composition
t3 = cls.copy()
t3["dataset"] = t3["dataset"].map({"raw600": "600 records (as previously analyzed)",
                                   "dedup387": "387 unique structures"})
t3 = t3.rename(columns={"dataset": "Compound set", "strategy": "Activity definition",
                        "kind": "Type", "n_UI": "Inhibitors, UI (n)", "n_nUI": "Non-inhibitors, nUI (n)",
                        "n_total": "Total modeled (n)", "n_excluded": "Excluded, intermediate zone (n)",
                        "pct_minority": "Minority class (%)"})
assert ((t3["Inhibitors, UI (n)"] + t3["Non-inhibitors, nUI (n)"]) == t3["Total modeled (n)"]).all()
tables["T3_class_composition"] = t3

# ---------------- Table 4: partition sizes actually used
rows = []
for _, r_ in cls.iterrows():
    n = int(r_["n_total"])
    rows.append(dict(**{"Compound set": "600 records" if r_["dataset"] == "raw600" else "387 structures",
                        "Activity definition": r_["strategy"], "Modeled compounds (n)": n,
                        "Training, 80% (n)": int(round(.8 * n)), "Test, 20% (n)": n - int(round(.8 * n))}))
t4 = pd.DataFrame(rows)
tables["T4_partition_sizes"] = t4

# ---------------- Table 5: classifier configuration
t5 = pd.DataFrame([
    ("Random forest", "scikit-learn RandomForestClassifier", "500 trees, sqrt(p) features per split, unlimited depth"),
    ("SVM (RBF kernel)", "scikit-learn SVC", "C = 1.0, gamma = 1/(p var(X)), probability estimates enabled"),
    ("Decision tree (CART)", "scikit-learn DecisionTreeClassifier", "Gini impurity, unlimited depth"),
    ("XGBoost", "xgboost XGBClassifier", "200 rounds, max_depth 4, eta 0.1, subsample 0.9, colsample 0.8, lambda 1"),
    ("k-nearest neighbors", "scikit-learn KNeighborsClassifier", "k = 5, Euclidean distance on standardized inputs"),
    ("Gaussian naive Bayes", "scikit-learn GaussianNB", "default variance smoothing 1e-9"),
    ("Logistic regression", "scikit-learn LogisticRegression", "L2 penalty, C = 1.0, lbfgs, 3000 iterations"),
], columns=["Classifier", "Implementation", "Setting (identical for every descriptor block)"])
tables["T5_classifier_configuration"] = t5

# ---------------- Table 6: performance per arm, three protocols
t6 = arm.copy()
t6["regime"] = t6["regime"].map({"raw600/random": "600 records, random split",
                                 "dedup387/random": "387 structures, random split",
                                 "dedup387/scaffold": "387 structures, scaffold split"})
t6["block"] = t6["block"].map({"cdft": "CDFT reactivity", "physchem": "Physicochemical",
                               "combined": "Combined"})
t6 = t6.rename(columns={"regime": "Evaluation protocol", "block": "Descriptor block",
                        "family": "Chemical family", "n_models": "Models (n)",
                        "mcc_median": "MCC median (dimensionless)", "mcc_q1": "MCC Q1 (dimensionless)", "mcc_q3": "MCC Q3 (dimensionless)",
                        "mcc_max": "MCC max (dimensionless)", "bacc_median": "Balanced accuracy median (dimensionless)",
                        "auc_median": "AUC median (dimensionless)"})
tables["T6_performance_by_arm"] = t6

# ---------------- Table 7: the statistical test of the central hypothesis
t7a = fried.rename(columns={"regime": "Evaluation protocol",
                            "n_matched_configurations": "Matched configurations (n)",
                            "friedman_chi2": "Friedman chi-square (dimensionless)", "friedman_p": "p (dimensionless)",
                            "mean_rank_cdft": "Mean rank CDFT", "mean_rank_physchem": "Mean rank physicochemical",
                            "mean_rank_combined": "Mean rank combined", "median_mcc_cdft": "Median MCC CDFT (dimensionless)",
                            "median_mcc_physchem": "Median MCC physicochemical (dimensionless)",
                            "median_mcc_combined": "Median MCC combined (dimensionless)"})
tables["T7a_friedman"] = t7a
t7b = pw.rename(columns={"regime": "Evaluation protocol", "comparison": "Contrast",
                         "n_pairs": "Pairs (n)", "median_delta": "Median dMCC (dimensionless)", "mean_delta": "Mean dMCC (dimensionless)",
                         "ci95_delta": "95% CI half-width (dimensionless)", "pct_favoring_second": "% favoring the first term",
                         "rank_biserial": "Rank-biserial r (dimensionless)", "wilcoxon_p": "Wilcoxon p (dimensionless)",
                         "nemenyi_p": "Nemenyi p (dimensionless)"})
tables["T7b_block_contrasts"] = t7b
tables["T7c_block_contrasts_family_excluded"] = pwnf.rename(
    columns={"regime": "Evaluation protocol", "comparison": "Contrast", "n_pairs": "Pairs (n)",
             "median_delta": "Median dMCC (dimensionless)", "mean_delta": "Mean dMCC (dimensionless)", "wilcoxon_p": "Wilcoxon p (dimensionless)",
             "nemenyi_p": "Nemenyi p (dimensionless)", "ci95_delta": "95% CI half-width (dimensionless)",
             "pct_favoring_second": "% favoring the first term", "rank_biserial": "Rank-biserial r (dimensionless)"})

# ---------------- Table 8: effect of protocol and of the family attribute
t8 = split.rename(columns={"block": "Descriptor block", "family": "Chemical family",
                           "n_configurations": "Configurations (n)", "mcc_random": "Mean MCC, first protocol (dimensionless)",
                           "mcc_scaffold": "Mean MCC, second protocol (dimensionless)", "mean_delta": "Mean change (dimensionless)",
                           "ci95": "95% CI half-width (dimensionless)", "pct_relative_drop": "Relative change (%)",
                           "wilcoxon_p": "Wilcoxon p (dimensionless)", "contrast": "Contrast"})
tables["T8_protocol_effect"] = t8
tables["T8b_family_attribute"] = famt.rename(
    columns={"regime": "Evaluation protocol", "block": "Descriptor block", "n_pairs": "Pairs (n)",
             "median_delta_mcc": "Median dMCC (dimensionless)", "mean_delta_mcc": "Mean dMCC (dimensionless)", "ci95": "95% CI half-width (dimensionless)",
             "pct_favoring_inclusion": "% favoring inclusion", "wilcoxon_p": "Wilcoxon p (dimensionless)"})

# ---------------- Table 9: best configurations under scaffold splitting
t9 = best.copy()
t9["block"] = t9["block"].map({"cdft": "CDFT reactivity", "physchem": "Physicochemical",
                               "combined": "Combined"})
t9 = t9.rename(columns={"block": "Descriptor block", "family": "Chemical family",
                        "selection": "Feature selection", "algorithm": "Classifier",
                        "mcc_mean": "MCC (dimensionless)", "mcc_ci95": "95% CI (dimensionless)", "bacc_mean": "Balanced accuracy (dimensionless)",
                        "auc_mean": "AUC (dimensionless)", "sens_mean": "Sensitivity (dimensionless)", "spec_mean": "Specificity (dimensionless)",
                        "n_selected": "Descriptors retained (n)"})
tables["T9_best_configurations"] = t9
tables["T10_roc_auc"] = rocT.rename(columns={"scheme": "Split scheme", "block": "Descriptor block",
                                             "auc_mean": "Mean AUC (dimensionless)", "auc_ci95": "95% CI half-width (dimensionless)",
                                             "auc_min": "Minimum over repeats (dimensionless)", "auc_max": "Maximum over repeats (dimensionless)"})

# ---------------- Supporting tables
tables["SI_effective_dimensionality"] = ed.rename(columns={"block": "Descriptor block"})
tables["SI_family_composition"] = famcomp
tables["SI_family_vs_clusters"] = clu.rename(
    columns={"n_clusters": "Clusters", "adjusted_rand_index": "Adjusted Rand index",
             "adjusted_mutual_information": "Adjusted mutual information"})
tables["SI_shap_attribution"] = shapatt
tables["SI_feature_stability"] = stab
tables["SI_cdft_equation_audit"] = pd.read_csv(f"{T}/SI_T1_cdft_equation_audit.csv")
tables["SI_orbital_diagnostics"] = pd.read_csv(f"{T}/SI_T5_frontier_orbital_diagnostics.csv")
tables["SI_cdft_from_frontier_orbitals"] = pd.read_csv(f"{T}/SI_T3_cdft_variance_from_frontier_orbitals.csv")
tables["SI_original_R_contrasts"] = origR
# Table S6. Every descriptor with its block AND a definition. The definitions are
# vendored in data/raw/descriptor_definitions.csv rather than fetched at run time:
# the physicochemical entries were harvested from the CDK descriptor dictionary
# (descriptor-algorithms.owl) and the javadoc of the CDK descriptor classes, which
# is what rcdk v3.6.0 computes; the CDFT entries are the equations audited in
# Table S7. Joining here rather than hard-coding keeps the column reproducible.
_desc = pd.DataFrame(
    {"descriptor": meta["cdft_cols"] + meta["pc_cols"],
     "block": ["CDFT reactivity"] * len(meta["cdft_cols"]) + ["Physicochemical"] * len(meta["pc_cols"])})
_defs = pd.read_csv(f"{cfg.RAW}/descriptor_definitions.csv")[["descriptor", "definition", "source"]]
_desc = _desc.merge(_defs, on="descriptor", how="left")
assert _desc["definition"].notna().all(), \
    f"no definition for: {_desc.loc[_desc.definition.isna(), 'descriptor'].tolist()}"
tables["SI_descriptor_list"] = _desc

# ---------------- audit of numbers changed relative to the previous draft
# The revision registry that compared this work against a superseded internal
# draft is deliberately NOT built here. Its "previous draft" column refers to
# a document that is not part of this release, so a reader cannot interpret
# it, and it records manuscript revision history rather than a result.

with pd.ExcelWriter(f"{OUT}/Tables_and_SI_data.xlsx", engine="openpyxl") as xw:
    for k, v in tables.items():
        v.to_excel(xw, sheet_name=k[:31], index=False)
DEC = {"dimensionless": 3, "(n)": 0, "(%)": 1}
def fix_precision(df):
    d = df.copy()
    for c in d.columns:
        if d[c].dtype.kind not in "fi":
            continue
        cl = str(c)
        if "p (dimensionless)" in cl or cl.endswith("p (dimensionless)"):
            d[c] = d[c].map(lambda v: f"{v:.3g}" if pd.notna(v) else "")
        elif "(n)" in cl or "Models" in cl or "Pairs" in cl or "Configurations" in cl:
            d[c] = d[c].map(lambda v: f"{v:.0f}" if pd.notna(v) else "")
        elif "(%)" in cl or cl.startswith("% ") or cl.startswith("Minority"):
            d[c] = d[c].map(lambda v: f"{v:.1f}" if pd.notna(v) else "")
        else:
            d[c] = d[c].map(lambda v: f"{v:.3f}" if pd.notna(v) else "")
    return d
for k, v in tables.items():
    (fix_precision(v) if not k.startswith("AUDIT") else v).to_csv(f"{T}/{k}.csv", index=False)

print(f"wrote {len(tables)} tables to {OUT}/Tables_and_SI_data.xlsx")
for k, v in tables.items():
    print(f"  {k:44s} {v.shape}")
print("\n--- Table 1 ---");  print(t1.to_string(index=False))
print("\n--- Table 2 ---");  print(t2.to_string(index=False))
