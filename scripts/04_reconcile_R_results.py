"""
04_reconcile_R_results.py — tidy the original R benchmark and apply the paired
statistics the manuscript needs.

The six workbooks hold 126 models each (7 algorithms x 6 categorization strategies
x 3 feature-selection regimes) for one (descriptor block, chemical-family) arm,
i.e. 756 models in total. Because the three blocks were run on different compound
sets (Set 1 n=677 from the 2024 study, Set 2 n=649, Set 3 n=600), the comparison
between blocks in these tables is not controlled; the purpose here is to quantify
what that uncontrolled comparison implies and to contrast it with the controlled
re-run in 03_benchmark.py.
"""
import os, itertools
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon
import scikit_posthocs as sp

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
BASE = str(cfg.R_RESULTS_DIR)
D, T = str(cfg.PROCESSED), str(cfg.TABLES)

FILES = {
    ("cdft", "included"): "resultados_ConFQ_solo_koopmans.xlsx",
    ("cdft", "excluded"): "resultados_SinFQ_solo_koopmans.xlsx",
    ("combined", "included"): "resultados_ConFQ_koopmans_mas_datos_ureasa.xlsx",
    ("combined", "excluded"): "resultados_SinFQ_koopmans_mas_datos_ureasa.xlsx",
    ("physchem", "included"): "Resultados_ureasa_ConFQ.xlsx",
    ("physchem", "excluded"): "Resultados_ureasa_SinFQ.xlsx",
}
ALGO = {"C50Tree": "Decision tree (C5.0)", "KNN": "k-nearest neighbors", "NB": "Naive Bayes",
        "RF": "Random forest", "RL": "Logistic regression", "Rlogistica": "Logistic regression",
        "SVM": "SVM (RBF)", "XGB": "XGBoost"}
SEL = {"NFS": "none", "SSA": "none", "BORUTA": "Boruta", "XGB": "XGBoost gain"}
NCOMP = {"physchem": 677, "cdft": 649, "combined": 600}

frames = []
for (blk, fam), fn in FILES.items():
    x = pd.ExcelFile(f"{BASE}/{fn}")
    df = x.parse(x.sheet_names[0])
    df.columns = [str(c).strip() for c in df.columns]
    df = df.rename(columns={"Algoritmo": "algorithm_raw", "Corte": "strategy_raw", "SA": "selection_raw",
                            "Balanced_Accuracy": "balanced_accuracy", "MCC": "mcc", "F1": "f1",
                            "Accuracy": "accuracy", "Sensitivity": "sensitivity",
                            "Specificity": "specificity", "Precision": "precision"})
    df["block"], df["family"], df["source_file"] = blk, fam, fn
    df["n_compounds"] = NCOMP[blk]
    frames.append(df)

r = pd.concat(frames, ignore_index=True)
r["algorithm"] = r["algorithm_raw"].astype(str).str.strip().map(ALGO)
r["selection"] = r["selection_raw"].astype(str).str.strip().str.upper().map(SEL)
# strategy strings appear as "5", "5 y 50", or "5 SinFQ" in different workbooks
s = r["strategy_raw"].astype(str).str.strip()
s = s.str.replace(r"\s*(SinFQ|ConFQ)\s*", "", regex=True).str.replace(" y ", " and ", regex=False)
r["strategy"] = s.where(s.str.contains("and"), s) + " uM"
r["strategy"] = r["strategy"].str.replace("and 50 uM", "and 50 uM", regex=False)
r["strategy"] = r["strategy"].str.replace(r"^(\d+) and 50 uM$", r"\1 and 50 uM", regex=True)

assert r["algorithm"].notna().all(), r.loc[r['algorithm'].isna(), 'algorithm_raw'].unique()
assert r["selection"].notna().all(), r.loc[r['selection'].isna(), 'selection_raw'].unique()
print(f"tidied {len(r)} original R models "
      f"({r.groupby(['block','family']).size().to_dict()})")
print("strategies found:", sorted(r['strategy'].unique()))
r.to_csv(f"{D}/original_R_results_tidy.csv", index=False)

# --------------------------------------------------- paired comparison of blocks
def paired_block_test(sub, metric="mcc"):
    piv = sub.pivot_table(index=["family", "strategy", "selection", "algorithm"],
                          columns="block", values=metric).dropna()
    cols = ["cdft", "physchem", "combined"]
    piv = piv[cols]
    st, p = friedmanchisquare(*[piv[c].values for c in cols])
    nem = sp.posthoc_nemenyi_friedman(piv.values)
    nem.index = nem.columns = cols
    ranks = piv.rank(axis=1, ascending=False).mean()
    pw = []
    for a, b in itertools.combinations(cols, 2):
        d = piv[b] - piv[a]
        w = wilcoxon(piv[b], piv[a])
        pw.append(dict(comparison=f"{b} - {a}", median_delta=round(float(d.median()), 4),
                       mean_delta=round(float(d.mean()), 4),
                       q1=round(float(d.quantile(.25)), 4), q3=round(float(d.quantile(.75)), 4),
                       n_pairs=len(d), n_favoring_second=int((d > 0).sum()),
                       wilcoxon_p=float(w.pvalue),
                       nemenyi_p=float(nem.loc[a, b])))
    return piv, st, p, ranks, pd.DataFrame(pw)

piv, st, p, ranks, pw = paired_block_test(r)
print(f"\n== ORIGINAL R BENCHMARK, paired across {len(piv)} matched "
      f"(family, strategy, selection, algorithm) configurations ==")
print(f"Friedman chi2 = {st:.2f}, p = {p:.3e}")
print("mean rank (1 = best):", ranks.round(2).to_dict())
print(pw.to_string(index=False))
print("\nmedian MCC by block:", piv.median().round(3).to_dict())
pw.to_csv(f"{T}/T_original_R_block_comparison.csv", index=False)
piv.to_csv(f"{D}/original_R_paired_mcc.csv")

# family effect
fam = r.pivot_table(index=["block", "strategy", "selection", "algorithm"],
                    columns="family", values="mcc").dropna()
d = fam["included"] - fam["excluded"]
w = wilcoxon(fam["included"], fam["excluded"])
print(f"\n== chemical-family attribute, original R benchmark ({len(d)} matched pairs) ==")
print(f"median dMCC = {d.median():.4f} (IQR {d.quantile(.25):.4f} to {d.quantile(.75):.4f}), "
      f"favoring inclusion in {int((d>0).sum())}/{len(d)}, Wilcoxon p = {w.pvalue:.3e}")
byblk = r.pivot_table(index=["strategy", "selection", "algorithm"], columns=["block", "family"], values="mcc")
rows = []
for blk in ["cdft", "physchem", "combined"]:
    sub = byblk[blk].dropna()
    dd = sub["included"] - sub["excluded"]
    rows.append(dict(block=blk, n_pairs=len(dd), median_delta_mcc=round(float(dd.median()), 4),
                     wilcoxon_p=float(wilcoxon(sub["included"], sub["excluded"]).pvalue)))
fam_tab = pd.DataFrame(rows)
fam_tab.to_csv(f"{T}/T_original_R_family_effect.csv", index=False)
print(fam_tab.to_string(index=False))

# maxima vs medians: the inflation the draft's Table 7 suffers from
mx = r.groupby(["block", "family"])["mcc"].agg(n="size", median="median", q90=lambda s: s.quantile(.9), max="max").round(3)
mx.to_csv(f"{T}/T_original_R_maxima_vs_median.csv")
print("\n== distribution of MCC per arm (shows how far the reported maxima sit above the typical model) ==")
print(mx.to_string())
