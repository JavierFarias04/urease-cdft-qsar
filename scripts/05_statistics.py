"""
05_statistics.py — paired statistics on the controlled benchmark.

Each (regime, family, strategy, selection, algorithm, repeat) cell is matched
across the three descriptor blocks, so the block comparison is paired on every
other factor. Tests: Friedman across blocks with post hoc Nemenyi; paired
Wilcoxon for the pairwise increments; Wilcoxon for the chemical-family attribute
and for the random-vs-scaffold contrast.
"""
import itertools, json, os
import numpy as np
import pandas as pd
from scipy.stats import friedmanchisquare, wilcoxon
import scikit_posthocs as sp

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
D, T = str(cfg.PROCESSED), str(cfg.TABLES)
BLOCKS = ["cdft", "physchem", "combined"]
BLABEL = {"cdft": "CDFT reactivity", "physchem": "Physicochemical", "combined": "Combined"}
RLABEL = {"raw600/random": "600 records, random split",
          "dedup387/random": "387 structures, random split",
          "dedup387/scaffold": "387 structures, scaffold split"}
KEY = ["family", "strategy", "selection", "algorithm", "repeat"]

r = pd.read_csv(f"{D}/benchmark_raw_results.csv")
r = r[r["mcc"].notna()].copy()
print(f"{len(r)} model results with metrics; "
      f"{r['error'].notna().sum() if 'error' in r else 0} failures")
print(r.groupby(["regime", "block"])["mcc"].agg(["count", "median", "max"]).round(3).to_string())


def ci95(x):
    x = np.asarray(x, float)
    return 1.96 * x.std(ddof=1) / np.sqrt(len(x))


# ---------------------------------------------- 1. headline metrics per arm
summ = (r.groupby(["regime", "block", "family", "selection", "algorithm"])
        .agg(n_repeats=("mcc", "size"), mcc_mean=("mcc", "mean"), mcc_ci95=("mcc", ci95),
             bacc_mean=("balanced_accuracy", "mean"), bacc_ci95=("balanced_accuracy", ci95),
             auc_mean=("auc", "mean"), auc_ci95=("auc", ci95),
             acc_mean=("accuracy", "mean"), sens_mean=("sensitivity", "mean"),
             spec_mean=("specificity", "mean"), n_selected=("n_selected", "mean"))
        .reset_index().round(4))
summ.to_csv(f"{T}/SI_T6_benchmark_summary_by_configuration.csv", index=False)

arm = (r.groupby(["regime", "block", "family"])
       .agg(n_models=("mcc", "size"), mcc_median=("mcc", "median"),
            mcc_q1=("mcc", lambda s: s.quantile(.25)), mcc_q3=("mcc", lambda s: s.quantile(.75)),
            mcc_max=("mcc", "max"), bacc_median=("balanced_accuracy", "median"),
            auc_median=("auc", "median"))
       .reset_index())
for c in ["mcc_median", "mcc_q1", "mcc_q3", "mcc_max", "bacc_median", "auc_median"]:
    arm[c] = arm[c].map(lambda v: f"{v:.3f}")
arm.to_csv(f"{T}/T_benchmark_by_arm.csv", index=False)
print("\n== MCC distribution per arm ==")
print(arm.to_string(index=False))


# ------------------------------------------- 2. block comparison, per regime
def block_stats(sub, metric="mcc"):
    piv = sub.pivot_table(index=KEY, columns="block", values=metric).dropna()[BLOCKS]
    st, p = friedmanchisquare(*[piv[c].values for c in BLOCKS])
    nem = sp.posthoc_nemenyi_friedman(piv.values)
    nem.index = nem.columns = BLOCKS
    ranks = piv.rank(axis=1, ascending=False).mean()
    rows = []
    for a, b in itertools.combinations(BLOCKS, 2):
        d = piv[b] - piv[a]
        w = wilcoxon(piv[b], piv[a])
        # rank-biserial effect size
        rb = float((d > 0).mean() - (d < 0).mean())
        rows.append(dict(comparison=f"{BLABEL[b]} - {BLABEL[a]}", n_pairs=len(d),
                         median_delta=round(float(d.median()), 4), mean_delta=round(float(d.mean()), 4),
                         ci95_delta=round(float(ci95(d)), 4),
                         pct_favoring_second=round(100 * float((d > 0).mean()), 1),
                         rank_biserial=round(rb, 3), wilcoxon_p=float(w.pvalue),
                         nemenyi_p=float(nem.loc[a, b])))
    return piv, st, p, ranks, pd.DataFrame(rows)


allpw, fried = [], []
for rg, sub in r.groupby("regime"):
    piv, st, p, ranks, pw = block_stats(sub)
    pw.insert(0, "regime", RLABEL[rg])
    allpw.append(pw)
    fried.append(dict(regime=RLABEL[rg], n_matched_configurations=len(piv),
                      friedman_chi2=round(float(st), 2), friedman_p=float(p),
                      **{f"mean_rank_{k}": round(float(v), 2) for k, v in ranks.items()},
                      **{f"median_mcc_{k}": round(float(piv[k].median()), 3) for k in BLOCKS}))
    piv.to_csv(f"{D}/paired_mcc_{rg.replace('/', '_')}.csv")
fried = pd.DataFrame(fried)
allpw = pd.concat(allpw, ignore_index=True)
fried.to_csv(f"{T}/T_friedman_by_regime.csv", index=False)
allpw.to_csv(f"{T}/T_block_pairwise_by_regime.csv", index=False)
print("\n== Friedman across descriptor blocks, by evaluation regime ==")
print(fried.to_string(index=False))
print("\n== pairwise block increments (paired on family x strategy x selection x algorithm x repeat) ==")
print(allpw.to_string(index=False))

# same test restricted to the configurations WITHOUT the family attribute,
# i.e. where the descriptors must carry the signal themselves
print("\n== block increments with the chemical-family attribute EXCLUDED ==")
sub_nf = []
for rg, sub in r[r["family"] == "excluded"].groupby("regime"):
    _, st, p, ranks, pw = block_stats(sub)
    pw.insert(0, "regime", RLABEL[rg])
    sub_nf.append(pw)
sub_nf = pd.concat(sub_nf, ignore_index=True)
sub_nf.to_csv(f"{T}/T_block_pairwise_family_excluded.csv", index=False)
print(sub_nf.to_string(index=False))


# --------------------------------------------- 3. chemical-family attribute
rows = []
for (rg, blk), sub in r.groupby(["regime", "block"]):
    piv = sub.pivot_table(index=["strategy", "selection", "algorithm", "repeat"],
                          columns="family", values="mcc").dropna()
    d = piv["included"] - piv["excluded"]
    rows.append(dict(regime=RLABEL[rg], block=BLABEL[blk], n_pairs=len(d),
                     median_delta_mcc=round(float(d.median()), 4),
                     mean_delta_mcc=round(float(d.mean()), 4), ci95=round(float(ci95(d)), 4),
                     pct_favoring_inclusion=round(100 * float((d > 0).mean()), 1),
                     wilcoxon_p=float(wilcoxon(piv["included"], piv["excluded"]).pvalue)))
fam = pd.DataFrame(rows)
fam.to_csv(f"{T}/T_family_attribute_effect.csv", index=False)
print("\n== effect of including the chemical family as an explicit attribute ==")
print(fam.to_string(index=False))


# ------------------------------------- 4. random vs scaffold split (leakage)
rows = []
for blk in BLOCKS:
    for famc in ["included", "excluded"]:
        a = r[(r["regime"] == "dedup387/random") & (r["block"] == blk) & (r["family"] == famc)]
        b = r[(r["regime"] == "dedup387/scaffold") & (r["block"] == blk) & (r["family"] == famc)]
        ka = ["strategy", "selection", "algorithm"]
        A = a.groupby(ka)["mcc"].mean()
        B = b.groupby(ka)["mcc"].mean()
        j = pd.concat([A.rename("random"), B.rename("scaffold")], axis=1).dropna()
        d = j["scaffold"] - j["random"]
        rows.append(dict(block=BLABEL[blk], family=famc, n_configurations=len(j),
                         mcc_random=round(float(j["random"].mean()), 3),
                         mcc_scaffold=round(float(j["scaffold"].mean()), 3),
                         mean_delta=round(float(d.mean()), 3), ci95=round(float(ci95(d)), 3),
                         pct_relative_drop=round(-100 * float(d.mean() / j["random"].mean()), 1),
                         wilcoxon_p=float(wilcoxon(j["scaffold"], j["random"]).pvalue)))
# also the raw600 -> dedup387 contrast, isolating the duplication effect
for blk in BLOCKS:
    for famc in ["included", "excluded"]:
        a = r[(r["regime"] == "raw600/random") & (r["block"] == blk) & (r["family"] == famc)]
        b = r[(r["regime"] == "dedup387/random") & (r["block"] == blk) & (r["family"] == famc)]
        ka = ["strategy", "selection", "algorithm"]
        j = pd.concat([a.groupby(ka)["mcc"].mean().rename("raw600"),
                       b.groupby(ka)["mcc"].mean().rename("dedup387")], axis=1).dropna()
        d = j["dedup387"] - j["raw600"]
        rows.append(dict(block=BLABEL[blk], family=famc, n_configurations=len(j),
                         mcc_random=round(float(j["raw600"].mean()), 3),
                         mcc_scaffold=round(float(j["dedup387"].mean()), 3),
                         mean_delta=round(float(d.mean()), 3), ci95=round(float(ci95(d)), 3),
                         pct_relative_drop=round(-100 * float(d.mean() / j["raw600"].mean()), 1),
                         wilcoxon_p=float(wilcoxon(j["dedup387"], j["raw600"]).pvalue),
                         contrast="duplicate removal (600 -> 387, both random)"))
leak = pd.DataFrame(rows)
leak["contrast"] = leak["contrast"].fillna("scaffold vs random split (387 structures)")
leak.to_csv(f"{T}/T_split_scheme_effect.csv", index=False)
print("\n== effect of the evaluation protocol on apparent performance ==")
print(leak.to_string(index=False))


# ------------------------------------------------ 5. top configurations table
best = (summ[summ["regime"] == "dedup387/scaffold"]
        .sort_values("mcc_mean", ascending=False).head(15)
        [["block", "family", "selection", "algorithm", "mcc_mean", "mcc_ci95",
          "bacc_mean", "auc_mean", "sens_mean", "spec_mean", "n_selected"]])
best.to_csv(f"{T}/T_top_configurations_scaffold.csv", index=False)
print("\n== 15 best configurations under scaffold splitting (mean over 10 repeats x 6 strategies) ==")
print(best.to_string(index=False))

json.dump(dict(n_results=int(len(r))), open(f"{D}/stats_meta.json", "w"))
