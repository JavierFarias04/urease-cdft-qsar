"""
06_interpretation.py — how much of the model's decision the reactivity block carries.

(1) SHAP attribution on the combined block, under scaffold splitting, partitioned
    by descriptor origin (physicochemical / CDFT / chemical family).
(2) Feature-selection stability: how often each descriptor survives Boruta and
    gain-based selection across strategies, family conditions and repeats.
"""
import json, os, warnings
import numpy as np
import pandas as pd
import shap
from joblib import Parallel, delayed
from sklearn.model_selection import GroupShuffleSplit, StratifiedShuffleSplit
from imblearn.over_sampling import SMOTE
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")

import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
from ureasecdft import modeling as bm
D, T = str(cfg.PROCESSED), str(cfg.TABLES)
CDFT, PC, FAMS = bm.CDFT, bm.PC, bm.FAMS
STRATEGIES = bm.STRATEGIES
SEED = bm.SEED
ded = bm.DATA["dedup387"]


def origin(col):
    if col.startswith("family_"):
        return "Chemical family"
    return "CDFT reactivity" if col in CDFT else "Physicochemical"


def build(strat_name, lo, hi, use_family, scheme, repeat):
    y_all, keep = bm.label(ded, lo, hi)
    sub = ded.loc[keep].reset_index(drop=True)
    y = y_all[keep]
    X = sub[[c for c in CDFT + PC if c in sub.columns]].astype(float).copy()
    if use_family:
        for f in FAMS:
            X[f"family_{f}"] = (sub["tipo_familia"] == f).astype(float)
    seed = SEED + repeat * 977
    if scheme == "scaffold":
        tr, te = next(GroupShuffleSplit(1, test_size=.2, random_state=seed)
                      .split(X, y, groups=sub["murcko_scaffold"].values))
    else:
        tr, te = next(StratifiedShuffleSplit(1, test_size=.2, random_state=seed).split(X, y))
    Xtr, Xte = bm.preprocess(X.iloc[tr], X.iloc[te])
    return Xtr, Xte, y[tr], y[te], seed


# --------------------------------------------------------------- (1) SHAP
def shap_cell(strat_name, lo, hi, use_family, repeat):
    Xtr, Xte, ytr, yte, seed = build(strat_name, lo, hi, use_family, "scaffold", repeat)
    if len(np.unique(ytr)) < 2 or len(np.unique(yte)) < 2:
        return None
    nmin = int(min(np.bincount(ytr)))
    A, y2 = (SMOTE(random_state=seed, k_neighbors=min(5, nmin - 1)).fit_resample(Xtr, ytr)
             if nmin >= 6 else (Xtr, ytr))
    m = XGBClassifier(n_estimators=200, max_depth=4, learning_rate=.1, subsample=.9,
                      colsample_bytree=.8, reg_lambda=1.0, eval_metric="logloss",
                      random_state=seed, n_jobs=1, verbosity=0).fit(A, y2)
    sv = shap.TreeExplainer(m).shap_values(Xte)
    mad = np.abs(sv).mean(axis=0)
    return pd.Series(mad, index=Xte.columns)


rows = []
for use_family in (True, False):
    jobs = [(sn, lo, hi, use_family, rp) for (sn, lo, hi) in STRATEGIES for rp in range(10)]
    out = Parallel(n_jobs=8, backend="threading")(delayed(shap_cell)(*j) for j in jobs)
    out = [o for o in out if o is not None]
    M = pd.DataFrame(out).fillna(0.0)
    share = M.div(M.sum(axis=1), axis=0)                      # normalize per model
    by_origin = share.T.groupby(origin).sum().T
    for k in by_origin.columns:
        rows.append(dict(family="included" if use_family else "excluded", origin=k,
                         n_models=len(M),
                         mean_share_pct=round(100 * float(by_origin[k].mean()), 2),
                         ci95_pct=round(100 * 1.96 * float(by_origin[k].std(ddof=1)) / np.sqrt(len(M)), 2)))
    tag = "family_included" if use_family else "family_excluded"
    top = (share.mean().sort_values(ascending=False).head(25).rename("mean_shap_share")
           .to_frame().assign(origin=lambda t: [origin(i) for i in t.index],
                              rank=lambda t: np.arange(1, len(t) + 1)))
    top["mean_shap_share"] = (100 * top["mean_shap_share"]).round(3)
    top.to_csv(f"{T}/T_shap_top25_{tag}.csv")
    # rank of the best CDFT descriptor
    allrank = share.mean().sort_values(ascending=False)
    ordr = [origin(i) for i in allrank.index]
    first_cdft = next((i + 1 for i, o in enumerate(ordr) if o == "CDFT reactivity"), None)
    n_cdft_top20 = sum(1 for o in ordr[:20] if o == "CDFT reactivity")
    print(f"\n== SHAP, combined block, scaffold split, family {tag.split('_')[1]} "
          f"({len(M)} models) ==")
    print(f"best-ranked CDFT descriptor sits at position {first_cdft} of {len(allrank)}; "
          f"{n_cdft_top20} CDFT descriptors in the top 20")
    print(top.head(15).to_string())
    share.mean().sort_values(ascending=False).to_csv(f"{T}/SI_T7_shap_all_{tag}.csv")

att = pd.DataFrame(rows)
att.to_csv(f"{T}/T_shap_attribution_by_origin.csv", index=False)
print("\n== share of total SHAP attribution by descriptor origin ==")
print(att.to_string(index=False))


# ----------------------------------------------- (2) feature-selection stability
def fs_cell(strat_name, lo, hi, use_family, scheme, repeat):
    Xtr, Xte, ytr, yte, seed = build(strat_name, lo, hi, use_family, scheme, repeat)
    if len(np.unique(ytr)) < 2:
        return []
    o = []
    for meth, fn in [("Boruta", bm.boruta), ("XGBoost gain", bm.xgb_gain)]:
        for f in fn(Xtr, ytr, seed):
            o.append(dict(method=meth, scheme=scheme, family="included" if use_family else "excluded",
                          strategy=strat_name, repeat=repeat, feature=f, n_candidates=Xtr.shape[1]))
    return o


jobs = [(sn, lo, hi, uf, sc, rp) for (sn, lo, hi) in STRATEGIES for uf in (True, False)
        for sc in ("random", "scaffold") for rp in range(10)]
res = Parallel(n_jobs=10, backend="threading", verbose=1)(delayed(fs_cell)(*j) for j in jobs)
sel = pd.DataFrame([x for s in res for x in s])
sel.to_csv(f"{D}/feature_selection_records.csv", index=False)

n_runs = sel.groupby(["method", "scheme", "family"])["repeat"].size()
freq = (sel.groupby(["method", "feature"]).size() /
        sel.groupby("method")[["strategy", "repeat", "scheme", "family"]]
        .apply(lambda g: len(g.drop_duplicates())))
stab = []
for meth, g in sel.groupby("method"):
    nruns = g[["scheme", "family", "strategy", "repeat"]].drop_duplicates().shape[0]
    c = g["feature"].value_counts() / nruns
    for f, v in c.items():
        stab.append(dict(method=meth, feature=f, origin=origin(f),
                         selection_frequency=round(float(v), 3), n_runs=nruns))
stab = pd.DataFrame(stab).sort_values(["method", "selection_frequency"], ascending=[True, False])
stab.to_csv(f"{T}/SI_T8_feature_selection_stability.csv", index=False)
print("\n== descriptors selected in at least half of runs ==")
for meth, g in stab.groupby("method"):
    h = g[g["selection_frequency"] >= 0.5]
    print(f"\n{meth}: {len(h)} descriptors above 0.5 "
          f"({(h['origin'] == 'CDFT reactivity').sum()} CDFT, "
          f"{(h['origin'] == 'Physicochemical').sum()} physicochemical, "
          f"{(h['origin'] == 'Chemical family').sum()} family indicators)")
    print(g.head(18)[["feature", "origin", "selection_frequency"]].to_string(index=False))

cd = stab[stab["origin"] == "CDFT reactivity"]
print("\n== retention of CDFT descriptors within the combined block ==")
print(cd.pivot_table(index="feature", columns="method", values="selection_frequency")
      .sort_values("Boruta", ascending=False).round(3).to_string())
