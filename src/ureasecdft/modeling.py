"""Shared modeling logic: labeling, preprocessing, feature selection, classifiers.

Used by scripts 03, 06 and 07. The function bodies are the ones that produced
the published results, moved here verbatim so the three scripts cannot drift
apart. Previously 06 and 07 loaded 03 through importlib from a hardcoded
relative path, which coupled them to the working layout of the original
analysis directory.
"""
import json

import numpy as np
import pandas as pd
from scipy.stats import binomtest
from sklearn.ensemble import RandomForestClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.naive_bayes import GaussianNB
from sklearn.neighbors import KNeighborsClassifier
from sklearn.preprocessing import StandardScaler
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from xgboost import XGBClassifier

from . import config as cfg

SEED = cfg.SEED

_meta = json.load(open(cfg.PROCESSED / "build_meta.json"))
PC = list(_meta["pc_cols"])
CDFT = list(_meta["cdft_cols"])
FAMS = ["Coumarin", "Tiazol", "Triazol", "Otros"]
STRATEGIES = [("5 uM", 5, None), ("10 uM", 10, None), ("25 uM", 25, None),
              ("5 and 50 uM", 5, 50), ("10 and 50 uM", 10, 50), ("25 and 50 uM", 25, 50)]
REGIMES = ["raw600/random", "dedup387/random", "dedup387/scaffold"]
BLOCKS = {"cdft": CDFT, "physchem": PC, "combined": CDFT + PC}
DATA = {"raw600": pd.read_csv(cfg.PROCESSED / "compounds_raw.csv"),
        "dedup387": pd.read_csv(cfg.PROCESSED / "compounds_dedup.csv")}


def classifiers(seed):
    return {
        "Random forest": RandomForestClassifier(n_estimators=500, random_state=seed, n_jobs=1),
        "SVM (RBF)": SVC(C=1.0, gamma="scale", probability=True, random_state=seed),
        "Decision tree (CART)": DecisionTreeClassifier(random_state=seed),
        "XGBoost": XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.1,
                                 subsample=0.9, colsample_bytree=0.8, reg_lambda=1.0,
                                 eval_metric="logloss", random_state=seed, n_jobs=1, verbosity=0),
        "k-nearest neighbors": KNeighborsClassifier(n_neighbors=5),
        "Naive Bayes": GaussianNB(),
        "Logistic regression": LogisticRegression(max_iter=3000, C=1.0, random_state=seed),
    }


# ------------------------------------------------------------------ labeling
def label(df, lo, hi):
    v = df["ic50_uM"].astype(float).values
    y = (v <= lo).astype(int)
    keep = np.ones(len(v), bool) if hi is None else ((v <= lo) | (v > hi))
    return y, keep


# ------------------------------------------------------------- preprocessing
def preprocess(Xtr, Xte):
    """Fit every step on the training partition only."""
    keep = Xtr.isna().mean() < 0.8
    Xtr, Xte = Xtr.loc[:, keep], Xte.loc[:, keep]
    med = Xtr.median()
    Xtr, Xte = Xtr.fillna(med), Xte.fillna(med)
    # near-zero variance: caret's default freqCut = 95/5, uniqueCut = 10
    ok = []
    n = len(Xtr)
    for c in Xtr.columns:
        vc = Xtr[c].value_counts()
        if len(vc) < 2:
            continue
        if (vc.iloc[0] / max(vc.iloc[1], 1) > 19) and (len(vc) / n * 100 < 10):
            continue
        ok.append(c)
    Xtr, Xte = Xtr[ok], Xte[ok]
    sc = StandardScaler().fit(Xtr)
    return (pd.DataFrame(sc.transform(Xtr), columns=ok, index=Xtr.index),
            pd.DataFrame(sc.transform(Xte), columns=ok, index=Xte.index))


# ---------------------------------------------------------- feature selection
def boruta(X, y, seed, max_iter=20, alpha=0.05):
    """Kursa & Rudnicki shadow-feature wrapper around a random forest."""
    rng = np.random.default_rng(seed)
    cols = list(X.columns)
    hits = np.zeros(len(cols))
    Xv = X.values
    for it in range(max_iter):
        sh = np.apply_along_axis(rng.permutation, 0, Xv)
        Z = np.hstack([Xv, sh])
        rf = RandomForestClassifier(n_estimators=200, random_state=seed + it, n_jobs=1).fit(Z, y)
        imp = rf.feature_importances_
        hits += (imp[:len(cols)] > imp[len(cols):].max()).astype(int)
    p = np.array([binomtest(int(h), max_iter, 0.5, alternative="greater").pvalue for h in hits])
    sel = [c for c, pv in zip(cols, p * len(cols)) if pv < alpha]      # Bonferroni
    return sel or [cols[int(np.argmax(hits))]]


def xgb_gain(X, y, seed, cum=0.95):
    m = XGBClassifier(n_estimators=200, max_depth=4, learning_rate=0.1, eval_metric="logloss",
                      random_state=seed, n_jobs=1, verbosity=0).fit(X, y)
    g = pd.Series(m.get_booster().get_score(importance_type="gain"))
    if g.empty:
        return list(X.columns)
    g = g.sort_values(ascending=False)
    keep = g.index[(g.cumsum() / g.sum()) <= cum].tolist()
    return (keep or [g.index[0]])


