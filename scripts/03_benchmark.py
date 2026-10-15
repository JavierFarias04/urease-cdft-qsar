"""
03_benchmark.py — controlled factorial benchmark of the three descriptor blocks.

Design
------
All three descriptor blocks are evaluated on IDENTICAL compounds, with identical
preprocessing, identical hyperparameters and identical partitions, so that the
only quantity varying between arms is descriptor content.

  evaluation regimes (3)
      raw600/random     600 records, stratified random 80/20  (replicates the original protocol)
      dedup387/random   387 unique structures, stratified random 80/20
      dedup387/scaffold 387 unique structures, Bemis-Murcko grouped 80/20 (no scaffold spans partitions)
  descriptor blocks (3)  cdft | physchem | combined
  chemical family (2)    included (one-hot) | excluded
  categorization (6)     single cutoff 5/10/25 uM; two cutoff 5,10,25 vs 50 uM (gray zone dropped)
  feature selection (3)  none | Boruta | XGBoost gain (95% cumulative gain)
  classifiers (7)        RF, SVM-RBF, decision tree (CART), XGBoost, kNN, Gaussian NB, logistic regression
  repeats (10)           paired across every arm by repeat index

Hyperparameters are FIXED and identical for every descriptor block: unequal tuning
effort between arms would confound the block comparison, which is the question here.
A separately tuned run is reported in 05_tuned_best.py.

Preprocessing is fitted on the training partition only: drop >=80% missing, median
impute, drop near-zero variance, standardize, then SMOTE on the training fold alone.
"""
import json, os, sys, warnings, time
import numpy as np
import pandas as pd
from joblib import Parallel, delayed
from sklearn.ensemble import RandomForestClassifier
from sklearn.svm import SVC
from sklearn.tree import DecisionTreeClassifier
from sklearn.neighbors import KNeighborsClassifier
from sklearn.naive_bayes import GaussianNB
from sklearn.linear_model import LogisticRegression
from sklearn.preprocessing import StandardScaler
from sklearn.model_selection import StratifiedShuffleSplit, GroupShuffleSplit
from sklearn.metrics import (accuracy_score, balanced_accuracy_score, matthews_corrcoef,
                             f1_score, roc_auc_score, confusion_matrix)
from imblearn.over_sampling import SMOTE
from xgboost import XGBClassifier
from scipy.stats import binomtest

warnings.filterwarnings("ignore")
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
D, T = str(cfg.PROCESSED), str(cfg.TABLES)
N_REPEATS = int(os.environ.get("N_REPEATS", 10))
N_JOBS = int(os.environ.get("N_JOBS", 10))

# Labeling, preprocessing, feature selection, the classifier set and the
# descriptor blocks all live in the shared module, so steps 06 and 07 use
# exactly the same code as this one.
from ureasecdft.modeling import (SEED, PC, CDFT, FAMS, STRATEGIES, REGIMES,
                                  BLOCKS, DATA, classifiers, label, preprocess,
                                  boruta, xgb_gain)


# --------------------------------------------------------------- one workload
def run_cell(regime, block, use_family, strat_name, lo, hi, repeat):
    tag = dict(regime=regime, block=block, family="included" if use_family else "excluded",
               strategy=strat_name, repeat=repeat)
    dataset, scheme = regime.split("/")
    df = DATA[dataset]
    y_all, keep = label(df, lo, hi)
    sub = df.loc[keep].reset_index(drop=True)
    y = y_all[keep]
    if len(np.unique(y)) < 2 or min(np.bincount(y)) < 12:
        return []

    cols = [c for c in BLOCKS[block] if c in sub.columns]
    X = sub[cols].astype(float).copy()
    if use_family:
        for f in FAMS:
            X[f"family_{f}"] = (sub["tipo_familia"] == f).astype(float)

    seed = SEED + repeat * 977
    if scheme == "random":
        spl = StratifiedShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
        tr, te = next(spl.split(X, y))
    else:
        spl = GroupShuffleSplit(n_splits=1, test_size=0.2, random_state=seed)
        tr, te = next(spl.split(X, y, groups=sub["murcko_scaffold"].values))
    if len(np.unique(y[te])) < 2 or len(np.unique(y[tr])) < 2:
        return []

    Xtr, Xte = preprocess(X.iloc[tr], X.iloc[te])
    ytr, yte = y[tr], y[te]

    out = []
    for fs in ("none", "Boruta", "XGBoost gain"):
        if fs == "none":
            sel = list(Xtr.columns)
        elif fs == "Boruta":
            sel = boruta(Xtr, ytr, seed)
        else:
            sel = xgb_gain(Xtr, ytr, seed)
        A, B = Xtr[sel], Xte[sel]
        nmin = int(min(np.bincount(ytr)))
        if nmin >= 6:
            A2, y2 = SMOTE(random_state=seed, k_neighbors=min(5, nmin - 1)).fit_resample(A, ytr)
        else:
            A2, y2 = A, ytr
        for name, clf in classifiers(seed).items():
            try:
                clf.fit(A2, y2)
                yp = clf.predict(B)
                try:
                    sc = clf.predict_proba(B)[:, 1]
                except Exception:
                    sc = clf.decision_function(B)
                tn, fp, fn, tp = confusion_matrix(yte, yp, labels=[0, 1]).ravel()
                out.append(dict(**tag, selection=fs, n_selected=len(sel), algorithm=name,
                                n_train=len(y2), n_test=len(yte),
                                accuracy=accuracy_score(yte, yp),
                                balanced_accuracy=balanced_accuracy_score(yte, yp),
                                mcc=matthews_corrcoef(yte, yp),
                                f1=f1_score(yte, yp, zero_division=0),
                                sensitivity=tp / (tp + fn) if tp + fn else np.nan,
                                specificity=tn / (tn + fp) if tn + fp else np.nan,
                                auc=roc_auc_score(yte, sc)))
            except Exception as e:
                out.append(dict(**tag, selection=fs, n_selected=len(sel), algorithm=name, error=str(e)[:120]))
    return out


if __name__ == "__main__":
    jobs = [(rg, bl, fam, sn, lo, hi, rp)
            for rg in REGIMES for bl in BLOCKS for fam in (True, False)
            for (sn, lo, hi) in STRATEGIES for rp in range(N_REPEATS)]
    print(f"{len(jobs)} workloads x 3 selection regimes x 7 classifiers "
          f"= {len(jobs) * 21} model fits", flush=True)
    t0 = time.time()
    res = Parallel(n_jobs=N_JOBS, verbose=5, backend="threading")(delayed(run_cell)(*j) for j in jobs)
    rows = [r for sub in res for r in sub]
    out = pd.DataFrame(rows)
    out.to_csv(f"{D}/benchmark_raw_results.csv", index=False)
    print(f"\n{len(out)} model results in {(time.time() - t0) / 60:.1f} min "
          f"({out['mcc'].notna().sum()} with metrics, {out.get('error', pd.Series(dtype=object)).notna().sum()} errors)")
    print(out.groupby(["regime", "block"])["mcc"].agg(["count", "mean", "max"]).round(3).to_string())
