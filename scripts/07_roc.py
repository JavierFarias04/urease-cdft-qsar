"""
07_roc.py — pooled ROC curves for the three descriptor blocks under both split
schemes, using the configuration that performed best under scaffold splitting
(XGBoost with gain-based selection, 10 uM single cutoff).
"""
import json, warnings
import numpy as np
import pandas as pd
from sklearn.metrics import roc_curve, roc_auc_score
from sklearn.model_selection import GroupShuffleSplit, StratifiedShuffleSplit
from imblearn.over_sampling import SMOTE
from xgboost import XGBClassifier

warnings.filterwarnings("ignore")
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
from ureasecdft import modeling as bm
D = str(cfg.PROCESSED)
ded = bm.DATA["dedup387"]
GRID = np.linspace(0, 1, 201)
rows, curves = [], {}

for scheme in ("random", "scaffold"):
    for block, cols in bm.BLOCKS.items():
        tprs, aucs = [], []
        for rp in range(10):
            y_all, keep = bm.label(ded, 10, None)
            sub = ded.loc[keep].reset_index(drop=True)
            y = y_all[keep]
            X = sub[[c for c in cols if c in sub.columns]].astype(float)
            seed = bm.SEED + rp * 977
            if scheme == "scaffold":
                tr, te = next(GroupShuffleSplit(1, test_size=.2, random_state=seed)
                              .split(X, y, groups=sub["murcko_scaffold"].values))
            else:
                tr, te = next(StratifiedShuffleSplit(1, test_size=.2, random_state=seed).split(X, y))
            Xtr, Xte = bm.preprocess(X.iloc[tr], X.iloc[te])
            ytr, yte = y[tr], y[te]
            sel = bm.xgb_gain(Xtr, ytr, seed)
            nmin = int(min(np.bincount(ytr)))
            A, y2 = (SMOTE(random_state=seed, k_neighbors=min(5, nmin - 1))
                     .fit_resample(Xtr[sel], ytr) if nmin >= 6 else (Xtr[sel], ytr))
            m = XGBClassifier(n_estimators=200, max_depth=4, learning_rate=.1, subsample=.9,
                              colsample_bytree=.8, reg_lambda=1.0, eval_metric="logloss",
                              random_state=seed, n_jobs=2, verbosity=0).fit(A, y2)
            p = m.predict_proba(Xte[sel])[:, 1]
            fpr, tpr, _ = roc_curve(yte, p)
            tprs.append(np.interp(GRID, fpr, tpr))
            aucs.append(roc_auc_score(yte, p))
        tprs = np.vstack(tprs)
        curves[f"{scheme}|{block}"] = dict(mean_tpr=tprs.mean(0).tolist(),
                                           lo=np.percentile(tprs, 5, axis=0).tolist(),
                                           hi=np.percentile(tprs, 95, axis=0).tolist())
        rows.append(dict(scheme=scheme, block=block, auc_mean=round(float(np.mean(aucs)), 4),
                         auc_ci95=round(float(1.96 * np.std(aucs, ddof=1) / np.sqrt(len(aucs))), 4),
                         auc_min=round(float(np.min(aucs)), 4), auc_max=round(float(np.max(aucs)), 4)))

json.dump(dict(grid=GRID.tolist(), curves=curves), open(f"{D}/roc_curves.json", "w"))
t = pd.DataFrame(rows)
t.to_csv(f"{cfg.TABLES}/T_roc_auc.csv", index=False)
print(t.to_string(index=False))
