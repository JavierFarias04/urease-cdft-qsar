# Figure 5 — what the evaluation protocol and the family attribute actually buy.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="open", sizes=(8, 7, 6))
from ureasecdft import config as cfg
D, T, F = str(cfg.PROCESSED), str(cfg.TABLES), str(cfg.FIGURES)
r = pd.read_csv(f"{D}/benchmark_raw_results.csv"); r = r[r["mcc"].notna()]
roc = json.load(open(f"{D}/roc_curves.json"))
auc = pd.read_csv(f"{T}/T_roc_auc.csv")
fam = pd.read_csv(f"{T}/T_family_attribute_effect.csv")

from ureasecdft.style import (BLOCK_COLOR as BC, BLOCK_LABEL as BLAB, BLOCK_HATCH, BLOCK_DASH,
                         SPLIT_COLOR, SPLIT_MARKER, ACCENT, GRAY, INK, FIG_W,
                         ZERO_STUB, panel_tag)
from matplotlib.lines import Line2D
ALGO = ["Random forest", "XGBoost", "SVM (RBF)", "Decision tree (CART)",
        "k-nearest neighbors", "Logistic regression", "Naive Bayes"]
ASHORT = {"Random forest": "Random forest", "XGBoost": "XGBoost", "SVM (RBF)": "SVM (RBF)",
          "Decision tree (CART)": "Decision tree", "k-nearest neighbors": "kNN",
          "Logistic regression": "Logistic regr.", "Naive Bayes": "Naïve Bayes"}

fig = plt.figure(figsize=(FIG_W, 3.16))
gs = fig.add_gridspec(1, 3, width_ratios=[1.26, 1.00, 1.06], wspace=0.62,
                      left=0.082, right=0.988, top=0.911, bottom=0.224)

# ---- (a) ROC under scaffold splitting -----------------------------------
ax = fig.add_subplot(gs[0, 0])
g = np.array(roc["grid"])
handles = []
for b in ["physchem", "combined", "cdft"]:
    c = roc["curves"][f"scaffold|{b}"]
    ax.fill_between(g, c["lo"], c["hi"], color=BC[b], alpha=.15, lw=0)
    a = auc[(auc["scheme"] == "scaffold") & (auc["block"] == b)].iloc[0]
    dash = BLOCK_DASH[b][1] or (None, None)
    ax.plot(g, c["mean_tpr"], color=BC[b], lw=1.7 if b == "cdft" else 1.4, zorder=4,
            dashes=dash)
    handles.append(Line2D([0], [0], color=BC[b], lw=1.4, dashes=dash,
                          label=f"{BLAB[b]}  {a['auc_mean']:.3f} $\\pm$ {a['auc_ci95']:.3f}"))
ax.plot([0, 1], [0, 1], color=GRAY, ls=":", lw=.8)
ax.set_xlabel("False positive rate"); ax.set_ylabel("True positive rate")
ax.set_xlim(0, 1); ax.set_ylim(0, 1.02)
# The AUC values were free text at (0.30, 0.10-0.38), which sat on the lower
# edge of the reactivity block's confidence band. A framed key in the
# conventional lower-right corner of a ROC panel reads cleanly over any band.
leg = ax.legend(handles=handles, loc="lower right", fontsize=6.0,
                title="AUC (mean $\\pm$ 95% CI)",
                handlelength=2.0, handletextpad=.5, labelspacing=.32,
                borderaxespad=.25, borderpad=.45, framealpha=.92,
                facecolor="white", edgecolor="none")
leg.get_title().set_fontsize(6.0)
leg.set_zorder(6)
panel_tag(ax, "a", pad_x_pt=-32)

# ---- (b) per-algorithm random -> scaffold, physicochemical block --------
ax = fig.add_subplot(gs[0, 1])
for i, a in enumerate(ALGO):
    v1 = r[(r["regime"] == "dedup387/random") & (r["block"] == "physchem") &
           (r["algorithm"] == a)]["mcc"].mean()
    v2 = r[(r["regime"] == "dedup387/scaffold") & (r["block"] == "physchem") &
           (r["algorithm"] == a)]["mcc"].mean()
    ax.plot([v1, v2], [i, i], color=GRAY, lw=.9, zorder=2)
    ax.scatter([v1], [i], s=24, marker=SPLIT_MARKER["random"], color=SPLIT_COLOR["random"], zorder=4)
    ax.scatter([v2], [i], s=24, marker=SPLIT_MARKER["scaffold"], color=SPLIT_COLOR["scaffold"], zorder=4)
    ax.text(v2 - 0.035, i, f"−{100 * (v1 - v2) / v1:.0f}%", fontsize=6.0,
            color=SPLIT_COLOR["scaffold"], ha="right", va="center")
ax.set_yticks(range(len(ALGO))); ax.set_yticklabels([ASHORT[a] for a in ALGO], fontsize=6)
ax.invert_yaxis()
ax.set_xlabel("Mean MCC, physicochemical block")
ax.set_xlim(0.12, 0.92)
ax.text(0.825, -0.62, "random split", color=SPLIT_COLOR["random"], fontsize=6.0,
        ha="center", va="center")
ax.text(0.455, -0.62, "scaffold split", color=SPLIT_COLOR["scaffold"], fontsize=6.0,
        ha="center", va="center")
ax.set_ylim(6.6, -1.05)
panel_tag(ax, "b", pad_x_pt=-52)

# ---- (c) the chemical-family attribute ---------------------------------
ax = fig.add_subplot(gs[0, 2])
REG = ["600 records, random split", "387 structures, random split", "387 structures, scaffold split"]
RS = ["600\nrandom", "387\nrandom", "387\nscaffold"]
BK = ["CDFT reactivity", "Physicochemical", "Combined"]
KEY = {"CDFT reactivity": "cdft", "Physicochemical": "physchem", "Combined": "combined"}
w = 0.26
for j, b in enumerate(BK):
    xs, ys, es = [], [], []
    for i, rg in enumerate(REG):
        row = fam[(fam["regime"] == rg) & (fam["block"] == b)].iloc[0]
        xs.append(i + (j - 1) * w); ys.append(row["mean_delta_mcc"]); es.append(row["ci95"])
    ax.bar(xs, ys, width=w * 0.92, color=BC[KEY[b]], alpha=.80, yerr=es, capsize=1.6,
           error_kw=dict(lw=.7), label=b, hatch=BLOCK_HATCH[KEY[b]], edgecolor="white", lw=.5)
    # A bar of 0.000 draws nothing and reads as missing data, so every value
    # too small to render gets an explicit hairline stub at the baseline.
    for xx, yy in zip(xs, ys):
        if abs(yy) < ZERO_STUB:
            ax.bar([xx], [ZERO_STUB], width=w * 0.92, color=INK, lw=0, zorder=4)
ax.axhline(0, color="k", lw=.9)
ax.set_xticks(range(3)); ax.set_xticklabels(RS, fontsize=6)
ax.set_xlabel("Compounds / split scheme", fontsize=6.5)
ax.set_ylabel("$\\Delta$MCC from adding\nthe chemical family")
ax.set_ylim(-0.045, 0.30)
ax.legend(loc="upper right", fontsize=6.0, handlelength=1.0, handletextpad=.4,
          labelspacing=.25, borderaxespad=.2)
panel_tag(ax, "c", pad_x_pt=-40)

fig.savefig(f"{F}/Figure6_protocol_effects.png", dpi=400)
fig.savefig(f"{F}/Figure6_protocol_effects.pdf")

