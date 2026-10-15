# Figure 2 — dataset structure: censoring, family/potency confound, analogue redundancy.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="open", sizes=(8, 7, 6))
from ureasecdft import config as cfg
D, T, F = str(cfg.PROCESSED), str(cfg.TABLES), str(cfg.FIGURES)
raw = pd.read_csv(f"{D}/compounds_raw.csv")
ded = pd.read_csv(f"{D}/compounds_dedup.csv")
strat = pd.read_csv(f"{T}/T_class_composition.csv")
S = np.load(f"{D}/tanimoto_similarity_dedup.npy")
audit = pd.read_csv(f"{D}/duplication_audit.csv")

from ureasecdft.style import (FAM_ORDER, FAM_LABEL, FAM_SHORT, FAM_COLOR as CB, FAM_MARKER,
                         FAM_HATCH, ACCENT, GRAY, INK, FIG_W, ZERO_STUB, panel_tag)

fig = plt.figure(figsize=(FIG_W, 5.48))
gs = fig.add_gridspec(2, 2, hspace=0.58, wspace=0.34, left=0.138, right=0.978, top=0.949, bottom=0.092)

# ---- (a) IC50 distribution by family, censored records marked -------------
ax = fig.add_subplot(gs[0, 0])
rng = np.random.default_rng(7)
for i, f in enumerate(FAM_ORDER):
    sub = raw[raw["tipo_familia"] == f]
    meas = sub.loc[~sub["censored"], "ic50_uM"].values
    cens = sub.loc[sub["censored"], "ic50_uM"].values
    ax.scatter(np.log10(meas), i + rng.uniform(-.28, .28, len(meas)), s=6,
               marker=FAM_MARKER[f], alpha=.6, color=CB[f], lw=0, zorder=3)
    ax.scatter(np.log10(cens), i + rng.uniform(-.28, .28, len(cens)), s=9, alpha=.75,
               facecolors="none", edgecolors=ACCENT, lw=.5, zorder=4)
    ax.plot([np.log10(np.median(meas))] * 2, [i - .34, i + .34], color="k", lw=1.6, zorder=6)
ax.vlines(np.log10(10), -0.45, 3.75, color=GRAY, ls="--", lw=.8, zorder=1)
ax.text(np.log10(10) + .13, -0.46, "10 µM", color=GRAY, ha="left", fontsize=6)
ax.set_yticks(range(4)); ax.set_yticklabels([FAM_SHORT[f] for f in FAM_ORDER], fontsize=6)
ax.set_xlabel("IC$_{50}$ (µM, log$_{10}$ scale)")
ax.set_xticks([-2, -1, 0, 1, 2, 3]); ax.set_xticklabels(["0.01", "0.1", "1", "10", "100", "1000"])
ax.set_ylim(-.6, 3.9); ax.margins(x=.05); ax.invert_yaxis()
# compact frameless key in the free band below the lowest family row (§7.3)
ax.legend(handles=[Line2D([], [], marker="o", ls="", ms=3.2, color=GRAY, label="measured IC$_{50}$"),
                   Line2D([], [], marker="o", ls="", ms=3.8, mfc="none", mec=ACCENT, lw=0,
                          label="censored at 1 mM"),
                   Line2D([], [], marker="|", ls="", ms=6, mew=1.6, color="k", label="family median")],
          loc="upper left", bbox_to_anchor=(-0.005, -0.275), ncol=3, fontsize=6,
          handletextpad=.3, columnspacing=.85, borderaxespad=0.0,
          handlelength=.9, frameon=False)
panel_tag(ax, "a", pad_x_pt=-40)

# ---- (b) fraction of potent compounds per family --------------------------
ax = fig.add_subplot(gs[0, 1])
fr = [100 * (raw.loc[raw["tipo_familia"] == f, "ic50_uM"] <= 10).mean() for f in FAM_ORDER]
nn = [int((raw["tipo_familia"] == f).sum()) for f in FAM_ORDER]
b = ax.barh(range(4), fr, color=[CB[f] for f in FAM_ORDER], height=.62, edgecolor="white", lw=.6)
for i, (v, n) in enumerate(zip(fr, nn)):
    ax.text(v + 1.4, i, f"{v:.1f}%" + (" (0/139)" if v == 0 else ""), va="center", fontsize=6)
    if v == 0:                                   # explicit zero stub (§D7)
        ax.barh([i], [0.95], height=.62, color=CB[FAM_ORDER[i]], edgecolor=INK, lw=.7, zorder=5)
ax.set_yticks(range(4)); ax.set_yticklabels([f"{FAM_SHORT[f]}\n(n = {n})" for f, n in zip(FAM_ORDER, nn)], fontsize=6)
ax.set_xlabel("Records with IC$_{50}$ ≤ 10 µM (%)")
ax.set_xlim(0, 68); ax.invert_yaxis()
panel_tag(ax, "b", pad_x_pt=-40)

# ---- (c) nearest-neighbor similarity: analogue redundancy ---------------
ax = fig.add_subplot(gs[1, 0])
n = S.shape[0]
nnsim = np.max(S - np.eye(n) * 2, axis=1)
iu = np.triu_indices(n, 1)
ax.hist(S[iu], bins=np.linspace(0, 1, 41), density=True, color=GRAY, alpha=.45,
        label="all compound pairs", lw=0)
ax.hist(nnsim, bins=np.linspace(0, 1, 41), density=True, color=ACCENT, alpha=.75,
        label="nearest neighbor of each compound", lw=0)
ax.axvline(0.7, color="k", ls=":", lw=.9)
ax.text(0.665, 6.2, f"{100 * (nnsim > 0.7).mean():.0f}% of structures have\na neighbor above 0.7", fontsize=6, ha="right", va="center", linespacing=1.4)
ax.set_ylim(0, 1.16 * max(np.histogram(nnsim, bins=np.linspace(0, 1, 41),
                                       density=True)[0].max(),
                          np.histogram(S[iu], bins=np.linspace(0, 1, 41),
                                       density=True)[0].max()))
ax.set_xlabel("ECFP4 Tanimoto similarity"); ax.set_ylabel("Density")
ax.legend(loc="upper center", handlelength=1.1, handletextpad=.5)
panel_tag(ax, "c", pad_x_pt=-40)

# ---- (d) structural redundancy ladder ------------------------------------
ax = fig.add_subplot(gs[1, 1])
lab = ["Records\nanalyzed", "Unique\nSMILES", "Unique\nInChIKeys", "Unique\nskeletons", "Murcko\nscaffolds"]
val = [len(raw), int(audit.loc[0, "n_unique"]), int(audit.loc[1, "n_unique"]),
       int(audit.loc[2, "n_unique"]), int(audit.loc[3, "n_unique"])]
cols = [GRAY] + ["#4a4a4a"] * 3 + [ACCENT]
ax.vlines(range(5), 0, val, color=cols, lw=1.0)
ax.scatter(range(5), val, s=34, color=cols, zorder=4)
for i, v in enumerate(val):
    ax.text(i, v + 26, str(v), ha="center", fontsize=6.5,
            fontweight="bold" if i in (0, 4) else "normal")
ax.set_xticks(range(5)); ax.set_xticklabels(lab)
ax.set_ylabel("Count"); ax.set_ylim(0, 700)
panel_tag(ax, "d", pad_x_pt=-40)

fig.savefig(f"{F}/Figure2_dataset_structure.png", dpi=400)
fig.savefig(f"{F}/Figure2_dataset_structure.pdf")

print("median nn Tanimoto %.3f | frac>0.7 %.3f" % (np.median(nnsim), (nnsim > 0.7).mean()))
