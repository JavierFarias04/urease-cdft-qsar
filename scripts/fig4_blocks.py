# Figure 4 — the controlled block comparison: the reactivity block adds nothing.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="open", sizes=(8, 7, 6))
from ureasecdft import config as cfg
D, T, F = str(cfg.PROCESSED), str(cfg.TABLES), str(cfg.FIGURES)
r = pd.read_csv(f"{D}/benchmark_raw_results.csv")
r = r[r["mcc"].notna()]
pw = pd.read_csv(f"{T}/T_block_pairwise_by_regime.csv")

REG = ["raw600/random", "dedup387/random", "dedup387/scaffold"]
RLAB = {"raw600/random": "600 records\nrandom", "dedup387/random": "387 structures\nrandom",
        "dedup387/scaffold": "387 structures\nscaffold"}
SHORT = {"raw600/random": "600\nrandom", "dedup387/random": "387\nrandom",
         "dedup387/scaffold": "387\nscaffold"}
from ureasecdft.style import (BLOCK_COLOR as BC, BLOCK_LABEL as BLAB, BLOCK_HATCH, BLOCK_MARKER,
                         BLOCK_DASH, ACCENT, GRAY, INK, FIG_W, panel_tag)
BLK = ["cdft", "physchem", "combined"]

# right=0.875 reserves the strip the panel (c) end-of-line labels occupy. The
# figure is saved at its declared width, so that strip stays reserved; under
# bbox_inches="tight" the canvas was cropped to the labels and they ended up
# flush against the image edge.
fig = plt.figure(figsize=(FIG_W, 3.17))
gs = fig.add_gridspec(1, 3, width_ratios=[1.50, 0.92, 1.18], wspace=0.52,
                      left=0.088, right=0.875, top=0.912, bottom=0.336)

# ---- (a) MCC distribution, every block under every protocol -------------
ax = fig.add_subplot(gs[0, 0])
from matplotlib.patches import Patch
for i, rg in enumerate(REG):
    for j, b in enumerate(BLK):
        x = i * 4.2 + j
        v = r[(r["regime"] == rg) & (r["block"] == b)]["mcc"].values
        bp = ax.boxplot([v], positions=[x], widths=.74, patch_artist=True, showfliers=False,
                        medianprops=dict(color="k", lw=1.3), whiskerprops=dict(lw=.7),
                        capprops=dict(lw=.7), boxprops=dict(lw=.6))
        bp["boxes"][0].set(facecolor=BC[b], alpha=.55, edgecolor=INK, hatch=BLOCK_HATCH[b])
    if i:
        ax.axvline(i * 4.2 - 1.6, color="#dddddd", lw=.8, zorder=0)
ax.set_xticks([i * 4.2 + 1 for i in range(3)])
ax.set_xticklabels([RLAB[x] for x in REG], fontsize=6)
ax.set_ylabel("Matthews correlation\ncoefficient")
ax.set_ylim(-0.30, 1.05)
ax.set_yticks([0, .2, .4, .6, .8, 1.0])
# Key below the axes, not inside them: at "lower center" it sat in the
# -0.30 to 0.0 band, where the scaffold-split reactivity whisker reaches.
ax.legend(handles=[Patch(facecolor=BC[b], alpha=.55, edgecolor=INK, hatch=BLOCK_HATCH[b],
                         label=BLAB[b]) for b in BLK],
          loc="upper left", bbox_to_anchor=(-0.01, -0.145), ncol=3, fontsize=6,
          handlelength=1.1, handletextpad=.4, columnspacing=1.4,
          borderaxespad=0.0, frameon=False)
ax.text(0.99, 0.02, "higher = better", transform=ax.transAxes,
        va="bottom", ha="right", fontsize=6.0, color=GRAY)
panel_tag(ax, "a", pad_x_pt=-32)

# ---- (b) the paired increment from adding the reactivity block ----------
ax = fig.add_subplot(gs[0, 1])
for i, rg in enumerate(REG):
    piv = pd.read_csv(f"{D}/paired_mcc_{rg.replace('/', '_')}.csv")
    d = (piv["combined"] - piv["physchem"]).values
    parts = ax.violinplot([d], positions=[i], widths=.75, showextrema=False)
    for b in parts["bodies"]:
        b.set(facecolor=BC["combined"], alpha=.30, edgecolor=INK, lw=.6)
    m, c = d.mean(), 1.96 * d.std(ddof=1) / np.sqrt(len(d))
    ax.errorbar([i], [m], yerr=[c], fmt="o", ms=4, color=ACCENT, lw=1.3, capsize=2.5, zorder=5)
    rlab = {"raw600/random": "600 records, random split",
            "dedup387/random": "387 structures, random split",
            "dedup387/scaffold": "387 structures, scaffold split"}[rg]
    row = pw[(pw["regime"] == rlab) & (pw["comparison"] == "Combined - Physicochemical")].iloc[0]
    ax.text(i, 0.52, f"{m:+.3f}", ha="center", fontsize=6.2, color=ACCENT, weight="bold")
ax.axhline(0, color="k", lw=1.0, zorder=1)
ax.set_xticks(range(3)); ax.set_xticklabels([SHORT[x] for x in REG], fontsize=6)
ax.set_xlabel("Compounds / split scheme", fontsize=6.5)
ax.set_ylabel("$\\Delta$MCC,\ncombined − physicochemical")
ax.set_ylim(-0.62, 0.66)
ax.text(2.44, 0.03, "no change", fontsize=6.0, color=INK, ha="right", va="bottom")
panel_tag(ax, "b", pad_x_pt=-34)

# ---- (c) mean Friedman ranks -------------------------------------------
ax = fig.add_subplot(gs[0, 2])
fr = pd.read_csv(f"{T}/T_friedman_by_regime.csv")
key = {"600 records, random split": "raw600/random",
       "387 structures, random split": "dedup387/random",
       "387 structures, scaffold split": "dedup387/scaffold"}
fr["rg"] = fr["regime"].map(key)
fr = fr.set_index("rg").loc[REG]
names = {"cdft": "CDFT reactivity", "physchem": "Physicochemical", "combined": "Combined"}
for j, b in enumerate(BLK):
    ax.plot(range(3), fr[f"mean_rank_{b}"].values, marker=BLOCK_MARKER[b], ms=4.2,
            color=BC[b], lw=1.4, ls="-", dashes=BLOCK_DASH[b][1] or (None, None))
end_of_line_labels(ax, [[0, 1, 2]] * 3, [fr[f"mean_rank_{b}"].values for b in BLK],
                   [names[b] for b in BLK], colors=[BC[b] for b in BLK], dx=0.02, fontsize=6)
ax.set_xticks(range(3)); ax.set_xticklabels([SHORT[x] for x in REG], fontsize=6)
ax.set_xlabel("Compounds / split scheme", fontsize=6.5)
ax.set_ylabel("Mean Friedman rank")
ax.set_ylim(1.2, 3.15); ax.invert_yaxis()
ax.set_xlim(-0.25, 3.05)
ax.text(-0.19, 1.28, "better", fontsize=6.0, color=GRAY, va="center")
panel_tag(ax, "c", pad_x_pt=-34)

fig.savefig(f"{F}/Figure4_block_comparison.png", dpi=400)
fig.savefig(f"{F}/Figure4_block_comparison.pdf")

