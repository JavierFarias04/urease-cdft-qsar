# Figure 3 — the CDFT block is a low-dimensional, partly unphysical description.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
import json
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Rectangle

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="open", sizes=(8, 7, 6))
from ureasecdft import config as cfg
D, T, F = str(cfg.PROCESSED), str(cfg.TABLES), str(cfg.FIGURES)
raw = pd.read_csv(f"{D}/compounds_raw.csv")
meta = json.load(open(f"{D}/build_meta.json"))
spec = np.load(f"{D}/pca_spectra.npy", allow_pickle=True).item()
CDFT = [c for c in meta["cdft_cols"] if raw[c].nunique() > 1]

from ureasecdft.style import (FAM_ORDER, FAM_SHORT as FAM_LABEL, FAM_COLOR as CB, FAM_MARKER,
                         BLOCK_COLOR, BLOCK_DASH, ACCENT, GRAY, INK, FIG_W,
                         CDFT_SYM, panel_tag)

fig = plt.figure(figsize=(FIG_W, 3.20))
gs = fig.add_gridspec(1, 3, width_ratios=[1.22, 1.0, 1.08], wspace=0.72,
                      left=0.082, right=0.988, top=0.913, bottom=0.189)

# ---- (a) |r| among the 18 informative CDFT descriptors -------------------
ax = fig.add_subplot(gs[0, 0])
C = raw[CDFT].astype(float).corr().abs()
# aspect="auto", not "equal": a square matrix shrinks its own axes inside the
# allotted column, which both wastes the panel and pulls the left-aligned title
# and the panel letter inward, out of line with the other two panels.
im = ax.imshow(C.values, cmap="magma_r", vmin=0, vmax=1, aspect="auto")
lab = [CDFT_SYM[c] for c in CDFT]
ax.set_xticks(range(len(lab))); ax.set_xticklabels(lab, rotation=90, fontsize=6)
ax.set_yticks(range(len(lab))); ax.set_yticklabels(lab, fontsize=6)
ax.tick_params(length=0)
for s in ax.spines.values():
    s.set_visible(False)
cb = fig.colorbar(im, ax=ax, fraction=.045, pad=.03,
                  ticks=[0, 0.25, 0.5, 0.75, 1.0])
cb.set_label("|Pearson $r$|", fontsize=6); cb.ax.tick_params(labelsize=6)
panel_tag(ax, "a", pad_x_pt=-34)

# ---- (b) cumulative PCA variance ----------------------------------------
ax = fig.add_subplot(gs[0, 1])
sty = {"CDFT reactivity": (BLOCK_COLOR["cdft"], 1.7, BLOCK_DASH["cdft"]),
       "Physicochemical": (BLOCK_COLOR["physchem"], 1.3, BLOCK_DASH["physchem"]),
       "Combined": (BLOCK_COLOR["combined"], 1.3, BLOCK_DASH["combined"])}
N95 = {}
for k, (c, lw, ls) in sty.items():
    cum = 100 * np.cumsum(spec[k])
    x = np.arange(1, len(cum) + 1)
    ax.plot(x, cum, color=c, lw=lw, ls=ls, zorder=3 if k.startswith("CDFT") else 2)
    n = int(np.argmax(cum >= 95) + 1)
    N95[k] = n
    ax.plot([n], [cum[n - 1]], marker="o", ms=3.2, color=c, mec="white", mew=.6,
            zorder=5)
ax.axhline(95, color="k", ls=":", lw=.8, zorder=1)
ax.set_xscale("log"); ax.set_xlim(1, 60); ax.set_ylim(0, 106)
ax.set_yticks([0, 20, 40, 60, 80, 95])   # the 95 tick labels the dotted line
ax.set_xticks([1, 3, 10, 30, 60]); ax.set_xticklabels(["1", "3", "10", "30", "60"])
ax.set_xlabel("Principal components"); ax.set_ylabel("Cumulative variance (%)")
# The three "components needed for 95%" counts are carried by a key in the
# empty lower-right of the panel. They were previously annotations tied to the
# curves by long diagonal leaders, which read as three further data series.
from matplotlib.lines import Line2D
ax.legend(handles=[Line2D([0], [0], color=sty[k][0], lw=sty[k][1], ls=sty[k][2],
                          marker="o", ms=3.2, mec="white", mew=.6,
                          label=f"{k}: {N95[k]}")
                   for k in ["CDFT reactivity", "Physicochemical", "Combined"]],
          loc="lower right", bbox_to_anchor=(1.015, -0.025), fontsize=6,
          frameon=False, handlelength=2.2, handletextpad=.5, labelspacing=.3)
panel_tag(ax, "b", pad_x_pt=-34)

# ---- (c) frontier orbital eigenvalues: bound vs unbound ------------------
ax = fig.add_subplot(gs[0, 2])
H = raw["homo_energy"].astype(float).values
L = raw["lumo_energy"].astype(float).values
ax.add_patch(Rectangle((0, -0.25), 0.30, 0.55, color=ACCENT, alpha=.10, lw=0, zorder=0))
for f in FAM_ORDER:
    m = raw["tipo_familia"].values == f
    ax.scatter(H[m], L[m], s=6.5, marker=FAM_MARKER[f], color=CB[f], alpha=.65, lw=0,
               label=FAM_LABEL[f], zorder=3)
ax.axvline(0, color="k", lw=.7, zorder=2)
ax.set_xlabel("$E_{HOMO}$ (hartree)"); ax.set_ylabel("$E_{LUMO}$ (hartree)")
ax.text(0.095, -0.145, f"$E_{{HOMO}}>0$\n{(H > 0).sum()} records\n(62 coumarins)",
        fontsize=6, color=ACCENT, ha="center", va="top", linespacing=1.3)
ax.legend(loc="upper left", fontsize=6, handletextpad=.25, borderaxespad=.35,
          labelspacing=.22, markerscale=1.4, frameon=False)
ax.set_ylim(-0.235, 0.30); ax.set_xlim(-0.345, 0.185)
panel_tag(ax, "c", pad_x_pt=-34)

fig.savefig(f"{F}/Figure3_cdft_block.png", dpi=400)
fig.savefig(f"{F}/Figure3_cdft_block.pdf")

