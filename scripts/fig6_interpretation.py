# Figure 5 — what the models actually use.
import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
import numpy as np
import pandas as pd
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import Patch

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="open", sizes=(8, 7, 6))
from ureasecdft import config as cfg
D, T, F = str(cfg.PROCESSED), str(cfg.TABLES), str(cfg.FIGURES)
att = pd.read_csv(f"{T}/T_shap_attribution_by_origin.csv")
top = pd.read_csv(f"{T}/T_shap_top25_family_excluded.csv", index_col=0)
stab = pd.read_csv(f"{T}/SI_T8_feature_selection_stability.csv")

from ureasecdft.style import (ORIGIN_COLOR as OC, ACCENT, GRAY, INK, FIG_W,
                         CDFT_NAME, panel_tag)

fig = plt.figure(figsize=(FIG_W, 3.74))
gs = fig.add_gridspec(1, 3, width_ratios=[0.70, 1.12, 1.06], wspace=0.86,
                      left=0.098, right=0.988, top=0.925, bottom=0.146)

# ---- (a) share of total SHAP attribution by descriptor origin -----------
# One stacked bar. The two minority segments are 5.1% and 0.27% wide, so they
# cannot carry their own inline labels; they are named by short horizontal
# leaders into a reserved strip to the right of the bar, which is why the x
# limit runs past 100 while the ticks stop at 100.
ax = fig.add_subplot(gs[0, 0])
BARS = [("included", 1, "family\nincluded"), ("excluded", 0, "family\nwithheld")]
ORDER = ["Physicochemical", "CDFT reactivity", "Chemical family"]
for cond, yy, _ in BARS:
    s = att[att["family"] == cond].set_index("origin")
    left = 0.0
    for o in ORDER:
        if o not in s.index:
            continue
        v = float(s.loc[o, "mean_share_pct"])
        ax.barh([yy], [v], left=left, height=.50, color=OC[o],
                edgecolor="white", lw=.7)
        left += v
    ax.text(float(s.loc["Physicochemical", "mean_share_pct"]) / 2, yy,
            f"{float(s.loc['Physicochemical', 'mean_share_pct']):.1f}%",
            ha="center", va="center", fontsize=7, color="white", weight="bold")
inc = att[att["family"] == "included"].set_index("origin")
pc = float(inc.loc["Physicochemical", "mean_share_pct"])
cd_share = float(inc.loc["CDFT reactivity", "mean_share_pct"])
fam_share = float(inc.loc["Chemical family", "mean_share_pct"])
ax.set_xlim(0, 152); ax.set_ylim(-0.56, 1.56)
ax.set_xticks([0, 25, 50, 75, 100])
ax.set_yticks([b[1] for b in BARS])
ax.set_yticklabels([b[2] for b in BARS], fontsize=6, linespacing=1.3)
ax.set_xlabel("Share of total SHAP attribution (%)")
# Only the upper bar is annotated: the lower bar repeats the same split, and
# that visual identity is the point being made about the family attribute.
ax.annotate(f"CDFT reactivity\n{cd_share:.1f}% $\\pm$ {float(inc.loc['CDFT reactivity', 'ci95_pct']):.1f}",
            xy=(pc + cd_share / 2, 1.21), xytext=(104, 1.45),
            fontsize=6, color=OC["CDFT reactivity"], ha="left", va="center",
            linespacing=1.35,
            arrowprops=dict(arrowstyle="-", lw=.6, color=OC["CDFT reactivity"],
                            shrinkA=0, shrinkB=1))
ax.annotate(f"chemical family\n{fam_share:.2f}%",
            xy=(100, 0.79), xytext=(104, 0.52),
            fontsize=6, color=OC["Chemical family"], ha="left", va="center",
            linespacing=1.35,
            arrowprops=dict(arrowstyle="-", lw=.6, color=OC["Chemical family"],
                            shrinkA=0, shrinkB=1))
ax.spines["left"].set_visible(False)
ax.tick_params(axis="y", length=0)
panel_tag(ax, "a", pad_x_pt=-30)

# ---- (b) top 20 descriptors ---------------------------------------------
# Twenty bars, so that the panel shows exactly the "top 20" its title claims.
# The note about the best CDFT descriptor carries NO leader: rank 25 is off the
# bottom of this list, and an arrow into the plot would land on a
# physicochemical bar and read as labeling it.
ax = fig.add_subplot(gs[0, 1])
t = top.head(20).iloc[::-1]
ax.barh(range(len(t)), t["mean_shap_share"].values,
        color=[OC[o] for o in t["origin"]], height=.70)
ax.set_yticks(range(len(t))); ax.set_yticklabels(t.index, fontsize=6.0)
ax.set_xlabel("Mean SHAP attribution (%)")
ax.text(t["mean_shap_share"].values[-1] - 0.25, len(t) - 1.0,
        f"{t['mean_shap_share'].values[-1]:.1f}%", fontsize=6, va="center",
        ha="right", color="white", weight="bold")
ax.text(2.7, 3.2, "all 20 bars are\nphysicochemical; the best\nCDFT descriptor ranks 25th",
        fontsize=6, color=OC["CDFT reactivity"], ha="left", va="center",
        linespacing=1.45)
ax.set_xlim(0, 11.5)
ax.margins(y=0.015)
panel_tag(ax, "b", pad_x_pt=-34)

# ---- (c) which reactivity descriptors survive selection ----------------
# Paired dots rather than paired hatched bars: with 18 rows the hatching was
# the dominant visual texture of the figure, and an exact zero needed a faked
# stub to show at all. A dot at zero is already a visible mark.
ax = fig.add_subplot(gs[0, 2])
cd = (stab[stab["origin"] == "CDFT reactivity"]
      .pivot_table(index="feature", columns="method", values="selection_frequency")
      .fillna(0.0).sort_values("Boruta"))
y = np.arange(len(cd))
ax.hlines(y, np.minimum(cd["Boruta"], cd["XGBoost gain"]),
          np.maximum(cd["Boruta"], cd["XGBoost gain"]),
          color=GRAY, lw=.8, zorder=2)
ax.scatter(cd["XGBoost gain"].values, y, s=14, marker="s", facecolor="white",
           edgecolor=INK, lw=.8, zorder=3, label="XGBoost gain")
ax.scatter(cd["Boruta"].values, y, s=16, marker="o", color=OC["CDFT reactivity"],
           lw=0, zorder=4, label="Boruta")
ax.vlines(0.5, -0.62, len(cd) + 0.9, color=INK, ls=":", lw=.8, zorder=1)
ax.set_yticks(y)
ax.set_yticklabels([CDFT_NAME.get(i, i) for i in cd.index], fontsize=6.0)
ax.set_xlabel("Fraction of runs retained")
ax.set_xlim(-0.04, 1.04)
ax.set_xticks([0, 0.25, 0.5, 0.75, 1.0])
# Reserve an empty row at the top for the key, and put the threshold label
# below the lowest data row, so neither sits on a datum.
ax.set_ylim(-1.5, len(cd) + 1.0)
ax.scatter([0.02], [len(cd) + 0.45], s=16, marker="o",
           color=OC["CDFT reactivity"], lw=0, clip_on=False)
ax.text(0.075, len(cd) + 0.45, "Boruta", fontsize=6, va="center",
        color=OC["CDFT reactivity"])
ax.scatter([0.40], [len(cd) + 0.45], s=14, marker="s", facecolor="white",
           edgecolor=INK, lw=.8, clip_on=False)
ax.text(0.455, len(cd) + 0.45, "XGBoost gain", fontsize=6, va="center", color=INK)
ax.text(0.5, -1.05, "half of runs", fontsize=6, color=GRAY,
        ha="center", va="center")
panel_tag(ax, "c", pad_x_pt=-86)

fig.savefig(f"{F}/Figure5_interpretation.png", dpi=400)
fig.savefig(f"{F}/Figure5_interpretation.pdf")

print("\nCDFT retention:\n", cd.round(3).to_string())
