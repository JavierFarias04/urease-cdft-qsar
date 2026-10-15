import sys as _sys, os as _os
_sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), "..", "src"))
from ureasecdft import config as cfg
# Figure 1 — study design: what is compared, and under which evaluation protocols.
import numpy as np
import matplotlib as mpl
import matplotlib.pyplot as plt
from matplotlib.patches import FancyBboxPatch, FancyArrowPatch

from ureasecdft.style import apply_figure_style, end_of_line_labels
apply_figure_style(frame="none", sizes=(8, 7, 6))
F = str(cfg.FIGURES)
from ureasecdft.style import BLOCK_COLOR, BLOCK_DASH, INK, GRAY, FIG_W
# Version 2. The figure previously ran two encodings in one hue family: color
# carried block identity in row 2 (orange = reactivity) and ALSO carried
# emphasis in rows 1 and 4 (vermillion = the protocol the conclusions rest on).
# Two warm oranges in adjacent rows invited the scaffold-split box to be read
# as a reactivity box. Here the two registers are separated: color encodes
# block identity and nothing else, and emphasis is carried by stroke weight and
# type weight alone. The combined block also takes the dotted stroke it has in
# Figures 4 and 6, since its near-black is otherwise the border color of every
# neutral box in the schematic.
EC_BLOCK = {"physchem": BLOCK_COLOR["physchem"], "cdft": BLOCK_COLOR["cdft"],
            "combined": BLOCK_COLOR["combined"]}
FC_BLOCK = {"physchem": "#eaf2f9", "cdft": "#fdf3e2", "combined": "#efefef"}

fig, ax = plt.subplots(figsize=(FIG_W, 4.9))
# The schematic is the whole figure, so the axes fills the canvas instead of
# sitting inside default subplot margins. Type is scaled by the same factor as
# the enlarged axes, so box-to-text proportions are unchanged and the printed
# schematic is more legible at column width.
fig.subplots_adjust(left=0.012, right=0.988, bottom=0.012, top=0.988)
# Limits run slightly past the 0-100 design grid on every side. The rounded
# boxes carry a 0.35-unit boxstyle pad, so a box declared flush to the grid
# edge would otherwise have its border stroke clipped by the axes clip box.
ax.set_xlim(-1.5, 101.5); ax.set_ylim(-4.0, 100.5); ax.axis("off")


def box(x, y, w, h, text, fc="white", ec=INK, fs=8.0, lw=.8, weight="normal", tc=INK,
        ls="-"):
    ax.add_patch(FancyBboxPatch((x, y), w, h, boxstyle="round,pad=0.35,rounding_size=1.4",
                                fc=fc, ec=ec, lw=lw, ls=ls, zorder=3))
    ax.text(x + w / 2, y + h / 2, text, ha="center", va="center", fontsize=fs,
            zorder=4, weight=weight, color=tc, linespacing=1.35)


def arrow(x1, y1, x2, y2, c=GRAY, lw=.9, style="-|>"):
    ax.add_patch(FancyArrowPatch((x1, y1), (x2, y2), arrowstyle=style, mutation_scale=7,
                                 color=c, lw=lw, shrinkA=1, shrinkB=1, zorder=2))


def stage(x, y, t):
    # Stage headings are neutral: they label structure, not a variable.
    ax.text(x, y, t, fontsize=7.9, color=INK, weight="bold", ha="left", va="center")


FS = 7.5  # schematic body type; FONT_MIN = 6.0 is the floor, not a target
# ---- row 1: curation --------------------------------------------------------
stage(0, 97.5, "DATA")
# Solid stroke = a processing step; dashed gray = an assumption imposed on the
# data rather than a step performed on it.
r1 = [("Urease inhibition data\nfrom the literature\n(IC$_{50}$, Web of Science)", INK, "-", INK),
      ("Structures built and\nprotonated at pH 7.2\n(2D sketch → Epik)", INK, "-", INK),
      ("600 records standardized\nand deduplicated →\n387 unique structures", INK, "-", INK),
      ("IC$_{50}$ not reported ⇒ no\ninhibition detected, censored\nat 1 mM (24% of records)",
       GRAY, (0, (3.2, 1.8)), GRAY)]
for i, (t, ec, ls, tc) in enumerate(r1):
    x = i * 25.5
    box(x, 84.0, 23.0, 10.0, t, fc="white", ec=ec, fs=FS, ls=ls, tc=tc)
    if i:
        arrow(x - 2.5, 89.0, x - 0.2, 89.0)

# ---- row 2: the three descriptor blocks ------------------------------------
stage(0, 78.0, "DESCRIPTOR BLOCKS")
ax.text(23.5, 78.0, "identical compounds in every block", fontsize=FS, color=INK, va="center")
r2 = [("physchem", "Physicochemical and topological\n287 CDK descriptors\nconstitution · topology · logP"),
      ("cdft", "CDFT global reactivity\n19 Koopmans descriptors\nfrom B3LYP-D3/6-31G"),
      ("combined", "Combined\n306 descriptors\nunion of both blocks")]
for i, (k, t) in enumerate(r2):
    x = 2 + i * 33.0
    box(x, 63.5, 31.0, 10.5, t, fc=FC_BLOCK[k], ec=EC_BLOCK[k], fs=FS, lw=1.1,
        ls=BLOCK_DASH[k] if k == "combined" else "-")
ax.text(50, 59.6, "chemical family (coumarin / thiazole / triazole / other) added to, "
                  "or withheld from, every block",
        ha="center", fontsize=FS, color=INK, style="italic")

# ---- row 3: modeling factorial -------------------------------------------
stage(0, 54.0, "MODELING")
ax.text(16.5, 54.0, "3 blocks × 2 family conditions × 6 activity definitions × 3 selection "
                    "regimes × 7 classifiers × 10 repeats", fontsize=FS, color=INK, va="center")
r3 = [("7 classifiers\nrandom forest · SVM (RBF)\ndecision tree · XGBoost\nkNN · naïve Bayes\n"
       "logistic regression", "white", INK),
      ("6 activity definitions\nsingle cutoff\n5 / 10 / 25 µM\ntwo cutoff 5, 10, 25 vs 50 µM\n"
       "(intermediate zone dropped)", "white", INK),
      ("3 feature-selection\nregimes\nnone\nBoruta\nXGBoost gain", "white", INK),
      ("Fixed hyperparameters,\nidentical across blocks\nSMOTE applied inside\ntraining folds only",
       "white", INK)]
for i, (t, fc, ec) in enumerate(r3):
    x = i * 25.5
    box(x, 33.0, 23.0, 16.5, t, fc=fc, ec=ec, fs=FS)

# ---- row 4: the three evaluation protocols --------------------------------
# Every model from the factorial above is judged under all three protocols, so
# the fan-out runs through one collecting bus. Drawing an arrow from under each
# modeling box, as before, implied that each protocol came from a different
# part of the factorial.
BUS_Y = 30.6
for cx in [i * 25.5 + 11.5 for i in range(4)]:
    ax.plot([cx, cx], [32.65, BUS_Y], color=GRAY, lw=.9, zorder=2)
ax.plot([11.5, 88.0], [BUS_Y, BUS_Y], color=GRAY, lw=.9, zorder=2)

stage(0, 26.8, "EVALUATION")
# Emphasis is monotone: gray and thin for the protocol shown only for
# comparison, black and thin for the intermediate control, black and heavy with
# bold type for the protocol the conclusions rest on. No hue is involved, so
# nothing here competes with the block colors one row above.
prot = [(0, "600 records\nrandom 80/20 split", "reproduces the\npublished protocol",
         GRAY, 0.8, "normal"),
        (33.5, "387 structures\nrandom 80/20 split", "isolates the effect of\nduplicate records",
         INK, 0.8, "normal"),
        (67, "387 structures\nBemis–Murcko scaffold split", "no scaffold shared between\ntraining and test",
         INK, 2.1, "bold")]
for x, t, sub, c, lw, wt in prot:
    box(x, 12.0, 32.0, 10.0, t, ec=c, lw=lw, fs=FS, fc="white", tc=c, weight=wt)
    ax.text(x + 16, 8.0, sub, ha="center", va="center", fontsize=FS, color=c, linespacing=1.3)
    arrow(x + 16, BUS_Y, x + 16, 22.4, c=c, lw=1.7 if wt == "bold" else .9)

# Footer. Previously a single line inside a gray panel: at 6 pt the string was
# about 6.7 in wide inside axes about 5.6 in wide, so it overflowed its own
# panel on the left and ran to the canvas edge on the right, which is what made
# the figure look cropped. Wrapped to two lines and left-aligned on the same
# x = 0 as the stage headings; the panel is gone because it framed nothing.
ax.text(0, 0.6,
        "95% confidence intervals over repeats · blocks compared by paired Friedman\n"
        "and Wilcoxon tests on matched configurations · SHAP attribution by descriptor origin",
        ha="left", va="center", fontsize=FS, color=INK, linespacing=1.55)

fig.savefig(f"{F}/Figure1_workflow_v2.png", dpi=400)
fig.savefig(f"{F}/Figure1_workflow_v2.pdf")

