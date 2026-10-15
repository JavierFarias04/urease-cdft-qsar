"""Figure style for the manuscript's figure set — self-contained.

The published figures were produced with an external styling layer. The exact
matplotlib rcParams it applied are reproduced verbatim below, so this
repository renders the same figures with no dependency beyond matplotlib.

One parameter is deliberately changed from the published run: `savefig.bbox`
was `"tight"`, which crops every canvas to its own content, so the six figures
were saved at widths between 6.67 and 7.26 in instead of the declared 7.20 in.
Here it is `None`, so FIG_W is honored and the set is dimensionally uniform.
See "Known differences on re-running" in the README.
"""
import matplotlib as mpl
import matplotlib.pyplot as plt

# Exact rcParams of the published run, plus the savefig.bbox correction.
RC = {
    "axes.labelsize": 8.0, "axes.linewidth": 0.6,
    "axes.spines.right": False, "axes.spines.top": False,
    "axes.titlelocation": "left", "axes.titlesize": 8.0,
    "figure.dpi": 200.0, "font.size": 8.0,
    "legend.fontsize": 7.0, "legend.frameon": False,
    "lines.linewidth": 1.2, "patch.linewidth": 0.6,
    "pdf.fonttype": 42, "ps.fonttype": 42,
    "savefig.bbox": None,          # published run used "tight"; see module docstring
    "savefig.dpi": 300.0,
    "xtick.labelsize": 6.0, "xtick.major.size": 3.0, "xtick.major.width": 0.6,
    "ytick.labelsize": 6.0, "ytick.major.size": 3.0, "ytick.major.width": 0.6,
}


def apply_figure_style(frame="open", sizes=(8.0, 7.0, 6.0)):
    """Apply the manuscript figure style. `frame`: 'open' hides the top and
    right spines, 'none' hides all four. `sizes` is the (base, secondary, tick)
    font ladder."""
    mpl.rcParams.update(RC)
    base, sec, tick = sizes
    mpl.rcParams.update({
        "font.size": base, "axes.labelsize": base, "axes.titlesize": base,
        "legend.fontsize": sec, "xtick.labelsize": tick, "ytick.labelsize": tick,
    })
    off = {"open": ("top", "right"), "none": ("top", "right", "left", "bottom")}[frame]
    mpl.rcParams.update({f"axes.spines.{s}": False for s in off})
    return dict(mpl.rcParams)


def end_of_line_labels(ax, xs, ys, labels, colors=None, dx=0.01, fontsize=None):
    """Label each line series at its right end instead of using a legend box."""
    if fontsize is None:
        fontsize = plt.rcParams["font.size"]
    if colors is None:
        colors = [None] * len(labels)
    span = ax.get_xlim()[1] - ax.get_xlim()[0]
    for x, y, lab, c in zip(xs, ys, labels, colors):
        ax.text(x[-1] + dx * span, y[-1], lab, color=c, va="center", ha="left",
                fontsize=fontsize)


FONT_MIN = 6.0          # absolute floor: no text in the set below this size
FS_BASE, FS_SEC, FS_TICK = 8.0, 7.0, 6.0
FS_PANEL = 9.0          # panel letters: the one permitted exception to the ladder
FIG_W = 7.20            # every figure in the set is this wide



def panel_tag(ax, letter, pad_x_pt=-30.0, fs=FS_PANEL):
    """Place a bold panel letter with its cap-top aligned to the title's top.

    The vertical offset is measured in POINTS from the rendered title box, not
    in axes fractions. An axes-fraction offset is the defect this replaces: the
    same fraction is a different physical gap in every panel whose height
    differs, so letters drifted relative to their titles both within a figure
    (unequal gridspec rows) and across the set (1x3 versus 2x2 layouts).
    """
    fig = ax.figure
    fig.canvas.draw()
    rr = fig.canvas.get_renderer()
    px_to_pt = 72.0 / fig.dpi
    if ax.title.get_text().strip():
        dy_pt = (ax.title.get_window_extent(rr).y1
                 - ax.get_window_extent(rr).y1) * px_to_pt
    else:
        dy_pt = 12.0
    return ax.annotate(letter, xy=(0, 1), xycoords="axes fraction",
                       xytext=(pad_x_pt, dy_pt), textcoords="offset points",
                       fontsize=fs, fontweight="bold", color="#2b2b2b",
                       va="top", ha="left", annotation_clip=False)


# Display names for the CDFT block. SYM is the symbol alone, for a matrix axis
# where the row and column order already identifies the quantity; NAME is the
# spelled-out name plus symbol, for any axis a reader meets cold. No figure in
# the set may fall through to a raw dataframe column name.
CDFT_SYM = {
    "homo_energy": "$E_{HOMO}$", "homo_.1_energy": "$E_{HOMO-1}$",
    "lumo_energy": "$E_{LUMO}$", "gap": "gap",
    "potential_ionization": "IP", "affinity_electronic": "EA",
    "chemical_potential": "$\\mu$", "hardness": "$\\eta$", "softness": "$S$",
    "electronegativity": "$\\chi$", "electrophilicity": "$\\omega$",
    "nucleophilicity": "$N$", "acceptor_power": "$\\omega^{+}$",
    "donor_power": "$\\omega^{-}$",
    "net_electrophilicity": "$\\Delta\\omega^{\\pm}$",
    "hyperhardness": "$\\gamma$",
    "chemical_potential_pos": "$\\mu^{+}$",
    "chemical_potential_neg": "$\\mu^{-}$",
    "dual_reactivity_index": "dual index",
}
CDFT_NAME = {
    "homo_energy": "HOMO energy  $E_{HOMO}$",
    "homo_.1_energy": "HOMO$-$1 energy  $E_{HOMO-1}$",
    "lumo_energy": "LUMO energy  $E_{LUMO}$",
    "gap": "HOMO$-$LUMO gap",
    "potential_ionization": "ionization potential  IP",
    "affinity_electronic": "electron affinity  EA",
    "chemical_potential": "chemical potential  $\\mu$",
    "hardness": "hardness  $\\eta$",
    "softness": "softness  $S$",
    "electronegativity": "electronegativity  $\\chi$",
    "electrophilicity": "electrophilicity  $\\omega$",
    "nucleophilicity": "nucleophilicity  $N$",
    "acceptor_power": "electroaccepting power  $\\omega^{+}$",
    "donor_power": "electrodonating power  $\\omega^{-}$",
    "net_electrophilicity": "net electrophilicity  $\\Delta\\omega^{\\pm}$",
    "hyperhardness": "hyperhardness  $\\gamma$",
    "chemical_potential_pos": "chemical potential  $\\mu^{+}$",
    "chemical_potential_neg": "chemical potential  $\\mu^{-}$",
    "dual_reactivity_index": "dual reactivity index",
}

# descriptor blocks
BLOCK_COLOR = {"physchem": "#0072B2", "cdft": "#E69F00", "combined": "#2b2b2b"}
BLOCK_LABEL = {"physchem": "Physicochemical", "cdft": "CDFT reactivity", "combined": "Combined"}
BLOCK_MARKER = {"physchem": "o", "cdft": "s", "combined": "D"}
BLOCK_HATCH = {"physchem": "", "cdft": "///", "combined": "..."}
BLOCK_DASH = {"physchem": (0, ()), "cdft": (0, (4, 1.6)), "combined": (0, (1, 1.2))}

# chemical families
FAM_ORDER = ["Tiazol", "Triazol", "Coumarin", "Otros"]
FAM_LABEL = {"Tiazol": "Thiazoles", "Triazol": "Triazoles",
             "Coumarin": "Coumarins", "Otros": "Other scaffolds"}
FAM_SHORT = {"Tiazol": "Thiazoles", "Triazol": "Triazoles",
             "Coumarin": "Coumarins", "Otros": "Other"}
FAM_COLOR = {"Tiazol": "#0072B2", "Triazol": "#E69F00",
             "Coumarin": "#2b2b2b", "Otros": "#9ad0f0"}
FAM_MARKER = {"Tiazol": "o", "Triazol": "s", "Coumarin": "^", "Otros": "D"}
FAM_HATCH = {"Tiazol": "", "Triazol": "///", "Coumarin": "...", "Otros": "xxx"}

# split schemes
SPLIT_COLOR = {"random": "#0072B2", "scaffold": "#D55E00"}
SPLIT_MARKER = {"random": "o", "scaffold": "s"}

# descriptor origin (Figure 6) reuses the block palette
ORIGIN_COLOR = {"Physicochemical": BLOCK_COLOR["physchem"],
                "CDFT reactivity": BLOCK_COLOR["cdft"],
                "Chemical family": "#7f7f7f"}
ORIGIN_HATCH = {"Physicochemical": "", "CDFT reactivity": "///", "Chemical family": "xxx"}

INK, GRAY, ACCENT = "#2b2b2b", "#6e6e6e", "#D55E00"
ZERO_STUB = 0.006       # hairline length for an exact-zero bar, in axis units
