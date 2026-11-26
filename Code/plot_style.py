"""Common plot style for AMBER_RadCor figures (large labels, friendly colours).

Usage:
    from plot_style import apply_style, COLORS
    apply_style()

Palette builds on the repo's existing colours (analyze_mcmule4vecs.py:
LO green #19866B, NLO pink #C879A3, R purple #9772B5, diff blue #477E9E),
slightly brighter, plus a warm orange and a soft highlight colour.
"""

import matplotlib
import matplotlib.pyplot as plt

COLORS = {
    "LO": "#19866B",        # green
    "NLO": "#C879A3",       # pink
    "NLO0": "#C879A3",
    "F": "#477E9E",         # blue
    "R": "#E6912E",         # warm orange
    "A": "#9772B5",         # purple
    "SUM": "#C2397F",       # strong pink (sums like F+R)
    "DIFF": "#477E9E",      # blue (differences, ratios)
    "PULL": "#19866B",      # green
    "BAND": "#9FB8C8",      # light blue-grey (error bands)
    "HIGHLIGHT": "#F6D55C", # soft yellow (highlighted regions)
    "REF": "#666666",       # reference lines (0, +-1 sigma)
}

# Two contrasting colours for "A vs B" comparisons
PAIR = ("#477E9E", "#C2397F")


def apply_style():
    matplotlib.rcParams["text.usetex"] = False
    plt.rcParams.update({
        "font.size": 15,
        "axes.titlesize": 17,
        "axes.labelsize": 17,
        "xtick.labelsize": 14,
        "ytick.labelsize": 14,
        "legend.fontsize": 14,
        "figure.titlesize": 18,
        "lines.linewidth": 2.0,
        "lines.markersize": 5,
        "axes.grid": True,
        "grid.alpha": 0.30,
        "grid.linestyle": ":",
        "axes.spines.top": False,
        "axes.spines.right": False,
        "legend.frameon": False,
        "savefig.dpi": 300,
        "savefig.bbox": "tight",
    })
