import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.ticker import ScalarFormatter
from pymule import *

# =========================
# Plot an error band for a histogram
# =========================
def plot_errorband(ax, hist, color):
    plt.sca(ax)
    for line in errorband(hist):
        line.set_color(color)

# =========================
# Plot LO, NLO, and full distributions with error bands
# =========================
def plot_lo_nlo_full(ax, lo, nlo, full, colors, labels=None):
    if lo is not None:
        plot_errorband(ax, lo, colors["lo"])
    if nlo is not None:
        plot_errorband(ax, nlo, colors["nlo"])
    if full is not None:
        plot_errorband(ax, full, colors["full"])

    if labels is not None:
        for key, label in labels.items():
            ax.plot([], [], color=colors[key], label=label)
        ax.legend(framealpha=0)

# =========================
# Plot x5/y5 bands in slices
# =========================
def plot_bands(bands_dict, *,
               xlabel, ylabel, title, savename, outdir, colors,
               slice_name="y5", slice_range=(-0.09, 0.09)):
    """
    Plottet Bänder einer Verteilung in Slices.
    Keine Hauptverteilung mehr — die Bänder sind die Story.
    """
    fig, ax = plt.subplots(figsize=(6, 5))

    if not bands_dict:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        save_figure(fig, savename, outdir=outdir)
        plt.close(fig)
        return

    band_indices = sorted(bands_dict.keys())
    n_bands = len(band_indices)
    min_val, max_val = slice_range
    edges = np.linspace(min_val, max_val, n_bands + 1)
    cmap = plt.cm.viridis

    for idx, i in enumerate(band_indices):
        band = bands_dict[i]
        if band is None or len(band) == 0:
            continue
        # Nur finite Werte
        mask = np.isfinite(band[:, 0]) & np.isfinite(band[:, 1])
        band = band[mask]
        if len(band) == 0:
            continue

        color = cmap(idx / max(n_bands - 1, 1))
        low, high = edges[idx], edges[idx + 1]
        label = f"{low:.3f} < {slice_name} < {high:.3f}"
        ax.plot(band[:, 0], band[:, 1], color=color, label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.legend(fontsize=8, ncol=2, framealpha=0)

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)

# =========================
# Styling for scientific x-axis
# =========================
def style_sci_x(ax, xlabel, ylabel, title, yscale="log", sharex=False):
    style_axis(ax, xlabel=xlabel, ylabel=ylabel, title=title, yscale=yscale, legend=False)

# =========================
# General axis styling
# =========================
def style_axis(ax, xlabel=None, ylabel=None, title=None, xscale="linear", yscale="linear", legend=False):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid(True, alpha=0.7, ls="dotted")

    if legend:
        ax.legend(framealpha=0)

# =========================
# Figure helpers
# =========================
def create_figure(nrows=1, ncols=1, figsize=(8,5), font_size=7, sharex=False, gridspec_kw=None):
    plt.rcParams.update({
        "font.size": font_size,
        "font.family": "serif",
        "text.usetex": False,
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, gridspec_kw=gridspec_kw)

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    return fig, axes

# =========================
# Add a secondary x-axis
# =========================
def add_secondary_xaxis(ax, transform_func, xlabel=None, scale="linear"):
    ax_sec = ax.twiny()
    ax_sec.set_xlim(transform_func(np.array(ax.get_xlim())))
    ax_sec.set_xscale(scale)
    if xlabel:
        ax_sec.set_xlabel(xlabel, fontsize=9)
    return ax_sec

# =========================
# Save figure in PDF and PNG
# =========================
def save_figure(fig, name, outdir, dpi=700):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    #fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches="tight")

# =========================
# Plot K-factor
# =========================
def plot_K(ax, K, colors_K, xlabel):
    plt.sca(ax)
    for line in errorband(K):
        line.set_color(colors_K)
    ax.plot([], [], color=colors_K, label="K-factor")
    style_sci_x(ax, xlabel, r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$", None, yscale="linear")
    ax.legend(framealpha=0)

# =========================
# Write values to file
# =========================
def write_file_with_values(filename, parameter_array, xlabel, parameter_label):
    with open(filename, "w") as f:
        f.write(f"{xlabel}    {parameter_label}    Error \n")
        for i, row in enumerate(parameter_array):
            bin_center = row[0]
            value      = row[1]
            error      = row[2]
            f.write(f"{bin_center: .6e}  {value: .6e}  {error: .6e}\n")