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
    plot_errorband(ax, lo, colors["lo"])
    plot_errorband(ax, nlo, colors["nlo"])
    plot_errorband(ax, full, colors["full"])

    # Add labels for legend
    if labels is not None:
        for key, label in labels.items():
            ax.plot([], [], color=colors[key], label=label)
        ax.legend(framealpha=0)

# =========================
# Plot x5/y5 bands in slices
# =========================
def plot_bands(full_base, bands_dict, *,
               xlabel, ylabel, title, savename, outdir, colors,
               slice_name="y5",
               slice_range=None):   # optional: (min, max)

    if full_base is None or len(bands_dict) == 0:
        return

    fig, ax = plt.subplots(figsize=(6,5))

    # Plot full distribution
    ax.plot(full_base[:,0], full_base[:,1],
            color=colors["full"], label="Full range")

    # Determine bands
    band_indices = sorted(bands_dict.keys())
    n_bands = len(band_indices)

    # Determine slice range
    if slice_range is not None:
        min_val, max_val = slice_range
    else:
        # automatic from data (robust)
        all_vals = np.concatenate([b[:,0] for b in bands_dict.values()])
        min_val, max_val = np.min(all_vals), np.max(all_vals)

    # Compute edges for labeling
    edges = np.linspace(min_val, max_val, n_bands + 1)

    # Color map
    cmap = plt.cm.viridis

    # Plot bands
    for idx, i in enumerate(band_indices):
        band = bands_dict[i]
        color = cmap(idx / max(n_bands - 1, 1))

        # Use robust indexing for edges
        low  = edges[idx]
        high = edges[idx+1]

        label = f"{low:.3f} < {slice_name} < {high:.3f}"

        ax.plot(band[:,0], band[:,1],
                color=color,
                label=label)

    # Labels, title, and axis limits
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)
    ax.set_xlim(full_base[:,0].min(), full_base[:,0].max())

    # Legend
    ax.legend(fontsize=8, ncol=2)

    # Save figure
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
    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches="tight")

# =========================
# Plot K-factor
# =========================
def plot_K(ax, K, colors_K, xlabel):
    plot_lo_nlo_full(ax, lo=K, nlo=K, full=K,
                     colors=dict(lo=colors_K, nlo=colors_K, full=colors_K),
                     labels=dict(full="K-factor"))
    style_sci_x(ax, xlabel, r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$", None, yscale="linear")

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