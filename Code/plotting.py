import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.ticker import ScalarFormatter
from pymule import *

def plot_errorband(ax, hist, color):
    plt.sca(ax)
    for line in errorband(hist):
        line.set_color(color)

def plot_lo_nlo_full(ax, lo, nlo, full, colors, labels=None):
    plot_errorband(ax, lo, colors["lo"])
    plot_errorband(ax, nlo, colors["nlo"])
    plot_errorband(ax, full, colors["full"])

    if labels is not None:
        for key, label in labels.items():
            ax.plot([], [], color=colors[key], label=label)
        ax.legend(framealpha=0)


# =========================
# Axis styling
# =========================
def style_sci_x(
    ax,
    xlabel,
    ylabel,
    title,
    yscale="log",
    sharex=False,
):
    style_axis(
        ax,
        xlabel=xlabel,
        ylabel=ylabel,
        title=title,
        yscale=yscale,
        legend=False,
    )

    ax.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
    ax.ticklabel_format(style="sci", axis="x", scilimits=(-3, 3))

    if sharex:
        ax.tick_params(labelbottom=False)


def style_axis(
    ax,
    xlabel=None,
    ylabel=None,
    title=None,
    xscale="linear",
    yscale="linear",
    legend=False,
):
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
def create_figure(
    nrows=1,
    ncols=1,
    figsize=(8, 5),
    font_size=7,
):
    plt.rcParams.update({
        "font.size": font_size,
        "font.family": "serif",
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])

    return fig, axes


def add_secondary_xaxis(ax, transform_func, xlabel=None, scale="linear"):
    ax_sec = ax.twiny()
    ax_sec.set_xlim(transform_func(np.array(ax.get_xlim())))
    ax_sec.set_xscale(scale)

    if xlabel:
        ax_sec.set_xlabel(xlabel, fontsize=9)

    return ax_sec


def save_figure(fig, name, outdir, dpi=700):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches="tight")

def plot_K(ax, K, colors_K, xlabel, title):
    plot_lo_nlo_full(
        ax,
        lo=K, nlo=K, full=K,
        colors=dict(lo=colors_K, nlo=colors_K, full=colors_K),
        labels=dict(full="K-factor"),
    )
    style_sci_x(ax, xlabel, r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$", title)


# =========================
# Write values in file
# =========================
def write_file_with_values(filename, parameter_array, xlabel, parameter_label):
    with open(filename, "w") as f:
        f.write(f"{xlabel}    {parameter_label}    Error \n")

        for i, row in enumerate(parameter_array):
            bin_center = row[0]
            value      = row[1]
            error      = row[2]
            f.write(f"{bin_center: .6e}  {value: .6e}  {error: .6e}\n")


