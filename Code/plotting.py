import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path

def style_axis(
    ax,
    xlabel=None,
    ylabel=None,
    title=None,       # <-- neu hinzufügen
    xscale="linear",
    yscale="linear",
    legend=True
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



def create_figure(
    nrows=1,
    ncols=1,
    figsize=(8, 5),
    font_size=7
):
    plt.rcParams.update({
        "font.size": font_size,
        "font.family": "serif"
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize)

    # Wichtig: axes immer als Array behandeln
    axes = np.atleast_1d(axes)

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
