"""Plotting helpers for McMule mu-p LO and NLO predictions.

Conventions:
    LO        = mp2mp0
    NLO       = mp2mpNLO0 + mp2mpR
    DeltaNLO  = NLO - LO
    K         = DeltaNLO / NLO   (user-chosen, non-standard)
"""

import matplotlib
matplotlib.use("Agg")

import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
import numpy as np
from pathlib import Path

from pymule import errorband, mergebins, scaleplot


DEFAULT_COLORS = {
    "lo": "#1f77b4",
    "nlo": "#2ca02c",
    "delta": "#ff7f0e",
    "k": "#d62728",
    "real": "#9467bd",
    "analytic": "#1f1f1f",
}

plt.rcParams.update({
    "font.family": "serif",
    "font.size": 12,
    "text.usetex": False,
    "axes.grid": True,
    "grid.alpha": 0.30,
    "grid.linestyle": ":",
    "savefig.dpi": 300,
})


def finite_bins(hist):
    arr = np.asarray(hist, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError("Expected histogram with columns [x, value, error].")
    return arr[np.all(np.isfinite(arr[:, :3]), axis=1), :3]


def write_file_with_values(filename, parameter_array, xlabel, parameter_label):
    path = Path(filename)
    path.parent.mkdir(parents=True, exist_ok=True)
    arr = finite_bins(parameter_array)
    np.savetxt(path, arr, fmt="%.12e",
               header=f"{xlabel}    {parameter_label}    MC_stat_error")


def calculate_rate(sigma_mb):
    ltarget = 80.0
    Hpres = 20.0
    Npvol = 2.687e19
    Ibeam = 2e6
    Np = Npvol * 2.0 * (Hpres / 1.013) * ltarget
    return float(sigma_mb) * 1e-27 * Np * Ibeam


def _save(fig, name, outdir):
    folder = Path(outdir)
    folder.mkdir(parents=True, exist_ok=True)
    fig.savefig(folder / f"{name}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def _draw_mc(ax, hist, color, label=None, zorder=3, linestyle="-"):
    arr = finite_bins(hist)
    npoints = len(arr)
    if npoints == 0:
        return
    if npoints < 3:
        ax.errorbar(
            arr[:, 0], arr[:, 1], yerr=arr[:, 2],
            fmt="o" if npoints == 1 else "o-",
            color=color, markersize=4.5, linewidth=1.6,
            elinewidth=1.2, capsize=3, zorder=zorder,
        )
    else:
        plt.sca(ax)
        artists = errorband(arr)
        for artist in artists:
            artist.set_color(color)
            artist.set_alpha(0.70 if isinstance(artist, Line2D) else 0.18)
            artist.set_zorder(zorder)
            if isinstance(artist, Line2D):
                artist.set_linestyle(linestyle)
                artist.set_linewidth(1.8)
    if label:
        ax.plot([], [], color=color, label=label, linewidth=2.0, linestyle=linestyle)


def _scale(hist, scale_factor):
    return finite_bins(scaleplot(hist, scale_factor))


def _estimate_edges(x):
    x = np.asarray(x, dtype=float)
    if x.ndim != 1 or len(x) < 2:
        raise ValueError("Need at least two bin centers to infer edges.")
    mids = 0.5 * (x[:-1] + x[1:])
    left = x[0] - 0.5 * (x[1] - x[0])
    right = x[-1] + 0.5 * (x[-1] - x[-2])
    return np.concatenate(([left], mids, [right]))


def rebin_variable(hist, new_edges):
    """Rebin a differential histogram to arbitrary new bin edges.

    Hist must have columns [x_center, dσ/dx, err(dσ/dx)].
    Errors are added in quadrature on the integrated bin content.
    """
    arr = finite_bins(hist)
    x, y, ey = arr.T
    old_edges = _estimate_edges(x)
    widths = np.diff(old_edges)

    new_edges = np.asarray(new_edges, dtype=float)
    if np.any(np.diff(new_edges) <= 0):
        raise ValueError("new_edges must be strictly increasing.")

    out = []
    idx = np.digitize(x, new_edges) - 1
    for i in range(len(new_edges) - 1):
        mask = idx == i
        if not np.any(mask):
            continue
        width_new = new_edges[i + 1] - new_edges[i]
        content = np.sum(y[mask] * widths[mask])
        err2 = np.sum((ey[mask] * widths[mask]) ** 2)
        out.append([
            0.5 * (new_edges[i] + new_edges[i + 1]),
            content / width_new,
            np.sqrt(err2) / width_new,
        ])
    if not out:
        return np.zeros((0, 3))
    return np.asarray(out, dtype=float)


def _ratio_delta_over_nlo_arrays(lo, nlo, k_xmin=None, k_xmax=None):
    x, l, sl = lo.T
    _, n, sn = nlo.T
    valid = np.isfinite(n) & np.isfinite(sn) & (np.abs(n) > 0)
    valid &= np.abs(n) > 2.0 * sn
    if k_xmin is not None:
        valid &= x >= k_xmin
    if k_xmax is not None:
        valid &= x <= k_xmax
    with np.errstate(divide="ignore", invalid="ignore"):
        k = 1.0 - l / n
        sk = np.sqrt((sl / n) ** 2 + ((l * sn) / (n * n)) ** 2)
    valid &= np.isfinite(k) & np.isfinite(sk)
    return np.column_stack((x[valid], k[valid], sk[valid]))


def _ratio_delta_over_nlo(lo_hist, nlo_hist, scale_factor, merge=5,
                          k_xmin=None, k_xmax=None):
    lo = _scale(lo_hist, scale_factor)
    nlo = _scale(nlo_hist, scale_factor)
    if merge > 1:
        lo = finite_bins(mergebins(lo, merge))
        nlo = finite_bins(mergebins(nlo, merge))
    return _ratio_delta_over_nlo_arrays(lo, nlo, k_xmin=k_xmin, k_xmax=k_xmax)


def save_lo_nlo_comparison(*, savename, lo_hist, nlo_hist, delta_hist,
                           scale_factor, x_label, y_label, main_title,
                           outdir, xlim=None, force_main_linear=False,
                           k_xmin=None, k_xmax=None, k_ylim=None, k_merge=5,
                           show_k=True, rebin_edges=None, note=None):
    if show_k:
        fig, (ax, ax_k) = plt.subplots(
            2, 1, figsize=(7.5, 6.0), sharex=True,
            gridspec_kw={"height_ratios": (3, 1.15), "hspace": 0.05},
        )
    else:
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax_k = None

    lo = _scale(lo_hist, scale_factor)
    nlo = _scale(nlo_hist, scale_factor)
    delta = _scale(delta_hist, scale_factor)
    if (lo.shape != nlo.shape or delta.shape != lo.shape
            or not np.allclose(lo[:, 1] + delta[:, 1], nlo[:, 1],
                               rtol=1e-7, atol=1e-12)):
        raise ValueError("Expected DeltaNLO = NLO - LO in every bin.")

    if rebin_edges is not None:
        lo = rebin_variable(lo, rebin_edges)
        nlo = rebin_variable(nlo, rebin_edges)
        delta = rebin_variable(delta, rebin_edges)

    _draw_mc(ax, nlo, DEFAULT_COLORS["nlo"], "NLO", zorder=3)
    _draw_mc(ax, lo, DEFAULT_COLORS["lo"], "LO", zorder=4, linestyle="--")
    ax.set_title(main_title)
    ax.set_ylabel(y_label)
    ax.legend(frameon=False)
    if note:
        ax.text(0.02, 0.02, note, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=9)

    positive = np.any(lo[:, 1] > 0) or np.any(nlo[:, 1] > 0)
    if not force_main_linear and positive:
        ax.set_yscale("log")
    if show_k:
        ax.tick_params(axis="x", labelbottom=False)
    else:
        ax.set_xlabel(x_label)

    if show_k:
        if rebin_edges is not None:
            k = _ratio_delta_over_nlo_arrays(lo, nlo, k_xmin=k_xmin, k_xmax=k_xmax)
        else:
            k = _ratio_delta_over_nlo(lo_hist, nlo_hist, scale_factor,
                                      merge=k_merge, k_xmin=k_xmin, k_xmax=k_xmax)
        if len(k):
            _draw_mc(ax_k, k, DEFAULT_COLORS["k"], zorder=4)
        if len(k) < 3:
            print(f"NOTE: {savename}: only {len(k)} K bins survived cuts; shown as individual point(s).")
        ax_k.axhline(0, color="0.35", linewidth=0.9, linestyle="--")
        ax_k.set_ylabel(r"$K = \Delta\mathrm{NLO}/\mathrm{NLO}$")
        ax_k.set_xlabel(x_label)
        if k_ylim is not None:
            ax_k.set_ylim(*k_ylim)
    else:
        k = None
    if xlim is not None:
        ax.set_xlim(*xlim)
        if ax_k is not None:
            ax_k.set_xlim(*xlim)

    _save(fig, savename, outdir)
    return k, lo, nlo, delta


def save_delta_nlo(*, savename, delta_hist, scale_factor, x_label, y_label,
                   title, outdir, xlim=None, ymin=None,
                   rebin_edges=None, note=None):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    h = _scale(delta_hist, scale_factor)
    if rebin_edges is not None:
        h = rebin_variable(h, rebin_edges)
    _draw_mc(ax, h, DEFAULT_COLORS["delta"], r"$\Delta\mathrm{NLO}$")
    ax.axhline(0, linewidth=0.9, linestyle="--", color="0.35")
    nonzero = np.abs(h[:, 1])
    nonzero = nonzero[nonzero > 0]
    if len(nonzero):
        ax.set_yscale("symlog", linthresh=max(float(np.max(nonzero)) * 0.002, 1e-14))
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    if note:
        ax.text(0.02, 0.02, note, transform=ax.transAxes,
                ha="left", va="bottom", fontsize=9)
    if xlim is not None:
        ax.set_xlim(*xlim)
    if ymin is not None:
        ax.set_ylim(bottom=ymin)
    ax.legend(frameon=False)
    _save(fig, savename, outdir)
    return h


def save_real_photon(*, savename, real_hist, scale_factor, x_label, y_label,
                     title, outdir, xlim=None):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    h = _scale(real_hist, scale_factor)
    _draw_mc(ax, h, DEFAULT_COLORS["real"], "R (FKS-subtracted)")
    positive = h[:, 1] > 0
    negative = h[:, 1] < 0
    if np.any(negative):
        absvals = np.abs(h[:, 1]); absvals = absvals[absvals > 0]
        ax.set_yscale("symlog", linthresh=max(float(np.max(absvals)) * 0.002, 1e-14))
        ax.axhline(0, linewidth=0.8, color="0.35", linestyle="--")
    elif np.any(positive):
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.text(0.02, 0.02, r"FKS-subtracted $R$ piece ($\xi$ dependent)",
            transform=ax.transAxes, ha="left", va="bottom", fontsize=9)
    if xlim is not None:
        ax.set_xlim(*xlim)
    _save(fig, savename, outdir)


def save_analytic_lo_check(*, savename, lo_hist, scale_factor, func,
                           x_label, y_label, title, outdir,
                           xlim=None, yscale="log"):
    h = _scale(lo_hist, scale_factor)
    x = h[:, 0]
    mask = np.isfinite(x)
    x = x[mask]
    with np.errstate(divide="ignore", invalid="ignore", over="ignore"):
        analytic = np.asarray([float(func(float(q))) for q in x], dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    _draw_mc(ax, h, DEFAULT_COLORS["lo"], "LO (numerical)")
    good = np.isfinite(analytic)
    if yscale == "log":
        good &= analytic > 0
    ax.plot(x[good], analytic[good], color=DEFAULT_COLORS["analytic"],
            linestyle="--", linewidth=1.4, label="LO analytic (at bin center)")
    ax.set_title(title)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_yscale(yscale)
    if xlim is not None:
        ax.set_xlim(*xlim)
    ax.legend(frameon=False)
    _save(fig, savename, outdir)
    return np.column_stack((x, analytic))


# ====================================================================
# Variable-width rebinning: original-bin-aligned grouping, no re-run.
# This is NOT a smoothing/interpolation procedure: integrated weights
# and their MC errors are summed over contiguous existing histogram bins.
# ====================================================================

def adaptive_grouped_triplet(lo_hist, nlo_hist, delta_hist, scale_factor,
                             regions, xmax=None):
    """Rebin LO, NLO and Delta with identical, exact original-bin edges.

    regions: ordered sequence (upper_x, bins_to_merge), in displayed units.
    Example for theta_mu: [(0.6,2), (1.4,3), (2.5,5), (np.inf,8)].
    McMule differential values are integrated over each original bin first.
    The new bin's differential cross section equals summed content / width.
    Statistical errors assume independent original histogram bins.

    Returns (lo, nlo, delta, edges), all in displayed x and y units.
    """
    lo = _scale(lo_hist, scale_factor)
    nlo = _scale(nlo_hist, scale_factor)
    delta = _scale(delta_hist, scale_factor)
    if len(lo) < 2 or lo.shape != nlo.shape or lo.shape != delta.shape:
        raise ValueError("Adaptive rebinning needs at least two matching bins.")
    x = lo[:, 0]
    if np.any(np.diff(x) <= 0):
        raise ValueError("Adaptive rebinning expects sorted, unique bin centers.")
    if not (np.allclose(x, nlo[:, 0], rtol=1e-9, atol=1e-12) and
            np.allclose(x, delta[:, 0], rtol=1e-9, atol=1e-12)):
        raise ValueError("LO/NLO/Delta x binning differs.")
    if not np.allclose(lo[:, 1] + delta[:, 1], nlo[:, 1], rtol=1e-6, atol=1e-12):
        raise ValueError("Expected DeltaNLO = NLO - LO before rebinning.")

    old_edges = _estimate_edges(x)
    dx = np.diff(old_edges)
    stop = len(x) if xmax is None else int(np.searchsorted(x, xmax, side="right"))
    if stop < 2:
        raise ValueError("Adaptive rebinning range contains fewer than two bins.")
    groups = []
    start = 0
    while start < stop:
        for bound, n_merge in regions:
            if x[start] <= bound:
                region_end = min(stop, int(np.searchsorted(x, bound, side="right")))
                end = min(start + int(n_merge), region_end)
                break
        else:
            raise ValueError("Adaptive regions must cover the entire x range.")
        if end <= start:
            raise RuntimeError("Invalid adaptive group boundary.")
        groups.append((start, end))
        start = end

    edges = np.asarray([old_edges[groups[0][0]]] +
                       [old_edges[j] for _, j in groups])

    def merge_grouped(a):
        rows = []
        for i, j in groups:
            new_width = float(old_edges[j] - old_edges[i])
            content = float(np.sum(a[i:j, 1] * dx[i:j]))
            error = float(np.sqrt(np.sum((a[i:j, 2] * dx[i:j]) ** 2)))
            rows.append([0.5 * (old_edges[i] + old_edges[j]),
                         content / new_width, error / new_width])
        return np.asarray(rows)

    lo_b, nlo_b, delta_b = map(merge_grouped, (lo, nlo, delta))
    if not np.allclose(lo_b[:, 1] + delta_b[:, 1], nlo_b[:, 1],
                       rtol=1e-6, atol=1e-12):
        raise RuntimeError("Adaptive rebinning violated DeltaNLO = NLO - LO.")

    # Verify that integrating the rebinned bins reproduces the selected input.
    for src, dst in ((lo, lo_b), (nlo, nlo_b), (delta, delta_b)):
        before = np.sum(src[:stop, 1] * dx[:stop])
        after = np.sum(dst[:, 1] * np.diff(edges))
        if not np.isclose(before, after, rtol=1e-11, atol=1e-10 * max(1.0, abs(before))):
            raise RuntimeError("Integrated cross section changed during rebinning.")
    return lo_b, nlo_b, delta_b, edges


def _draw_variable_steps(ax, h, edges, color, label, linestyle="-", zorder=3):
    """Show actual variable-width bins and a 1-sigma MC error band."""
    ax.stairs(h[:, 1], edges, label=label, color=color, linewidth=1.8,
              linestyle=linestyle, zorder=zorder)
    low = h[:, 1] - h[:, 2]
    high = h[:, 1] + h[:, 2]
    ax.fill_between(edges, np.r_[low, low[-1]], np.r_[high, high[-1]],
                    step="post", alpha=0.14, color=color, linewidth=0,
                    zorder=zorder-1)


def save_adaptive_pair(*, savename, lo, nlo, delta, edges,
                       x_label, y_label, title, outdir, xlim=None,
                       show_k=True, k_xmax=None, note="variable-width bins"):
    if show_k:
        fig, (ax, ax_k) = plt.subplots(2, 1, figsize=(7.5, 6.0), sharex=True,
                                      gridspec_kw={"height_ratios": (3, 1.15), "hspace": 0.05})
    else:
        fig, ax = plt.subplots(figsize=(7.5, 4.5))
        ax_k = None
    _draw_variable_steps(ax, nlo, edges, DEFAULT_COLORS["nlo"], "NLO")
    _draw_variable_steps(ax, lo, edges, DEFAULT_COLORS["lo"], "LO", "--", 4)
    ax.set_title(title)
    ax.set_ylabel(y_label)
    ax.legend(frameon=False)
    if np.any(lo[:, 1] > 0) or np.any(nlo[:, 1] > 0):
        ax.set_yscale("log")
    ax.text(0.02, 0.02, note, fontsize=9, va="bottom", transform=ax.transAxes)
    if not show_k:
        ax.set_xlabel(x_label)
    else:
        k = _ratio_delta_over_nlo_arrays(lo, nlo, k_xmax=k_xmax)
        if len(k):
            # Draw individual valid K-bins, never connect across invalid bins.
            ax_k.errorbar(k[:, 0], k[:, 1], yerr=k[:, 2], fmt="o", markersize=3,
                          color=DEFAULT_COLORS["k"], linewidth=1, capsize=2)
        ax_k.axhline(0, color="0.35", linewidth=0.9, linestyle="--")
        ax_k.set_ylabel(r"$K = \Delta\mathrm{NLO}/\mathrm{NLO}$")
        ax_k.set_xlabel(x_label)
    if xlim is not None:
        ax.set_xlim(*xlim)
    _save(fig, savename, outdir)
    return k if show_k else None


def save_adaptive_delta(*, savename, delta, edges, x_label, y_label,
                        title, outdir, xlim=None, note="variable-width bins"):
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    _draw_variable_steps(ax, delta, edges, DEFAULT_COLORS["delta"],
                         r"$\Delta\mathrm{NLO}$")
    # Signed corrections require a signed axis. If all values are positive,
    # a log plot is clearer and avoids artificial negative tick labels.
    y, dy = delta[:, 1], delta[:, 2]
    if np.any(y <= 0) or np.any(y - dy <= 0):
        amp = max(float(np.max(np.abs(y))), float(np.max(dy)), 1e-12)
        ax.set_yscale("symlog", linthresh=0.03 * amp)
        ax.axhline(0, color="0.35", linewidth=0.9, linestyle="--")
    else:
        ax.set_yscale("log")
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_title(title)
    ax.legend(frameon=False)
    ax.text(0.02, 0.02, note, fontsize=9, va="bottom", transform=ax.transAxes)
    if xlim is not None:
        ax.set_xlim(*xlim)
    _save(fig, savename, outdir)
