#!/usr/bin/env python3
"""F + R differential distributions and a check of the fluctuations around
theta_mu ~ 1.4 mrad.

F: devel branch (evtgen_24_11_F_devel, McMule devel, alpha = 1 internally ->
   multiplied by alpha^3 here, see compare_F_branches.py).
R: event-generator branch (evtgen_18_11_NLO0R0, alpha physical inside McMule).

Rules (as in analyze_mcmule_meeting.py): integrated values from Record.value,
histograms added only if all bin centres are identical, no rebinning; errors
added in quadrature (independent runs).

Fluctuation metric: in a region, fit a weighted cubic polynomial to the
histogram and report chi2/ndf of the residuals with respect to the McMule bin
errors (~1 if the scatter is purely statistical and the errors are right) and
the RMS relative residual vs the median relative bin error.  A sliding version
(window of 31 bins) is plotted versus theta.
"""

import argparse
import re
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymule import setup, sigma, mergefks, conv
from pymule.constants import alpha
from scipy.stats import chi2 as chi2_dist

import sys as _sys; _sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, COLORS, PAIR
apply_style()

BASE = Path("/nfs/momos/user/mleibelt")
MOMOS = BASE / "McMule_Results"
TAG = "26_11_meeting"
F_RUN, R_RUN = "evtgen_24_11_F_devel", "evtgen_18_11_NLO0R0"
R_OLD_RUN = "evtgen_18_11_complete"   # older R (different binning), context only
HISTS = ["th3", "Emu", "costh3", "Qsq", "th5", "Eph", "phi5", "x5", "y5"]
PRETTY = {"th3": r"$\theta_\mu$", "Qsq": r"$Q^2$", "Emu": r"$E_\mu$"}
OBS = {"th3": (1e-3, r"Muon angle $\theta_\mu$ [mrad]"), "Qsq": (1e6, r"$Q^2$ [GeV$^2$]")}
COL = {"F": COLORS["F"], "R": COLORS["R"], "F+R": COLORS["SUM"]}
ZOOM = (1.0, 1.8)          # mrad
REGIONS = [(0.6, 0.8), (0.8, 1.0), (1.0, 1.2), (1.2, 1.3), (1.3, 1.5),
           (1.5, 1.6), (1.6, 1.8), (1.8, 2.0), (2.0, 2.4)]  # mrad
WINDOW = 31


def load(run, piece, factor=1.0):
    setup(folder=str(BASE / run / "out"))
    rec = mergefks(sigma(piece)) * (conv * factor)
    return rec


def worker_header(run, piece):
    """Build + userfile cut lines printed by McMule in the worker log."""
    files = sorted((BASE / run).glob(f"worker-{piece}-*.txt"))
    if not files:
        return ["no worker log found"]
    keep, inside = [], False
    for line in files[0].read_text(errors="replace").splitlines()[:20]:
        s = line.strip()
        if re.search(r"SHA|branch", s):
            keep.append(s)
        elif "Welcome" in s:
            inside = True
        elif inside and s.startswith("* ") and not s.startswith("**"):
            keep.append(s)
        elif inside:
            inside = False
    return keep


def h(rec, key):
    return np.asarray(rec[key], dtype=float)[:, :3]


def same_bins(a, b):
    return a.shape == b.shape and np.array_equal(a[:, 0], b[:, 0])


def add(a, b):
    return np.column_stack((a[:, 0], a[:, 1] + b[:, 1], np.hypot(a[:, 2], b[:, 2])))


def region_stats(x, y, e, lo, hi, deg=3):
    m = (x >= lo) & (x < hi) & (e > 0)
    n = int(m.sum())
    if n <= deg + 2:
        return None
    xs, ys, es = x[m], y[m], e[m]
    xc = (xs - xs.mean()) / (xs.ptp() or 1)
    coef = np.polyfit(xc, ys, deg, w=1/es)
    fit = np.polyval(coef, xc)
    res = ys - fit
    chi2 = float(np.sum((res/es)**2))
    ndf = n - deg - 1
    return dict(n=n, chi2ndf=chi2/ndf, pval=float(chi2_dist.sf(chi2, ndf)), rms_rel_res=float(np.sqrt(np.mean((res/fit)**2))),
                med_rel_err=float(np.median(es/np.abs(ys))), maxpull=float(np.max(np.abs(res/es))),
                resid=(xs, res/es))


def sliding(x, y, e, w=WINDOW, deg=3):
    half = w // 2
    xc, c2 = [], []
    for i in range(half, len(x) - half):
        s = slice(i - half, i + half + 1)
        if np.any(e[s] <= 0):
            continue
        xx = (x[s] - x[s].mean()) / x[s].ptp()
        coef = np.polyfit(xx, y[s], deg, w=1/e[s])
        r = (y[s] - np.polyval(coef, xx)) / e[s]
        xc.append(x[i]); c2.append(np.sum(r**2) / (w - deg - 1))
    return np.array(xc), np.array(c2)


def draw_dist(path, hs, key, zoom=None, title=""):
    scale, xlabel = OBS[key]
    unit = xlabel.split("[")[-1].rstrip("]")
    nrow = 2 if zoom else 1
    fig, axes = plt.subplots(nrow, 1, figsize=(10, 6.0 * nrow))
    axes = np.atleast_1d(axes)
    for ax, rng in zip(axes, [None, zoom][:nrow]):
        for lab, hh in hs.items():
            x = hh[1:-1, 0] / scale
            y, e = hh[1:-1, 1] * scale, hh[1:-1, 2] * scale
            m = slice(None) if rng is None else (x >= rng[0]) & (x <= rng[1])
            if rng is None:
                ax.step(x[m], y[m], where="mid", color=COL[lab], lw=1.3, label=lab)
                ax.fill_between(x[m], (y-e)[m], (y+e)[m], step="mid", color=COL[lab], alpha=.22)
            else:
                ax.errorbar(x[m], y[m], yerr=e[m], fmt="o", ms=2.5, lw=.8, color=COL[lab],
                            label=lab, capsize=0)
        ax.axhline(0, color=COLORS["REF"], lw=.7)
        ax.grid(alpha=.2)
        ax.set(xlabel=xlabel, ylabel=r"$d\sigma/dx$ [$\mu$b / " + unit + "]")
        if rng is None:
            ax.set_title(title)
            if key == "th3":
                ax.axvspan(*ZOOM, color=COLORS["HIGHLIGHT"], alpha=.18, lw=0)
        else:
            ax.axvspan(1.3, 1.5, color=COLORS["HIGHLIGHT"], alpha=.35, lw=0, label="1.3-1.5 mrad")
            ax.set_title(f"Zoom {rng[0]}-{rng[1]} mrad (points with McMule bin errors)")
        ax.legend(frameon=False)
    if key == "Qsq":
        axes[0].set_xscale("log")
        axes[0].set_yscale("symlog", linthresh=1e-3)
    if key == "th3":
        axes[0].set_xscale("log")
        axes[0].set_yscale("symlog", linthresh=1e-3)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description="F + R plots and 1.4 mrad fluctuation check")
    p.add_argument("--outdir", default=str(MOMOS / "Figures" / TAG))
    p.add_argument("--outdir-vals", default=str(MOMOS / "Vals" / TAG))
    args = p.parse_args()
    outfig, outval = Path(args.outdir), Path(args.outdir_vals)
    outfig.mkdir(parents=True, exist_ok=True)
    outval.mkdir(parents=True, exist_ok=True)

    F = load(F_RUN, "mp2mpF", alpha**3)
    R = load(R_RUN, "mp2mpR")
    Rold = load(R_OLD_RUN, "mp2mpR")
    L = []
    out = L.append
    out("TASK 2: F (devel, x alpha^3) + R (event-generator) — units microbarn")
    out(f"F: {F_RUN}/out mp2mpF; R: {R_RUN}/out mp2mpR (xi_c = 1.0 for both)")
    out("")
    out("COMPATIBILITY")
    for run, piece in ((F_RUN, "mp2mpF"), (R_RUN, "mp2mpR")):
        out(f"  {run} {piece}: " + " | ".join(worker_header(run, piece)))
    diff = [k for k in HISTS if not same_bins(h(F, k), h(R, k))]
    out("  histogram bin centres F vs R: " + ("identical for all 9 histograms" if not diff
                                              else "DIFFERENT: " + ", ".join(diff)))
    out("  NOTE: user.so differs (compiled against different McMule builds); the printed cuts agree.")
    out("")
    vF, eF = map(float, F.value[:2]); vR, eR = map(float, R.value[:2])
    out("INTEGRATED (Record.value)")
    out(f"  F      = {vF:.10e} +- {eF:.3e} ub")
    out(f"  R      = {vR:.10e} +- {eR:.3e} ub")
    out(f"  F + R  = {vF+vR:.10e} +- {np.hypot(eF, eR):.3e} ub")
    out("")

    hs = {}
    for key in ("th3", "Qsq"):
        hf, hr = h(F, key), h(R, key)
        if not same_bins(hf, hr):
            out(f"  {key}: F and R binning differ -> F+R NOT built")
            continue
        hs[key] = {"F": hf, "R": hr, "F+R": add(hf, hr)}
        path = outfig / f"{TAG}_F_R_FplusR_{key}.png"
        draw_dist(path, hs[key], key, zoom=ZOOM if key == "th3" else None,
                  title=f"F (devel), R (event-gen.) and F + R, {PRETTY.get(key, key)}")
        out(f"  saved {path}")
    out("")

    # ---------------- fluctuation metric on th3 ----------------
    sc = 1e-3
    series = dict(hs["th3"])
    series["R_old(18_11_complete, old binning)"] = h(Rold, "th3")
    out(f"FLUCTUATIONS in th3: weighted cubic fit per region; chi2/ndf of residuals w.r.t. bin errors")
    out(f"  {'series':36s} {'region[mrad]':13s} {'nbins':>5s} {'chi2/ndf':>8s} {'p-value':>8s} {'max|pull|':>9s} "
        f"{'rms rel resid':>13s} {'median rel err':>14s}")
    table = []
    for lab, hh in series.items():
        x, y, e = hh[1:-1, 0]/sc, hh[1:-1, 1]*sc, hh[1:-1, 2]*sc
        for lo, hi in REGIONS:
            s = region_stats(x, y, e, lo, hi)
            if s is None:
                continue
            table.append((lab, lo, hi, s))
            flag = "  <-- 1.4 mrad region" if (lo, hi) == (1.3, 1.5) else ""
            out(f"  {lab:36s} {lo:4.1f}-{hi:<4.1f}     {s['n']:5d} {s['chi2ndf']:8.2f} {s['pval']:8.3f} "
                f"{s['maxpull']:9.2f} {s['rms_rel_res']:13.3e} {s['med_rel_err']:14.3e}{flag}")
    out("")
    out("  Interpretation: chi2/ndf ~ 1 (p not small) -> scatter consistent with the quoted bin errors;")
    out("  a localised structure at 1.4 mrad shows up as a small p-value in 1.3-1.5 compared to the")
    out("  neighbouring regions. 'rms rel resid' is unstable where the fit crosses zero (R near 1.5-1.6).")
    out("  Regions are not a trials-corrected test: with ~9 regions per series one p<0.1 is expected.")
    out("")

    # bins whose error is far above the local error level (single large vegas
    # weights / FKS event-counterevent pairs split over neighbouring bins)
    from scipy.ndimage import median_filter
    out("ERROR-OUTLIER BINS in th3 (bin error > 4 x running median of errors over 15 bins)")
    others = {"R_24_11_NLO0corrected (context)": h(load("evtgen_24_11_NLO0corrected", "mp2mpR"), "th3")}
    for lab, hh in {**series, **others}.items():
        x, y, e = hh[1:-1, 0]/sc, hh[1:-1, 1]*sc, hh[1:-1, 2]*sc
        ratio = e / np.maximum(median_filter(e, size=15), 1e-300)
        idx = np.where(ratio > 4)[0]
        txt_b = ", ".join(f"{x[i]:.3f} mrad: {y[i]:+.4f}+-{e[i]:.4f} ({ratio[i]:.0f}x)" for i in idx)
        out(f"  {lab}: {txt_b if len(idx) else 'none'}")

    # sliding chi2/ndf plot + residual zoom
    fig, ax = plt.subplots(2, 1, figsize=(10, 9.5))
    for lab, hh in series.items():
        x, y, e = hh[1:-1, 0]/sc, hh[1:-1, 1]*sc, hh[1:-1, 2]*sc
        xc, c2 = sliding(x, y, e)
        c = COL.get(lab, COLORS["DIFF"])
        ls = "-" if lab in COL else ":"
        ax[0].plot(xc, c2, ls, color=c, lw=1.2, label=lab)
    ax[0].axhline(1, color=COLORS["REF"], lw=.7)
    ax[0].axvspan(1.3, 1.5, color=COLORS["HIGHLIGHT"], alpha=.35, lw=0)
    ax[0].set(xlabel=r"Muon angle $\theta_\mu$ [mrad]", ylabel=f"chi2/ndf (cubic fit, {WINDOW}-bin window)",
              title="Local smoothness test: scatter vs McMule bin errors", yscale="log")
    ax[0].legend(frameon=False)
    ax[0].grid(alpha=.2)
    hh = series["F+R"]
    x, y, e = hh[1:-1, 0]/sc, hh[1:-1, 1]*sc, hh[1:-1, 2]*sc
    s = region_stats(x, y, e, *ZOOM)
    ax[1].plot(*s["resid"], "o", ms=3, color=COL["F+R"], label="F+R: (data - cubic fit)/error")
    for k in (-2, 0, 2):
        ax[1].axhline(k, color=COLORS["REF"], lw=.7, ls="-" if k == 0 else "--")
    ax[1].axvspan(1.3, 1.5, color=COLORS["HIGHLIGHT"], alpha=.35, lw=0)
    ax[1].set(xlabel=r"Muon angle $\theta_\mu$ [mrad]", ylabel=r"residual / $\sigma$",
              title=f"F+R residuals {ZOOM[0]}-{ZOOM[1]} mrad, chi2/ndf = {s['chi2ndf']:.2f}")
    ax[1].legend(frameon=False)
    ax[1].grid(alpha=.2)
    fig.tight_layout()
    path = outfig / f"{TAG}_FplusR_th3_fluctuation_check.png"
    fig.savefig(path, dpi=170)
    plt.close(fig)
    out(f"  saved {path}")

    txt = outval / f"{TAG}_F_plus_R_report.txt"
    txt.write_text("\n".join(L) + "\n", encoding="utf-8")
    tsv = outval / f"{TAG}_F_plus_R_th3_regions.tsv"
    with tsv.open("w", encoding="utf-8") as f:
        f.write("series\tlo_mrad\thi_mrad\tnbins\tchi2ndf\tp_value\tmax_abs_pull\trms_rel_resid\tmedian_rel_err\n")
        for lab, lo, hi, s in table:
            f.write(f"{lab}\t{lo}\t{hi}\t{s['n']}\t{s['chi2ndf']:.4f}\t{s['pval']:.4f}\t{s['maxpull']:.3f}\t"
                    f"{s['rms_rel_res']:.4e}\t{s['med_rel_err']:.4e}\n")
    for key, d in hs.items():
        p = outval / f"{TAG}_F_R_FplusR_{key}_hist.tsv"
        with p.open("w", encoding="utf-8") as f:
            f.write(f"# {key}: raw McMule x (centre) and dsigma/dx in ub per raw unit; rows 0/-1 = under/overflow\n")
            f.write("x\tF\tF_err\tR\tR_err\tFR\tFR_err\n")
            for i in range(len(d["F"])):
                f.write("\t".join(f"{v:.10e}" for v in (d["F"][i, 0], d["F"][i, 1], d["F"][i, 2],
                                                         d["R"][i, 1], d["R"][i, 2],
                                                         d["F+R"][i, 1], d["F+R"][i, 2])) + "\n")
        L.append(f"  saved {p}")
    print("\n".join(L))
    print(f"Saved {txt}\nSaved {tsv}")


if __name__ == "__main__":
    main()
