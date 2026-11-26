#!/usr/bin/env python3
"""Compare McMule F from the event-generator branch and from the devel branch,
and check NLO0 - LO against F (and F + A).

Conventions (same as analyze_mcmule_meeting.py):
  * integrated cross sections come from pymule's per-piece Record.value
    (vegas integral), NEVER from histogram sums;
  * histograms are compared bin-by-bin ONLY if all bin centres are identical;
    nothing is rebinned or interpolated;
  * errors of differences/sums are added in quadrature (independent seeds,
    covariance ignored).

Units: mergefks(...) * conv gives microbarn (raw McMule unit MeV^-2).

Coupling: the event-generator branch sets alpha = 1/137.035999084 inside
McMule (mcmule_header.f95), the devel branch keeps alpha = 1 (pre/global_def.f95)
as in standard McMule. F is O(alpha^3) (Born alpha^2 x one-loop alpha), so the
devel F is multiplied by pymule.constants.alpha**3 (= 1/137.035999084**3).

Outputs go to McMule_Results/{Figures,Vals}/26_11_meeting.
"""

import argparse
import math
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from pymule import setup, sigma, mergefks, conv
from pymule.constants import alpha

import sys as _sys; _sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, COLORS, PAIR
apply_style()

BASE = Path("/nfs/momos/user/mleibelt")
MOMOS = BASE / "McMule_Results"
TAG = "26_11_meeting"

# (label, run folder, McMule piece, extra factor)
SOURCES = [
    ("LO",       "evtgen_18_11_complete",       "mp2mp0",    1.0),
    ("NLO0",     "evtgen_18_11_complete",       "mp2mpNLO0", 1.0),
    ("F_evtgen", "evtgen_18_11_complete",       "mp2mpF",    1.0),
    ("A",        "evtgen_18_11_complete",       "mp2mpA",    1.0),
    ("R",        "evtgen_18_11_complete",       "mp2mpR",    1.0),
    ("F_devel",  "evtgen_24_11_F_devel",        "mp2mpF",    alpha**3),
    # independent replica of LO / NLO0 (same build + userfile as F_evtgen, A)
    ("LO_rep",   "evtgen_18_11_NLO0R0",         "mp2mp0",    1.0),
    ("NLO0_rep", "evtgen_18_11_NLO0R0",         "mp2mpNLO0", 1.0),
    # new event-generator build (cda91c8): NLO0 now includes A
    ("NLO0_newA", "evtgen_24_11_NLO0corrected", "mp2mpNLO0", 1.0),
]

HISTS = ["th3", "Emu", "costh3", "Qsq", "th5", "Eph", "phi5", "x5", "y5"]
PRETTY = {"th3": r"$\theta_\mu$", "Qsq": r"$Q^2$", "Emu": r"$E_\mu$"}
OBS = {  # key: (scale, x label)  raw/scale -> display, density*scale
    "th3": (1e-3, r"Muon angle $\theta_\mu$ [mrad]"),
    "Qsq": (1e6, r"$Q^2$ [GeV$^2$]"),
    "Emu": (1e3, "Muon energy [GeV]"),
}


def load(run, piece, factor):
    setup(folder=str(BASE / run / "out"))
    xs = sigma(piece)
    meta = []
    for key, rec in xs._children.items():
        meta.append(f"seed={key[0]} xi={key[1]} file-format={type(rec).__name__} "
                    f"sha={getattr(rec, 'sha', '?')} chi2a={getattr(rec, '_chi2a', '?')}")
    rec = mergefks(xs) * (conv * factor)
    return rec, "; ".join(meta)


def val(rec):
    v, e = map(float, rec.value[:2])
    return v, abs(e)


def comb(terms):
    """terms: list of (sign, (v, e)) -> (value, quadrature error)."""
    v = sum(s * x[0] for s, x in terms)
    e = math.sqrt(sum(x[1] ** 2 for _, x in terms))
    return v, e


def hist(rec, key):
    return np.asarray(rec[key], dtype=float)[:, :3]


def same_bins(h1, h2):
    return h1.shape == h2.shape and np.array_equal(h1[:, 0], h2[:, 0])


def perbin_stats(h1, h2):
    """Pulls (h1-h2)/sqrt(e1^2+e2^2) over bins where the combined error > 0."""
    d = h1[:, 1] - h2[:, 1]
    e = np.hypot(h1[:, 2], h2[:, 2])
    m = e > 0
    pulls = d[m] / e[m]
    with np.errstate(divide="ignore", invalid="ignore"):
        rel = np.where(h2[:, 1] != 0, d / h2[:, 1], np.nan)
    return pulls, rel, m


def plot_compare(path, ha, hb, la, lb, key, title, sa='A', sb='B'):
    """sa, sb: short names used in the axis labels of the ratio and pull panels."""
    scale, xlabel = OBS[key]
    inner = slice(1, -1)  # drop under/overflow (centres -inf/+inf)
    x = ha[inner, 0] / scale
    ya, ea = ha[inner, 1] * scale, ha[inner, 2] * scale
    yb, eb = hb[inner, 1] * scale, hb[inner, 2] * scale
    e = np.hypot(ea, eb)
    with np.errstate(divide="ignore", invalid="ignore"):
        pull = np.where(e > 0, (ya - yb) / e, np.nan)
        rel = np.where(yb != 0, 100 * (ya - yb) / yb, np.nan)        # in percent
        rel_e = np.where(yb != 0, 100 * e / np.abs(yb), np.nan)
    fig, ax = plt.subplots(3, 1, figsize=(10, 10), sharex=True,
                           gridspec_kw={"height_ratios": [2.2, 1, 1]})
    for y, err, lab, c in ((ya, ea, la, PAIR[0]), (yb, eb, lb, PAIR[1])):
        ax[0].step(x, y, where="mid", color=c, lw=1.4, label=lab)
        ax[0].fill_between(x, y - err, y + err, step="mid", color=c, alpha=.2)
    unit = xlabel.split("[")[-1].rstrip("]")
    ax[0].set(ylabel=r"$d\sigma/dx$ [$\mu$b / " + unit + "]", title=title)
    ax[0].legend(frameon=False)
    ax[1].step(x, rel, where="mid", color=COLORS["DIFF"], lw=1.5)
    ax[1].fill_between(x, -rel_e, rel_e, step="mid", color=COLORS["BAND"], alpha=.45,
                       label=r"$\pm 1\sigma$ (combined)")
    ax[1].axhline(0, color=COLORS["REF"], lw=.7)
    ax[1].set(ylabel="relative difference [%]")
    ax[1].text(0.01, 0.95, f"$100\\cdot({sa} - {sb})\\,/\\,{sb}$", transform=ax[1].transAxes, va="top", fontsize=13)
    ax[1].legend(frameon=False)
    ax[2].plot(x, pull, ".", ms=5, color=COLORS["PULL"])
    for k in (-2, 0, 2):
        ax[2].axhline(k, color=COLORS["REF"], lw=.7, ls="-" if k == 0 else "--")
    ax[2].set(ylabel=r"difference / $\sigma$", xlabel=xlabel)
    ax[2].text(0.01, 0.95, f"$({sa} - {sb})\\,/\\,\\sigma$,   $\\sigma$ = combined MC error per bin",
               transform=ax[2].transAxes, va="top", fontsize=13)
    for a in ax:
        a.grid(alpha=.2)
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def main():
    p = argparse.ArgumentParser(description=__doc__.split("\n")[0])
    p.add_argument("--outdir", default=str(MOMOS / "Figures" / TAG))
    p.add_argument("--outdir-vals", default=str(MOMOS / "Vals" / TAG))
    args = p.parse_args()
    outfig, outval = Path(args.outdir), Path(args.outdir_vals)
    outfig.mkdir(parents=True, exist_ok=True)
    outval.mkdir(parents=True, exist_ok=True)

    rec, meta, V = {}, {}, {}
    for lab, run, piece, fac in SOURCES:
        rec[lab], meta[lab] = load(run, piece, fac)
        V[lab] = val(rec[lab])

    lines = []
    out = lines.append
    out("TASK 1: F event-generator branch vs F devel branch; NLO0 - LO vs F (+A)")
    out("Integrated values = pymule Record.value (vegas integral) in microbarn.")
    out(f"F_devel multiplied by alpha^3 = {alpha**3:.10e} (devel: alpha=1 in McMule).")
    out("Errors of combinations: quadrature (independent seeds, no covariance).")
    out("")
    out("INPUT PIECES")
    for lab, run, piece, fac in SOURCES:
        v, e = V[lab]
        out(f"  {lab:10s} {run:28s} {piece:10s} {v:.12e} +- {e:.3e} ub  [{meta[lab]}]")
    out("")

    F_raw_ratio = V["F_devel"][0] / alpha**3 / V["F_evtgen"][0]
    out(f"Raw ratio F_devel(alpha=1)/F_evtgen = {F_raw_ratio:.8e};  1/alpha^3 = {alpha**-3:.8e}")
    out("")

    rows = []  # (name, value, error, pull, rel, comment)

    def check(name, terms, ref=None, comment=""):
        v, e = comb(terms)
        pull = v / e if e > 0 else float("nan")
        rel = v / ref if ref else float("nan")
        rows.append((name, v, e, pull, rel, comment))

    d = comb([(1, V["NLO0"]), (-1, V["LO"])])
    rows.append(("NLO0 - LO [complete run]", d[0], d[1], float("nan"), float("nan"), ""))
    for F in ("F_evtgen", "F_devel"):
        check(f"(a) (NLO0-LO) - {F}", [(1, V["NLO0"]), (-1, V["LO"]), (-1, V[F])],
              V[F][0], "relative to F")
        check(f"(b) (NLO0-LO) - ({F}+A)",
              [(1, V["NLO0"]), (-1, V["LO"]), (-1, V[F]), (-1, V["A"])],
              V[F][0] + V["A"][0], "relative to F+A; NLO0 of this run has NO A")
    check("F_devel - F_evtgen", [(1, V["F_devel"]), (-1, V["F_evtgen"])],
          V["F_evtgen"][0], "relative to F_evtgen")
    # replicas: which LO/NLO0 combination is used matters
    for nl in ("NLO0", "NLO0_rep"):
        for lo in ("LO", "LO_rep"):
            if (nl, lo) == ("NLO0", "LO"):
                continue
            check(f"(a') ({nl}-{lo}) - F_evtgen",
                  [(1, V[nl]), (-1, V[lo]), (-1, V["F_evtgen"])], V["F_evtgen"][0],
                  "replica combination")
    check("LO_rep - LO", [(1, V["LO_rep"]), (-1, V["LO"])], V["LO"][0],
          "two independent LO runs")
    check("NLO0_rep - NLO0", [(1, V["NLO0_rep"]), (-1, V["NLO0"])], V["NLO0"][0],
          "two independent NLO0 runs")
    for lo in ("LO", "LO_rep"):
        for F in ("F_evtgen", "F_devel"):
            check(f"(b') (NLO0_newA-{lo}) - ({F}+A)",
                  [(1, V["NLO0_newA"]), (-1, V[lo]), (-1, V[F]), (-1, V["A"])],
                  V[F][0] + V["A"][0], "new NLO0 incl. A (cda91c8), cross-run")

    out("INTEGRATED CHECKS (value, error in ub; pull = value/error; rel = value/reference)")
    out(f"  {'check':44s} {'value [ub]':>14s} {'error [ub]':>11s} {'value [mb]':>12s} {'pull':>8s} {'rel':>11s}")
    for name, v, e, pull, rel, c in rows:
        out(f"  {name:44s} {v:14.6e} {e:11.3e} {v/1000:12.3e} {pull:8.2f} {rel:11.3e}  {c}")
    out("")

    # ---------------- binning ----------------
    out("HISTOGRAM BINNING (identical bin centres as F_evtgen?)")
    ref = rec["F_evtgen"]
    for lab in rec:
        diff = [k for k in HISTS if not same_bins(hist(ref, k), hist(rec[lab], k))]
        out(f"  {lab:10s} " + ("all 9 identical" if not diff else "DIFFERENT: " + ", ".join(diff)))
    for lab in ("LO", "LO_rep"):
        h = hist(rec[lab], "th3")
        out(f"  {lab} th3: {len(h)-2} bins, first/last centre {h[1,0]:.6g}/{h[-2,0]:.6g} rad, "
            f"underflow {h[0,1]:.4g} ub/rad-bin")
    out("")

    # ---------------- per-bin F_devel vs F_evtgen ----------------
    out("PER-BIN F_devel vs F_evtgen (all bins incl. under/overflow with error>0)")
    for k in HISTS:
        ha, hb = hist(rec["F_devel"], k), hist(rec["F_evtgen"], k)
        if not same_bins(ha, hb):
            out(f"  {k}: SKIPPED, different binning")
            continue
        pulls, rel, m = perbin_stats(ha, hb)
        n = len(pulls)
        out(f"  {k:6s} nbins={n:4d} chi2/ndf={np.sum(pulls**2)/n:6.3f} "
            f"mean pull={pulls.mean():+.3f} rms pull={pulls.std():.3f} "
            f"max|pull|={np.abs(pulls).max():.2f} frac(|pull|>3)={np.mean(np.abs(pulls)>3):.4f} "
            f"median|rel diff|={np.nanmedian(np.abs(rel[m])):.2e}")
    out("")

    # ---------------- per-bin (NLO0-LO) vs F ----------------
    out("PER-BIN (NLO0 - LO) - F  (only where all three have identical binning)")
    for nl, lo in (("NLO0", "LO"), ("NLO0_rep", "LO_rep")):
        for F in ("F_evtgen", "F_devel"):
            for k in ("th3", "Qsq", "Emu"):
                hn, hl, hf = hist(rec[nl], k), hist(rec[lo], k), hist(rec[F], k)
                if not (same_bins(hn, hl) and same_bins(hn, hf)):
                    out(f"  {nl}-{lo} vs {F} {k}: SKIPPED (different bin grids, no rebinning)")
                    continue
                d = np.column_stack((hn[:, 0], hn[:, 1] - hl[:, 1], np.hypot(hn[:, 2], hl[:, 2])))
                pulls, rel, m = perbin_stats(d, hf)
                out(f"  {nl}-{lo} vs {F} {k:4s}: nbins={len(pulls)} chi2/ndf={np.sum(pulls**2)/len(pulls):.3f} "
                    f"mean pull={pulls.mean():+.3f}; median per-bin error of (NLO0-LO) / |F| = "
                    f"{np.nanmedian(d[m,2]/np.abs(hf[m,1])):.2e}")
    out("")

    txt = outval / f"{TAG}_F_branches_comparison.txt"
    txt.write_text("\n".join(lines) + "\n", encoding="utf-8")
    tsv = outval / f"{TAG}_F_branches_checks.tsv"
    with tsv.open("w", encoding="utf-8") as f:
        f.write("check\tvalue_ub\terror_ub\tpull\trel\tcomment\n")
        for lab in rec:
            f.write(f"{lab}\t{V[lab][0]:.16e}\t{V[lab][1]:.16e}\t\t\tinput piece\n")
        for name, v, e, pull, rel, c in rows:
            f.write(f"{name}\t{v:.16e}\t{e:.16e}\t{pull:.6g}\t{rel:.6g}\t{c}\n")
    print("\n".join(lines))
    print(f"Saved {txt}\nSaved {tsv}")

    for k in ("th3", "Qsq", "Emu"):
        ha, hb = hist(rec["F_devel"], k), hist(rec["F_evtgen"], k)
        if same_bins(ha, hb):
            path = outfig / f"{TAG}_F_devel_vs_F_evtgen_{k}.png"
            plot_compare(path, ha, hb, r"F devel ($\times\alpha^3$)", "F event-generator", k,
                         f"F: devel vs event-generator branch, {PRETTY[k]}",
                         r"F_\mathrm{devel}", r"F_\mathrm{evtgen}")
            print(f"Saved {path}")
        hn, hl, hf = (hist(rec["NLO0_rep"], k), hist(rec["LO_rep"], k),
                      hist(rec["F_evtgen"], k))
        if same_bins(hn, hl) and same_bins(hn, hf):
            d = np.column_stack((hn[:, 0], hn[:, 1] - hl[:, 1], np.hypot(hn[:, 2], hl[:, 2])))
            path = outfig / f"{TAG}_NLO0rep_minus_LOrep_vs_F_{k}.png"
            # NLO0 of this (old-build) run does NOT contain A -> NLO0 - LO must equal F.
            plot_compare(path, d, hf, "NLO0 - LO (18_11_NLO0R0, NLO0 without A)", "F event-generator", k,
                         f"NLO0 − LO vs F, {PRETTY[k]}",
                         r"(\mathrm{NLO0}-\mathrm{LO})", "F")
            print(f"Saved {path}")
        # Corrected NLO0 (24_11 run, contains A) -> NLO0 - LO must equal F + A.
        hn2, ha2 = hist(rec["NLO0_newA"], k), hist(rec["A"], k)
        if same_bins(hn2, hl) and same_bins(hn2, hf) and same_bins(hn2, ha2):
            d2 = np.column_stack((hn2[:, 0], hn2[:, 1] - hl[:, 1], np.hypot(hn2[:, 2], hl[:, 2])))
            fa = np.column_stack((hf[:, 0], hf[:, 1] + ha2[:, 1], np.hypot(hf[:, 2], ha2[:, 2])))
            path = outfig / f"{TAG}_NLO0withA_minus_LO_vs_FplusA_{k}.png"
            plot_compare(path, d2, fa, "NLO0 (24_11, with A) - LO (18_11_NLO0R0)", "F + A", k,
                         f"NLO0 − LO vs F + A, {PRETTY[k]}",
                         r"(\mathrm{NLO0}-\mathrm{LO})", r"(F+A)")
            print(f"Saved {path}")
        else:
            print(f"SKIP NLO0withA plot for {k}: binning differs")


if __name__ == "__main__":
    main()
