#!/usr/bin/env python3
"""Relative MC error per theta_mu bin for LO / NLO0 / R in several McMule runs.

Shows where VEGAS stops sampling densely: the per-bin relative error of the
2->2 pieces (LO, NLO0) jumps by a factor 20-70 somewhere between 1.3 and
1.6 mrad, at a different place in every run. That region is where the
"fluctuations around 1.4 mrad" in NLO0 - LO and F + R come from.

Output: McMule_Results/Figures/26_11_meeting/26_11_meeting_vegas_relerr_th3.png
        McMule_Results/Vals/26_11_meeting/26_11_meeting_vegas_relerr_jumps.tsv
"""
import sys
from pathlib import Path

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, str(Path(__file__).resolve().parent))
from plot_style import apply_style, COLORS  # noqa: E402
import compare_F_branches as C  # noqa: E402

apply_style()

RUNS = [
    ("evtgen_18_11_complete", "old run (18_11_complete)", "-"),
    ("evtgen_18_11_NLO0R0", "replica (18_11_NLO0R0)", "--"),
    ("evtgen_24_11_NLO0corrected", "new NLO0 with A (24_11)", ":"),
]
PIECES = [("mp2mp0", "LO", COLORS["LO"]), ("mp2mpNLO0", "NLO0", COLORS["NLO"]),
          ("mp2mpR", "R", COLORS["R"])]
OUT_FIG = C.MOMOS / "Figures" / C.TAG / f"{C.TAG}_vegas_relerr_th3.png"
OUT_TSV = C.MOMOS / "Vals" / C.TAG / f"{C.TAG}_vegas_relerr_jumps.tsv"


def first_jump(x, r, factor=5.0, window=20, xmin=0.2):
    for i in range(window, len(r)):
        if x[i] > xmin and r[i] > factor * np.median(r[i - window:i]):
            return x[i]
    return np.nan


def main():
    fig, ax = plt.subplots(figsize=(11, 6.5))
    rows = []
    for run, run_label, ls in RUNS:
        for piece, plabel, color in PIECES:
            try:
                rec = C.load(run, piece, 1.0)[0]
            except Exception:
                continue
            h = C.hist(rec, "th3")[1:-1]
            x, y, e = h[:, 0] * 1e3, h[:, 1], h[:, 2]
            with np.errstate(all="ignore"):
                r = np.abs(e / y)
            ok = np.isfinite(r) & (r > 0)
            x, r = x[ok], r[ok]
            jump = first_jump(x, r)
            lo = np.median(r[(x > 0.8) & (x < 1.2)])
            hi = np.median(r[(x > 1.6) & (x < 2.0)])
            rows.append((run, plabel, jump, lo, hi, hi / lo))
            ax.step(x, r, where="mid", color=color, ls=ls, lw=1.6,
                    label=f"{plabel}, {run_label}")
    ax.axvspan(1.3, 1.6, color=COLORS["HIGHLIGHT"], alpha=0.35, lw=0,
               label="jump region 1.3-1.6 mrad")
    ax.set_yscale("log")
    ax.set_xlim(0, 4.05)
    ax.set_xlabel(r"Muon angle $\theta_\mu$ [mrad]")
    ax.set_ylabel("relative MC error per bin")
    ax.set_title("Where VEGAS stops sampling densely")
    ax.legend(ncol=2, fontsize=11, loc="lower right")
    OUT_FIG.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(OUT_FIG)
    plt.close(fig)
    OUT_TSV.parent.mkdir(parents=True, exist_ok=True)
    with OUT_TSV.open("w") as f:
        f.write("run\tpiece\tfirst_jump_mrad\trelerr_0.8-1.2\trelerr_1.6-2.0\tratio\n")
        for row in rows:
            f.write("\t".join(str(v) if isinstance(v, str) else f"{v:.4g}" for v in row) + "\n")
    print("Saved", OUT_FIG)
    print("Saved", OUT_TSV)
    for row in rows:
        print("\t".join(str(v) if isinstance(v, str) else f"{v:.3g}" for v in row))


if __name__ == "__main__":
    main()
