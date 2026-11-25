#!/usr/bin/env python3
"""Binning-safe McMule analysis for an interim meeting.

NLO0 contains Born + F, but excludes the separately calculated R and A.
Full NLO = NLO0 + R + A; Delta NLO = NLO0 + R + A - LO.
Independent checks: NLO0 == LO+F, and NLO0+R == LO+R+F.

Integrated cross sections are read from the per-piece Record.value, NEVER
from histogram sums. This is safe with incompatible histogram grids provided
that physics cuts and normalisations are otherwise identical.

Only compare/add histograms after requiring equal bin centres; never rebin.
"""

import argparse
import csv
import math
from array import array
from pathlib import Path

import numpy as np
from pymule import setup, sigma, mergefks, conv

PIECES = {"LO": "mp2mp0", "NLO0": "mp2mpNLO0", "R": "mp2mpR",
          "F": "mp2mpF", "A": "mp2mpA"}
# Units: McMule raw x / scale -> x units for display; density * scale
OBS = [
    ("th3", "Muon scattering angle", 1e-3, "Muon angle [mrad]", None),
    ("Emu", "Scattered muon energy", 1e3, "Muon energy [GeV]", None),
    ("costh3", "Cosine of muon angle", 1., "cos(theta_mu)", None),
    ("Qsq", "Momentum transfer", 1e6, "Q^2 [GeV^2]", None),
]
PHOTONS = [
    ("th5", "Photon polar angle", 1e-3, "Photon angle [mrad]", None),
    ("Eph", "Photon energy", 1e3, "Photon energy [GeV]", None),
    ("phi5", "Photon azimuth", 1., "Photon phi [rad]", None),
]
COLORS = {"LO": "#19866B", "NLO0": "#C879A3", "A": "#477E9E",
          "F": "#9772B5", "R": "#C18A42", "Delta_NLO0": "#477E9E",
          "NLO0_plus_R": "#C879A3", "NLO_full": "#C879A3",
          "Delta_NLO": "#C879A3", "Delta_NLO_no_A": "#477E9E"}


def hist_rows(record, key):
    """Return raw finite Nx3 histogram; report failure without masking it."""
    data = np.asarray(record[key], dtype=float)
    if data.ndim != 2 or data.shape[1] < 3:
        raise ValueError("expected Nx3 histogram")
    # Keep bin centres and physical weights; don't silently drop non-finite data
    finite = np.isfinite(data[:, 0])
    data = data[finite, :3]
    if len(data) < 2 or not np.all(np.isfinite(data)):
        raise ValueError(f"only {len(data)} finite-centred bins or non-finite weights")
    if not np.all(np.diff(data[:, 0]) > 0):
        raise ValueError("bin centres not strictly increasing")
    return data


def compatible(h1, h2):
    """Strict compatibility: shape and ALL bin centres must agree."""
    return h1.shape == h2.shape and np.allclose(
        h1[:, 0], h2[:, 0], rtol=1e-12, atol=1e-14)


def add_hists(a, b, sign=1):
    if not compatible(a, b):
        raise ValueError("different bin centres (no rebinning or interpolation)")
    return np.column_stack((a[:, 0], a[:, 1] + sign*b[:, 1],
                            np.hypot(a[:, 2], b[:, 2])))


def scalar(results, key):
    return results[key][:2]


def add_cs(results, keys, signs=None):
    signs = signs or [1] * len(keys)
    value = sum(results[k][0]*sgn for k, sgn in zip(keys, signs))
    error = math.sqrt(sum(results[k][1]**2 for k in keys))
    return value, error


def fmt(value, error, unit="mb"):
    return f"{value:.12g} +/- {error:.3g} {unit}"


def draw_single(path, h, info, label):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["text.usetex"] = False  # support Greek symbols without external LaTeX
    key, title, scale, xlabel, _ = info
    x, y, err = h[:, 0]/scale, h[:, 1]*scale, np.abs(h[:, 2]*scale)
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    c = COLORS.get(label, "#477E9E")
    ax.step(x, y, color=c, where="mid", lw=1.6, label=label)
    ax.fill_between(x, y-err, y+err, color=c, alpha=.17, step="mid")
    ax.axhline(0, color="gray", lw=.7)
    xunits = xlabel.split("[")[-1].rstrip("]") if "[" in xlabel else ""
    ylab = r"$d\sigma/dx$ [$\mu$b" + (f" / {xunits}]" if xunits else "]")
    ax.set(title=f"{title} — {label}", xlabel=xlabel, ylabel=ylab)
    ax.grid(alpha=.20)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def draw_pair(path, h_lo, h_nlo0, info):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    plt.rcParams["text.usetex"] = False  # support Greek symbols without external LaTeX
    key, title, scale, xlabel, _ = info
    x = h_lo[:, 0]/scale
    fig, ax = plt.subplots(figsize=(8.6, 5.3))
    for label, h in (("LO", h_lo), ("NLO0 (without A)", h_nlo0)):
        c = COLORS["LO"] if label == "LO" else COLORS["NLO0"]
        ax.step(x, h[:, 1]*scale, where="mid", color=c, lw=1.65, label=label)
    ax.set(title=f"{title} — LO versus NLO0 (no A)", xlabel=xlabel,
           ylabel="Differential cross section (scaled to x units)")
    ax.grid(alpha=.2)
    ax.legend()
    fig.tight_layout()
    fig.savefig(path, dpi=170)
    plt.close(fig)


def root_write(path, record_rows, log):
    """Optional export of independently valid histograms without rebinning."""
    try:
        import ROOT
    except ImportError:
        log.append("ROOT unavailable: skipped ROOT export (PNG and CS unaffected)")
        return
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(str(path), "RECREATE")
    if not f or f.IsZombie():
        log.append(f"Could not create ROOT file: {path}")
        return
    try:
        for name, (h, info) in record_rows.items():
            key, _, scale, xlabel, _ = info
            x = h[:, 0]/scale
            mids = (x[1:] + x[:-1])/2
            edges = np.r_[2*x[0]-mids[0], mids, 2*x[-1]-mids[-1]]
            if not np.all(np.diff(edges) > 0):
                log.append(f"ROOT skip {name}: invalid edges")
                continue
            th = ROOT.TH1D(name, name, len(x), array("d", edges.tolist()))
            th.Sumw2()
            th.GetXaxis().SetTitle(xlabel)
            for i, row in enumerate(h, 1):
                th.SetBinContent(i, float(row[1]*scale))
                th.SetBinError(i, float(abs(row[2]*scale)))
            f.cd()
            th.Write()
            log.append("ROOT saved " + name)
    finally:
        f.Close()


def main():
    parser = argparse.ArgumentParser(description="Meeting-safe McMule CS + valid standalone distributions")
    parser.add_argument("--rundir", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--outdir-vals", required=True)
    parser.add_argument("--savename", required=True)
    parser.add_argument("--pieces", nargs="+", default=list(PIECES))
    parser.add_argument("--cs-only", action="store_true", help="Compute integrated CS without any plotting")
    args = parser.parse_args()
    directory = Path(args.rundir) / args.run / "out"
    outfig = Path(args.outdir)
    outval = Path(args.outdir_vals)
    outfig.mkdir(parents=True, exist_ok=True)
    outval.mkdir(parents=True, exist_ok=True)
    basename = f"{args.savename}_{args.run}"
    print(f"Reading McMule results from: {directory}", flush=True)
    if not directory.is_dir():
        parser.error(f"input folder does not exist: {directory}")
    setup(folder=str(directory))
    log = [f"Input: {directory}", "No rebinning or interpolation was performed."]
    results = {}  # in mb
    records = {}
    selected = list(dict.fromkeys("LO" if p == "0" else p for p in args.pieces))
    for key in selected:
        if key not in PIECES:
            log.append(f"UNKNOWN piece {key}; skipped")
            continue
        try:
            # Merge ONE physical piece at a time: never align A/F histograms
            # with NLO0/LO to obtain the integrated CS.
            rec = mergefks(sigma(PIECES[key])) * conv
            cs_microbarn, err_microbarn = map(float, rec.value[:2])
            if not (np.isfinite(cs_microbarn) and np.isfinite(err_microbarn)):
                raise ValueError("non-finite integrated cross section")
            results[key] = (cs_microbarn/1000, abs(err_microbarn)/1000)
            records[key] = rec
            log.append(f"LOADED {key} ({PIECES[key]}): {fmt(*results[key])}")
        except Exception as exc:
            log.append(f"FAILED {key}: {type(exc).__name__}: {exc}")

    derived = {}
    # Integrated values are independent of histogram binning.  F is already
    # inside NLO0, so adding F to NLO0 would double count it.
    if all(k in results for k in ("NLO0", "R")):
        derived["NLO_no_A = NLO0 + R"] = add_cs(results, ["NLO0", "R"])
    if all(k in results for k in ("NLO0", "R", "LO")):
        derived["Delta_NLO_no_A = NLO0+R-LO"] = add_cs(
            results, ["NLO0", "R", "LO"], [1, 1, -1])
    if all(k in results for k in ("NLO0", "R", "A")):
        derived["NLO_total = NLO0 + R + A"] = add_cs(results, ["NLO0", "R", "A"])
    if all(k in results for k in ("NLO0", "R", "A", "LO")):
        derived["Delta_NLO = NLO0+R+A-LO"] = add_cs(
            results, ["NLO0", "R", "A", "LO"], [1, 1, 1, -1])
        lo, elo = results["LO"]
        nlo, enlo = derived["NLO_total = NLO0 + R + A"]
        if lo != 0:
            frac = 100*(nlo/lo - 1)
            # Approximation: the integrated piece estimates are independent.
            # Covariance, if any, must be evaluated separately.
            frac_err = 100*math.sqrt((enlo/lo)**2 + (nlo*elo/lo**2)**2)
            derived["Relative Delta_NLO / LO [%]"] = (frac, frac_err)
            derived["NLO / LO"] = (nlo/lo, math.sqrt((enlo/lo)**2 + (nlo*elo/lo**2)**2))
    # NLO0 has Born + F but no R; these are independent cross-checks.
    if all(k in results for k in ("NLO0", "LO", "F")):
        derived["CHECK F-(NLO0-LO)"] = add_cs(
            results, ["F", "NLO0", "LO"], [1, -1, 1])
    if all(k in results for k in ("NLO0", "LO")):
        derived["Delta_NLO0 = NLO0 - LO"] = add_cs(results, ["NLO0", "LO"], [1, -1])
    if all(k in results for k in ("R", "F")):
        derived["R + F (integrated)"] = add_cs(results, ["R", "F"])
    if all(k in results for k in ("LO", "R", "F")):
        derived["0 + R + F (integrated)"] = add_cs(results, ["LO", "R", "F"])
    if all(k in results for k in ("LO", "R", "F", "NLO0")):
        # Both R terms cancel algebraically, including their statistical
        # fluctuations. Don't double-count R's uncertainty here.
        check_value, check_error = add_cs(
            results, ["F", "NLO0", "LO"], [1, -1, 1])
        derived["CHECK (R+F)-(NLO0+R-LO)"] = (check_value, check_error)
        derived["CHECK (0+R+F)-(NLO0+R)"] = (check_value, check_error)

    rows = [(k, *v) for k, v in {**results, **derived}.items()]
    dest_tsv = outval / f"{basename}_MEETING_cross_sections.tsv"
    dest_txt = outval / f"{basename}_MEETING_cross_sections.txt"
    with dest_tsv.open("w", newline="", encoding="utf-8") as f:
        writer = csv.writer(f, delimiter="\t")
        writer.writerow(["contribution", "value", "uncertainty", "units"])
        for key, value, error in rows:
            writer.writerow([key, f"{value:.16e}", f"{error:.16e}",
                             "%" if "[%]" in key else ("dimensionless" if key=="NLO / LO" else "mb")])
    with dest_txt.open("w", encoding="utf-8") as f:
        f.write(f"McMule integrated cross sections — {args.run}\n")
        f.write("Full NLO = NLO0 + R + A; Delta_NLO = NLO0 + R + A - LO.\n")
        f.write("NLO0 already contains Born + F, but excludes R and A.\n")
        f.write("CHECK: NLO0 == LO + F, and NLO0 + R == LO + R + F.\n")
        f.write("Unit is mb unless indicated otherwise. Coupling already included in integration.\n")
        f.write("Error estimates for sums and ratios assume independent samples; covariance not included.\n")
        f.write("Histogram bin changes do not alter these CS if physical cuts/normalisations match.\n\n")
        for key, value, error in rows:
            unit = "%" if "[%]" in key else ("" if key=="NLO / LO" else "mb")
            f.write(f"{key:37} {fmt(value, error, unit)}\n")
    print("\nINTEGRATED CROSS SECTIONS / CHECKS")
    print("Definition: NLO_total = NLO0 + R + A (NLO0 already includes LO + F).")
    for key, value, error in rows:
        unit = "%" if "[%]" in key else ("" if key=="NLO / LO" else "mb")
        print(f"{key:37} {fmt(value, error, unit)}")
    print(f"\nSaved: {dest_txt}\nSaved: {dest_tsv}", flush=True)

    if not args.cs_only:
        root_rows = {}
        hist_dict = {}
        for info in OBS + PHOTONS:
            key = info[0]
            for piece, rec in records.items():
                if key in ("th5", "Eph", "phi5") and piece != "R":
                    continue  # photon plots: FKS real piece only
                try:
                    h = hist_rows(rec, key)
                    hist_dict[piece, key] = h
                    out = outfig / f"{basename}_{piece}_{key}_MEETING.png"
                    draw_single(out, h, info, piece)
                    root_rows[f"{piece}_{key}"] = (h, info)
                    log.append(f"PLOTTED standalone {piece}:{key} ({len(h)} bins)")
                except Exception as exc:
                    log.append(f"SKIPPED standalone {piece}:{key}: {exc}")
            if key in ("th5", "Eph", "phi5"):
                continue
            lo, nlo = hist_dict.get(("LO", key)), hist_dict.get(("NLO0", key))
            if lo is not None and nlo is not None and compatible(lo, nlo):
                try:
                    draw_pair(outfig / f"{basename}_LO_vs_NLO0_{key}_MEETING.png", lo, nlo, info)
                    d0 = add_hists(nlo, lo, sign=-1)
                    draw_single(outfig / f"{basename}_Delta_NLO0_{key}_MEETING.png",
                                d0, info, "Delta_NLO0")
                    root_rows[f"Delta_NLO0_{key}"] = (d0, info)
                    log.append(f"PLOTTED LO vs NLO0 + Delta_NLO0 for {key}")
                except Exception as exc:
                    log.append(f"SKIPPED LO/NLO0 for {key}: {exc}")
            else:
                log.append(f"SKIPPED LO vs NLO0 and Delta_NLO0 for {key}: missing or incompatible grids")
            # Compose only if all original bin centres agree.  Unlike the
            # integrated cross sections, these curves require matching grids.
            combos = [
                ("R_plus_F", "R", "F"),
                ("NLO0_plus_R", "NLO0", "R"),
            ]
            for name, p1, p2 in combos:
                h1, h2 = hist_dict.get((p1, key)), hist_dict.get((p2, key))
                if h1 is None or h2 is None or not compatible(h1, h2):
                    log.append(f"SKIPPED {name}:{key}: missing/different bin grids")
                    continue
                try:
                    res = add_hists(h1, h2)
                    draw_single(outfig / f"{basename}_{name}_{key}_MEETING.png", res, info, name)
                    root_rows[f"{name}_{key}"] = (res, info)
                    log.append(f"PLOTTED {name}:{key} (original bins)")
                except Exception as exc:
                    log.append(f"SKIPPED {name}:{key}: {type(exc).__name__}: {exc}")

            # Differential delta without A is still useful for the meeting.
            nr = root_rows.get(f"NLO0_plus_R_{key}")
            if nr is not None and lo is not None and compatible(nr[0], lo):
                try:
                    d_no_a = add_hists(nr[0], lo, -1)
                    draw_single(outfig / f"{basename}_Delta_NLO_no_A_{key}_MEETING.png",
                                d_no_a, info, "Delta_NLO_no_A")
                    root_rows[f"Delta_NLO_no_A_{key}"] = (d_no_a, info)
                    log.append(f"PLOTTED Delta_NLO_no_A:{key}")
                except Exception as exc:
                    log.append(f"SKIPPED Delta_NLO_no_A:{key}: {type(exc).__name__}: {exc}")
            else:
                log.append(f"SKIPPED Delta_NLO_no_A:{key}: missing/incompatible bins")

            # Full NLO = NLO0 + R + A; DO NOT add the separate F again.
            a = hist_dict.get(("A", key))
            if nr is not None and a is not None and compatible(nr[0], a):
                try:
                    nf = add_hists(nr[0], a)
                    draw_single(outfig / f"{basename}_NLO_full_{key}_MEETING.png",
                                nf, info, "NLO_full")
                    root_rows[f"NLO_full_{key}"] = (nf, info)
                    log.append(f"PLOTTED NLO_full:{key}")
                    if lo is not None and compatible(nf, lo):
                        d = add_hists(nf, lo, -1)
                        draw_single(outfig / f"{basename}_Delta_NLO_{key}_MEETING.png",
                                    d, info, "Delta_NLO")
                        root_rows[f"Delta_NLO_{key}"] = (d, info)
                        log.append(f"PLOTTED Delta_NLO:{key}")
                    else:
                        log.append(f"SKIPPED Delta_NLO:{key}: LO grid incompatible")
                except Exception as exc:
                    log.append(f"SKIPPED NLO_full/Delta_NLO:{key}: {type(exc).__name__}: {exc}")
            else:
                log.append(f"SKIPPED NLO_full and Delta_NLO:{key}: A/NLO0/R grids incompatible or missing")
        root_write(outfig / f"{basename}_MEETING_histograms.root", root_rows, log)
    else:
        log.append("--cs-only: no histograms attempted")
    logfile = outval / f"{basename}_MEETING_report.log"
    logfile.write_text("\n".join(log)+"\n", encoding="utf-8")
    print(f"Saved report: {logfile}")

if __name__ == "__main__":
    main()