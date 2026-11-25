#!/usr/bin/env python3
"""Analyze weighted elastic mu-p McMule distributions, with optional pieces.

Run from command_analyze_mcmule4vecs.sh with PIECES set in the sourced .conf.
No adaptive rebinning is performed; all McMule histogram bins are retained.

NLO0 contains Born + F + A in the new event-generator branch (cda91c88).
Full NLO = NLO0 + R; Delta NLO = NLO0 + R - LO.
R is FKS-subtracted, not a standalone physical photon event rate.
Use only with the new NLO0 definition.
"""

import os
os.environ.setdefault("MPLBACKEND", "Agg")
import matplotlib
matplotlib.use("Agg")

import argparse
from array import array
import math
from pathlib import Path
import sys
import traceback

import numpy as np
from pymule import setup, sigma, mergefks, conv
from plotting import (
    write_file_with_values, save_real_photon,
    save_delta_nlo, save_analytic_lo_check,
)

# Avoid external LaTeX dependencies (e.g. Unicode sigma errors on cluster).
matplotlib.rcParams["text.usetex"] = False

PIECE_NAMES = {
    "LO": "mp2mp0",
    "NLO0": "mp2mpNLO0",
    "R": "mp2mpR",
}


# Comparison palette: green LO / muted pink NLO.
COLOR_LO = "#19866B"
COLOR_NLO = "#C879A3"
COLOR_R = "#9772B5"
COLOR_DIFF = "#477E9E"


class Tee:
    """Keep the console output, while writing it to a file as well."""

    def __init__(self, path):
        self.terminal = sys.stdout
        self.file = Path(path).open("w", encoding="utf-8")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


class Report:
    """Track declared, successfully loaded, missing and skipped outputs."""

    def __init__(self, selected):
        self.selected = set(selected)
        self.loaded = []
        self.failures = {}
        self.created = []
        self.skipped = []
        self.messages = []
        self.raw_records = {}
        self.combo_sources = {}

    def ok(self, output):
        self.created.append(str(output))

    def skip(self, output, reason):
        self.skipped.append((str(output), str(reason)))
        print(f"SKIPPED {output}: {reason}")

    def note(self, message):
        self.messages.append(str(message))
        print(message)

    def write(self, path, *, error=None):
        path = Path(path)
        with path.open("w", encoding="utf-8") as f:
            f.write("McMule piece availability and plot status\n")
            f.write("=" * 65 + "\n")
            for piece, process in PIECE_NAMES.items():
                if piece in self.loaded:
                    status = "AVAILABLE (loaded)"
                elif piece not in self.selected:
                    status = "NOT SELECTED in config"
                else:
                    status = "REQUESTED BUT UNAVAILABLE: " + self.failures.get(piece, "unknown")
                f.write(f"{piece:5s} {process:13s}  {status}\n")
            f.write("\nImportant definitions:\n")
            f.write("LO = mp2mp0; NLO0 = LO+F+A; R = mp2mpR.\n")
            f.write("NLO_full = NLO0 + R; Delta_NLO = NLO0 + R - LO.\n")
            f.write("R is FKS-subtracted; no event-level data are read.\n")
            f.write("Histogram uncertainty bands for differences ignore covariance.\n")
            f.write("No adaptive rebinning; original McMule bins are used.\n")
            for msg in self.messages:
                f.write("NOTE: " + msg + "\n")
            f.write(f"\nCreated plots/histograms: {len(self.created)}\n")
            for item in self.created:
                f.write("  + " + item + "\n")
            f.write(f"\nSkipped plots/histograms: {len(self.skipped)}\n")
            for name, reason in self.skipped:
                f.write(f"  - {name}: {reason}\n")
            if error is not None:
                f.write("\nANALYSIS ERROR:\n" + str(error) + "\n")
        print(f"Piece/plot report: {path}")


def merged(*records):
    """Combine FKS records and convert to microbarn (physical alpha already included)."""
    return mergefks(*records) * conv


def build_contributions(args, report):
    """Only load selected pieces. A missing piece does not stop all plots."""
    raw = {}
    for key, process in PIECE_NAMES.items():
        if key not in report.selected:
            continue
        try:
            raw[key] = sigma(process)
            report.loaded.append(key)
            report.note(f"Loaded {key} ({process})")
        except Exception as exc:
            report.failures[key] = f"{type(exc).__name__}: {exc}"
            report.note(f"WARNING: {key} ({process}) could not be loaded: {exc}")

    # Keep exact source pieces so that composite histograms can be checked
    # against their native bin centres (mergefks may silently drop bins).
    report.raw_records = raw
    built = {}

    def make(name, required, fn):
        missing = [p for p in required if p not in raw]
        if missing:
            report.note(f"Unavailable {name}: missing {', '.join(missing)}")
            return
        try:
            built[name] = fn()
            report.combo_sources[name] = tuple(required)
        except Exception as exc:
            report.note(f"WARNING: could not build {name}: {type(exc).__name__}: {exc}")

    make("LO", ("LO",), lambda: merged(raw["LO"]))
    make("NLO0", ("NLO0",), lambda: merged(raw["NLO0"]))
    make("R", ("R",), lambda: merged(raw["R"]))
    make("NLO_full", ("NLO0", "R"),
         lambda: merged(raw["NLO0"], raw["R"]))

    def difference(name, positive, negative):
        if positive not in built or negative not in built:
            report.note(f"Unavailable {name}: needs {positive}, {negative}")
            return
        try:
            built[name] = built[positive] - built[negative]
            report.combo_sources[name] = tuple(dict.fromkeys(
                report.combo_sources.get(positive, ()) +
                report.combo_sources.get(negative, ())))
        except Exception as exc:
            report.note(f"WARNING: could not calculate {name}: {exc}")

    difference("Delta_NLO0", "NLO0", "LO")  # F+A
    difference("Delta_NLO", "NLO_full", "LO")  # F+A+R
    return built


def bin_rows(hist, cfg):
    """Convert McMule [raw center, density, error] to displayed units."""
    arr = np.asarray(hist, dtype=float)
    if arr.ndim != 2 or arr.shape[1] < 3:
        raise ValueError("Expected Nx3 [bin center, density, statistical error]")
    rows = arr[np.isfinite(arr[:, 0])]
    if len(rows) < 2:
        raise ValueError("Fewer than two finite histogram bins (possibly incompatible McMule binning)")
    if not np.all(np.isfinite(rows[:, :3])):
        raise ValueError("Non-finite histogram weight or error")
    if not np.all(np.diff(rows[:, 0]) > 0):
        raise ValueError("Non-increasing histogram bin centers")
    scale = float(cfg["scale"])
    return rows[:, 0] / scale, rows[:, 1] * scale, np.abs(rows[:, 2] * scale)


def get_hist(record, cfg, report, file_label, record_name):
    if record is None:
        return None
    try:
        hist = record[cfg["key"]]
        x, _, _ = bin_rows(hist, cfg)
        # Never trust a merged McMule histogram unless all underlying pieces
        # have *identical original bin centres*. Otherwise mergefks can give
        # a truncated histogram even when >=2 common centres happen to remain.
        sources = report.combo_sources.get(record_name, ())
        for piece in sources:
            original = report.raw_records[piece][cfg["key"]]
            xp, _, _ = bin_rows(original, cfg)
            if len(xp) != len(x) or not np.allclose(xp, x, rtol=1e-12, atol=1e-13):
                raise ValueError(
                    f"incompatible original bin centres: {piece} vs "
                    f"{record_name} (source {len(xp)} bins; combined {len(x)} bins)")
        return hist
    except Exception as exc:
        report.skip(file_label, f"histogram {cfg['key']} invalid or unavailable: {exc}")
        return None


def save_single(path, hist, cfg, label, color, *, title=None):
    import matplotlib.pyplot as plt
    x, y, e = bin_rows(hist, cfg)
    fig, ax = plt.subplots(figsize=(9, 5.5))
    ax.step(x, y, where="mid", color=color, linewidth=1.8, label=label)
    ax.fill_between(x, y-e, y+e, step="mid", color=color, alpha=0.18)
    ax.axhline(0, color="#777777", lw=0.7, alpha=0.6)
    ax.set(xlabel=cfg["xlabel"], ylabel=cfg["ylabel"],
           title=title if title else cfg["title"])
    if cfg.get("xlim"):
        ax.set_xlim(*cfg["xlim"])
    ax.legend()
    ax.grid(alpha=0.22)
    fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)


def aligned(h1, h2, cfg):
    x1, y1, e1 = bin_rows(h1, cfg)
    x2, y2, e2 = bin_rows(h2, cfg)
    if len(x1) != len(x2) or not np.allclose(x1, x2, rtol=1e-11, atol=1e-13):
        raise ValueError(f"{cfg['key']}: incompatible histogram bin centers")
    return x1, y1, e1, y2, e2


def pair_and_k(path, lo_hist, nlo_hist, cfg, label_nlo):
    """LO green, NLO muted pink; K = (NLO - LO)/NLO (user convention)."""
    import matplotlib.pyplot as plt
    x, yl, el, yn, en = aligned(lo_hist, nlo_hist, cfg)
    # Report ratio data in displayed x-units, even if normal histogram uses raw x.
    k = np.full_like(yn, np.nan)
    ek = np.full_like(yn, np.nan)
    good = np.isfinite(yn) & np.isfinite(yl) & (yn != 0)
    with np.errstate(divide="ignore", invalid="ignore"):
        k[good] = (yn[good] - yl[good]) / yn[good]
        # LO and the full NLO are treated as statistically independent.
        ek[good] = np.sqrt((el[good]/yn[good])**2 +
                           (yl[good]*en[good]/yn[good]**2)**2)
    kr = np.column_stack((x, k, ek))

    show_k = cfg.get("show_k", True)
    if show_k:
        fig, (ax, axk) = plt.subplots(2, 1, figsize=(9, 7.1), sharex=False,
                                    gridspec_kw={"height_ratios": [2.5, 1]},
                                    constrained_layout=True)
    else:
        fig, ax = plt.subplots(figsize=(9, 5.5))
        axk = None
    ax.step(x, yl, where="mid", lw=1.7, color=COLOR_LO, label="LO")
    ax.step(x, yn, where="mid", lw=1.7, color=COLOR_NLO, label=label_nlo)
    ax.set_ylabel(cfg["ylabel"])
    ax.set_title(cfg["title"])
    ax.legend()
    ax.grid(alpha=0.22)
    if cfg.get("xlim"):
        ax.set_xlim(*cfg["xlim"])
    if axk is not None:
        axk.step(x, np.ma.masked_invalid(k), where="mid", color=COLOR_DIFF,
                 linewidth=1.4)
        axk.fill_between(x, np.ma.masked_invalid(k-ek),
                         np.ma.masked_invalid(k+ek), step="mid",
                         color=COLOR_DIFF, alpha=0.16)
        axk.axhline(0, color="#777777", linewidth=0.8)
        axk.set_ylabel(r"$(\mathrm{NLO}-\mathrm{LO})/\mathrm{NLO}$")
        axk.set_xlabel(cfg["xlabel"])
        axk.grid(alpha=0.22)
        left, right = ax.get_xlim()
        if cfg.get("k_xmax") is not None:
            right = min(right, cfg["k_xmax"])
        axk.set_xlim(left, right)
    else:
        ax.set_xlabel(cfg["xlabel"])
        fig.tight_layout()
    fig.savefig(path, dpi=180)
    plt.close(fig)
    return kr if show_k else None


def root_axes(cfg, *, ratio=False):
    titles = {
        "th3": ("#theta_{#mu} (mrad)", "d#sigma/d#theta_{#mu} (#mub/mrad)"),
        "Emu": ("E_{#mu} (GeV)", "d#sigma/dE_{#mu} (#mub/GeV)"),
        "costh3": ("cos#theta_{#mu}", "d#sigma/dcos#theta_{#mu} (#mub)"),
        "Qsq": ("Q^{2} (GeV^{2})", "d#sigma/dQ^{2} (#mub/GeV^{2})"),
        "th5": ("#theta_{#gamma} (mrad)", "d#sigma_{R}/d#theta_{#gamma} (#mub/mrad)"),
        "Eph": ("E_{#gamma} (GeV)", "d#sigma_{R}/dE_{#gamma} (#mub/GeV)"),
        "phi5": ("#phi_{#gamma} (rad)", "d#sigma_{R}/d#phi_{#gamma} (#mub/rad)"),
    }
    x, y = titles[cfg["key"]]
    return x, "(NLO-LO)/NLO" if ratio else y


def queue_root(entries, name, hist, cfg, *, ratio=False):
    entries.append(dict(name=name, hist=np.asarray(hist), cfg=cfg, ratio=ratio))


def save_root(path, histograms, report):
    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError("PyROOT unavailable; source your ROOT environment") from exc
    ROOT.gROOT.SetBatch(True)
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    f = ROOT.TFile.Open(str(path), "RECREATE")
    if not f or f.IsZombie():
        raise OSError(f"Could not create ROOT file {path}")
    count = 0
    try:
        for item in histograms:
            name, cfg, ratio = item["name"], item["cfg"], item["ratio"]
            hdata = item["hist"]
            if hdata.ndim != 2 or hdata.shape[1] < 3:
                report.skip(f"ROOT {name}", "invalid Nx3 histogram")
                continue
            raw = hdata[np.isfinite(hdata[:, 0])]
            if len(raw) < 2:
                report.skip(f"ROOT {name}", "fewer than 2 finite bins")
                continue
            scale = 1.0 if ratio else float(cfg["scale"])
            x = raw[:, 0] / scale
            mid = 0.5 * (x[1:] + x[:-1])
            edges = np.r_[2*x[0]-mid[0], mid, 2*x[-1]-mid[-1]]
            if not np.all(np.isfinite(edges)) or not np.all(np.diff(edges) > 0):
                report.skip(f"ROOT {name}", "invalid bin edges")
                continue
            x_title, y_title = root_axes(cfg, ratio=ratio)
            h = ROOT.TH1D(name, f"{cfg['title']} - {name}", len(edges)-1,
                          array("d", edges.tolist()))
            h.Sumw2()
            h.GetXaxis().SetTitle(x_title)
            h.GetYaxis().SetTitle(y_title)
            for i, (_, value, err, *_) in enumerate(raw, start=1):
                if np.isfinite(value) and np.isfinite(err):
                    h.SetBinContent(i, float(value) * scale)
                    h.SetBinError(i, abs(float(err)) * scale)
            if not ratio:
                underflow = hdata[hdata[:, 0] == -np.inf]
                overflow = hdata[hdata[:, 0] == np.inf]
                for rows, idx in ((underflow, 0), (overflow, len(edges))):
                    if len(rows) and np.all(np.isfinite(rows[0, 1:3])):
                        h.SetBinContent(idx, float(rows[0, 1]) * scale)
                        h.SetBinError(idx, abs(float(rows[0, 2])) * scale)
            f.cd()
            h.Write()
            report.ok(f"ROOT {name}")
            count += 1
    finally:
        f.Close()
    print(f"Saved {count} TH1D histograms to {path}")


def save_cross_sections(path_base, built, report):
    """Calculate all integrated CS directly from piece.value, not histograms.

    Histogram grids can be different (e.g. A or F with modified x limits).
    The integrated scalar CS still permit independent combinations, provided
    that the cuts and normalisation are identical for all pieces.
    """
    scalar = {}
    for name in ("LO", "NLO0", "R"):
        if name not in built:
            continue
        try:
            val, err = map(float, built[name].value[:2])
            if not (np.isfinite(val) and np.isfinite(err)):
                raise ValueError("non-finite integrated CS")
            scalar[name] = (val / 1000.0, abs(err) / 1000.0)
        except Exception as exc:
            report.note(f"WARNING: cross section missing for {name}: {exc}")

    def combination(label, signed_pieces):
        if any(name not in scalar for name, _ in signed_pieces):
            return
        val = sum(sign * scalar[name][0] for name, sign in signed_pieces)
        err = math.sqrt(sum(scalar[name][1] ** 2 for name, _ in signed_pieces))
        scalar[label] = (val, err)

    combination("NLO_full", (("NLO0", 1), ("R", 1)))
    combination("Delta_NLO0", (("NLO0", 1), ("LO", -1)))
    combination("Delta_NLO", (("NLO0", 1), ("R", 1), ("LO", -1)))

    if "NLO_full" in scalar and "LO" in scalar and scalar["LO"][0] != 0:
        nlo, enlo = scalar["NLO_full"]
        lo, elo = scalar["LO"]
        frac = 100 * (nlo / lo - 1)
        efrac = 100 * math.hypot(enlo / lo, nlo * elo / lo ** 2)
        scalar["Relative_Delta_NLO_over_LO_percent"] = (frac, efrac)
        scalar["NLO_over_LO"] = (nlo / lo, efrac / 100)

    rows = list(scalar.items())
    for name, (val, err) in rows:
        units = "%" if name.endswith("_percent") else ("" if name == "NLO_over_LO" else "mb")
        print(f"{name:39s} {val:.12g} +/- {err:.3g} {units}")
    for suffix in ("txt", "tsv"):
        path = path_base.with_suffix("." + suffix)
        with path.open("w", encoding="utf-8") as f:
            if suffix == "tsv":
                f.write("contribution\tvalue\tstat_error\tunit\n")
                for name, (val, err) in rows:
                    units = "%" if name.endswith("_percent") else ("dimensionless" if name == "NLO_over_LO" else "mb")
                    f.write(f"{name}\t{val:.16e}\t{err:.16e}\t{units}\n")
            else:
                f.write("McMule integrated cross sections; physical alpha included.\n")
                f.write("NLO0 = LO+F+A; FULL NLO = NLO0+R.\n")
                f.write("Delta_NLO = NLO0+R-LO. Neither F nor A is added separately.\n")
                f.write("CS calculated directly from per-piece integrated values, independent of histograms.\n")
                f.write("All error combinations assume independent statistical errors; covariance ignored.\n")
                f.write("Integrals require the same physical cuts and normalisation across pieces.\n\n")
                for name, (val, err) in rows:
                    units = "%" if name.endswith("_percent") else ("" if name == "NLO_over_LO" else "mb")
                    f.write(f"{name:39s} {val:.16e} +/- {err:.8e} {units}\n")
        print(f"Cross-section table: {path}")


def run_analysis(args):
    outdir = Path(args.outdir)
    outdir_vals = Path(args.outdir_vals)
    outdir.mkdir(parents=True, exist_ok=True)
    outdir_vals.mkdir(parents=True, exist_ok=True)
    basename = f"{args.savename}_{args.run}"
    logfile = outdir_vals / f"{basename}_pieces_and_skipped.log"
    tee = Tee(outdir_vals / f"{basename}_output.txt")
    sys.stdout = tee
    # Both "0" and "LO" refer to the same mp2mp0 piece.
    report = Report("LO" if name == "0" else name for name in args.pieces)
    failed = None
    try:
        run_folder = Path(args.rundir) / args.run / "out"
        print(f"Run: {basename}\nInput folder: {run_folder}")
        print(f"Declared pieces: {', '.join(args.pieces)}")
        if not run_folder.is_dir():
            raise FileNotFoundError(f"Missing McMule folder {run_folder}")
        setup(folder=str(run_folder))
        built = build_contributions(args, report)
        # Integrated CS do not depend on histogram bin compatibility; always
        # save these BEFORE plotting so an incompatible bin cannot hide them.
        save_cross_sections(outdir_vals / f"{basename}_cross_sections", built, report)

        muon_cfg = [
            dict(key="th3", label="th3", scale=1e-3,
                 xlabel=r"$\theta_\mu\ (\mathrm{mrad})$",
                 ylabel=r"$d\sigma/d\theta_\mu\ (\mu\mathrm{b}/\mathrm{mrad})$",
                 title="Muon Scattering Angle", xlim=(-0.25, 4.5), k_xmax=3.9),
            dict(key="Emu", label="Emu", scale=1e3,
                 xlabel=r"$E_\mu\ (\mathrm{GeV})$",
                 ylabel=r"$d\sigma/dE_\mu\ (\mu\mathrm{b}/\mathrm{GeV})$",
                 title="Scattered Muon Energy", show_k=False, delta_ymin=-100.0),
            dict(key="costh3", label="costh3", scale=1.0,
                 xlabel=r"$\cos\theta_\mu$",
                 ylabel=r"$d\sigma/d\cos\theta_\mu\ (\mu\mathrm{b})$",
                 title="Muon Scattering Angle Cosine", xlim=(0.98, 1.0)),
            dict(key="Qsq", label="Q2", scale=1e6,
                 xlabel=r"$Q^2\ (\mathrm{GeV}^2)$",
                 ylabel=r"$d\sigma/dQ^2\ (\mu\mathrm{b}/\mathrm{GeV}^2)$",
                 title="Momentum Transfer"),
        ]
        photon_cfg = [
            dict(key="th5", label="th5", scale=1e-3,
                 xlabel=r"$\theta_\gamma\ (\mathrm{mrad})$",
                 ylabel=r"$d\sigma_R/d\theta_\gamma\ (\mu\mathrm{b}/\mathrm{mrad})$",
                 title="Photon Polar Angle", xlim=(-1.0, 13.5)),
            dict(key="Eph", label="Eph", scale=1e3,
                 xlabel=r"$E_\gamma\ (\mathrm{GeV})$",
                 ylabel=r"$d\sigma_R/dE_\gamma\ (\mu\mathrm{b}/\mathrm{GeV})$",
                 title="Photon Energy"),
            dict(key="phi5", label="phi5", scale=1.0,
                 xlabel=r"$\phi_\gamma\ (\mathrm{rad})$",
                 ylabel=r"$d\sigma_R/d\phi_\gamma\ (\mu\mathrm{b}/\mathrm{rad})$",
                 title="Photon Azimuth"),
        ]
        root_queue = []
        root_names = {
            "LO": "lo", "NLO0": "nlo0", "R": "real_fks",
            "NLO_full": "nlo_full", "Delta_NLO": "delta_nlo",
            "Delta_NLO0": "delta_nlo0",
        }

        for cfg in muon_cfg:
            label = cfg["label"]
            records = {}
            for key, root_name in root_names.items():
                if key not in built:
                    report.skip(f"ROOT {root_name}_{label}",
                                "required McMule piece(s) not available")
            for key, record in built.items():
                h = get_hist(record, cfg, report, f"ROOT {root_names[key]}_{label}", key)
                if h is not None:
                    records[key] = h
                    queue_root(root_queue, f"{root_names[key]}_{label}", h, cfg)
                    path = outdir_vals / f"{root_names[key]}_{label}_{basename}_lab.txt"
                    write_file_with_values(path, h, f"{label} (McMule raw x units)", key)

            # Each plot is independent: incompatible grids must not abort the run.
            def plot_safely(path, action):
                try:
                    action()
                    report.ok(path)
                    return True
                except Exception as exc:
                    report.skip(path.name, f"{type(exc).__name__}: {exc}")
                    path.unlink(missing_ok=True)
                    return False

            pair_path = outdir / f"{basename}_lab_{label}_pair.png"
            delta_path = outdir / f"{basename}_lab_{label}_delta_nlo.png"
            k_path = outdir_vals / f"k_delta_over_nlo_{label}_{basename}_lab.txt"
            if all(k in records for k in ("LO", "NLO_full", "Delta_NLO")):
                def draw_main_pair():
                    k = pair_and_k(pair_path, records["LO"], records["NLO_full"],
                                   cfg, "NLO0 + R")
                    if k is not None:
                        write_file_with_values(k_path, k,
                                               f"{label} (display x units)",
                                               "K=(NLO-LO)/NLO")
                        queue_root(root_queue, f"k_delta_over_nlo_{label}",
                                   k, cfg, ratio=True)
                    else:
                        k_path.unlink(missing_ok=True)
                if not plot_safely(pair_path, draw_main_pair):
                    k_path.unlink(missing_ok=True)
                def draw_main_delta():
                    save_delta_nlo(
                        savename=f"{basename}_lab_{label}_delta_nlo",
                        delta_hist=records["Delta_NLO"],
                        title=cfg["title"], ymin=cfg.get("delta_ymin"),
                        scale_factor=cfg["scale"], x_label=cfg["xlabel"],
                        y_label=cfg["ylabel"], outdir=outdir,
                        xlim=cfg.get("xlim"),
                    )
                plot_safely(delta_path, draw_main_delta)
            else:
                missing = [x for x in ("LO", "NLO0", "R") if x not in built]
                reason = ("missing piece(s): " + ", ".join(missing)
                          if missing else "NLO0+R histogram has incompatible/missing bins")
                report.skip(pair_path.name, reason)
                report.skip(delta_path.name, reason)
                report.skip(f"k_delta_over_nlo_{label}", reason)
                for stale in (pair_path, delta_path, k_path):
                    stale.unlink(missing_ok=True)
                for key in ("LO", "NLO0", "R", "NLO_full"):
                    if key in records:
                        out = outdir / f"{basename}_lab_{label}_{key.lower()}.png"
                        plot_safely(out, lambda key=key, out=out: save_single(
                            out, records[key], cfg, key,
                            COLOR_LO if key == "LO" else COLOR_NLO))

            delta0_single = outdir / f"{basename}_lab_{label}_delta_nlo0.png"
            if "Delta_NLO0" in records:
                plot_safely(delta0_single, lambda: save_single(
                    delta0_single, records["Delta_NLO0"], cfg,
                    "Delta NLO0 = NLO0 - LO (F + A)", COLOR_DIFF,
                    title=f"{cfg['title']}: Delta NLO0 (F + A)"))
            else:
                report.skip(delta0_single.name, "requires NLO0 and LO with compatible bins")
                delta0_single.unlink(missing_ok=True)

        for cfg in photon_cfg:
            label = cfg["label"]
            name = f"{basename}_lab_{label}_real_photon"
            if "R" not in built:
                report.skip(name, "R not available")
                report.skip(f"ROOT real_fks_{label}", "R not available")
                continue
            hist = get_hist(built["R"], cfg, report, name, "R")
            if hist is None:
                continue
            save_real_photon(savename=name, real_hist=hist,
                             scale_factor=cfg["scale"], x_label=cfg["xlabel"],
                             y_label=cfg["ylabel"], title=cfg["title"],
                             outdir=outdir, xlim=cfg.get("xlim"))
            report.ok(name)
            queue_root(root_queue, f"real_fks_{label}", hist, cfg)
            write_file_with_values(
                outdir_vals / f"real_fks_{label}_{basename}_lab.txt",
                hist, f"{label} (McMule raw x units)", "R (FKS-subtracted)")

        if args.with_analytic:
            if "LO" in built:
                from theo_calc import dsigma_dQ2, dsigma_dcosth
                save_analytic_lo_check(
                    savename=f"{basename}_lab_Q2_analytic",
                    lo_hist=built["LO"]["Qsq"], scale_factor=1e6,
                    func=dsigma_dQ2, outdir=outdir,
                    x_label=r"$Q^2\ (\mathrm{GeV}^2)$",
                    y_label=r"$d\sigma/dQ^2\ (\mu\mathrm{b}/\mathrm{GeV}^2)$",
                    title=r"LO vs. Analytic $Q^2$")
                save_analytic_lo_check(
                    savename=f"{basename}_lab_costh3_analytic",
                    lo_hist=built["LO"]["costh3"], scale_factor=1.0,
                    func=lambda c: dsigma_dcosth(np.arccos(np.clip(c, -1.0, 1.0))),
                    outdir=outdir, x_label=r"$\cos\theta_\mu$",
                    y_label=r"$d\sigma/d\cos\theta_\mu\ (\mu\mathrm{b})$",
                    title=r"LO vs. Analytic $\cos\theta_\mu$",
                    yscale="linear", xlim=(0.98, 1.0))
                report.ok("analytic LO Q2 and costh3")
            else:
                report.skip("analytic LO Q2 and costh3", "LO not available")

        root_path = Path(args.root_file) if args.root_file else outdir / f"{basename}_histograms.root"
        try:
            save_root(root_path, root_queue, report)
        except Exception as exc:
            report.skip(root_path.name, f"ROOT export: {type(exc).__name__}: {exc}")
        if args.split_root:
            for suffix, prefix in (("nlo0", "nlo0_"), ("real_fks", "real_fks_")):
                subset = [item for item in root_queue if item["name"].startswith(prefix)]
                if not subset:
                    report.skip(f"ROOT {suffix}", "no histograms available")
                    continue
                split_path = root_path.with_name(f"{root_path.stem}_{suffix}.root")
                try:
                    save_root(split_path, subset, report)
                except Exception as exc:
                    report.skip(split_path.name, f"ROOT export: {type(exc).__name__}: {exc}")
        report.note("No adaptive rebinning was performed.")
        if "NLO_full" not in built:
            report.note("NLO0+R combined record unavailable; integrated CS were saved "
                        "independently when all input pieces existed.")
    except Exception as exc:
        failed = traceback.format_exc()
        print(f"ERROR: {type(exc).__name__}: {exc}")
        raise
    finally:
        try:
            report.write(logfile, error=failed)
        finally:
            sys.stdout = tee.terminal
            tee.close()


def main():
    parser = argparse.ArgumentParser(description="McMule LO, NLO0 (including F+A), R and full NLO")
    parser.add_argument("--momosdir", required=True,
                        help="Retained for compatibility with your command file")
    parser.add_argument("--outdir", required=True)
    parser.add_argument("--outdir-vals", required=True)
    parser.add_argument("--run", required=True)
    parser.add_argument("--rundir", required=True)
    parser.add_argument("--savename", required=True)
    parser.add_argument("--nbins", type=int, required=True,
                        help="Retained for compatibility; original bins are not rebinned")
    parser.add_argument("--pieces", nargs="+", required=True,
                        choices=(*PIECE_NAMES, "0"),
                        help="Pieces from .conf, e.g. --pieces 0 NLO0 R; LO is an alias for 0")
    parser.add_argument("--root-file", default=None)
    parser.add_argument("--split-root", action="store_true",
                        help="Also export separate NLO0 and R histogram ROOT files (not event trees)")
    parser.add_argument("--with-analytic", action="store_true")
    run_analysis(parser.parse_args())


if __name__ == "__main__":
    main()