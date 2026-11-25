#!/usr/bin/env python3
"""Separate McMule-style PNGs of TGEANT histograms and a summary text file.

No histogram is rebinned or normalized. The ROOT archive also contains a
TGEANTTruth event tree with parsed primary-particle info, NOT detector hits.
"""
import argparse
import re
import math
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import matplotlib.colors as mcolors
import numpy as np

plt.rcParams.update({
    "font.family": "serif", "font.size": 12, "text.usetex": False,
    "axes.grid": True, "grid.alpha": 0.28, "grid.linestyle": ":",
    "savefig.dpi": 300,
})
COLOR = "#2ca02c"   # same green used for NLO in McMule plots

# Plot titles and axes match the conventions used for McMule ROOT event plots.
PLOT_INFO = {
    "h_muon_energy": ("Scattered Muon Energy", r"$E_\mu^{\mathrm{lab}}\ (\mathrm{GeV})$", "Events / bin"),
    "h_muon_theta": ("Muon Scattering Angle", r"$\theta_\mu^{\mathrm{lab}}\ (\mathrm{mrad})$", "Events / bin"),
    "h_Q2": ("Momentum Transfer", r"$Q^2\ (\mathrm{GeV}^2)$", "Events / bin"),
    "h_muon_energy_loss": ("Muon Energy Loss", r"$\Delta E_\mu\ (\mathrm{GeV})$", "Events / bin"),
    "h_photon_energy": ("Real Photon Energy", r"$E_\gamma^{\mathrm{lab}}\ (\mathrm{GeV})$", "Photons / bin"),
    "h_photon_theta": ("Real Photon Polar Angle", r"$\theta_\gamma^{\mathrm{lab}}\ (\mathrm{mrad})$", "Photons / bin"),
    "h_photon_theta_full": ("Real Photon Polar Angle (Full Range)", r"$\theta_\gamma^{\mathrm{lab}}\ (\mathrm{mrad})$", "Photons / bin"),
}



def clean_text(s):
    """Translate simple ROOT TLatex tokens to readable Unicode for mpl."""
    s = str(s or "")
    for key, value in (
        ("#theta", "θ"), ("#phi", "φ"), ("#mu", "μ"),
        ("#gamma", "γ"), ("#Delta", "Δ"), ("#sigma", "σ"),
        ("#pi", "π"), ("#eta", "η"), ("#sqrt", "√"),
    ):
        s = s.replace(key, value)
    s = s.replace("^{2}", "²").replace("^{0}", "⁰")
    return s


def filename_for(key):
    safe = re.sub(r"[^A-Za-z0-9_.-]+", "_", key.replace("/", "__"))
    return safe[:180].strip("._") or "hist"


def obj_walk(directory, prefix=""):
    """Yield the latest cycle for each named ROOT object, recursively."""
    keys = {}
    for k in directory.GetListOfKeys():
        name = str(k.GetName())
        if name not in keys or k.GetCycle() > keys[name].GetCycle():
            keys[name] = k
    for name, k in sorted(keys.items()):
        obj = k.ReadObj()
        fullname = f"{prefix}/{name}" if prefix else name
        if obj.InheritsFrom("TDirectory"):
            yield from obj_walk(obj, fullname)
        else:
            yield fullname, obj


def log_or_linear(ax, y):
    valid = y[np.isfinite(y)]
    if len(valid) == 0:
        return
    pos = valid[valid > 0]
    if not len(pos):
        return
    dynamic_range = float(np.max(pos) / np.min(pos))
    if np.any(valid < 0):
        ax.set_yscale("symlog", linthresh=max(float(np.max(np.abs(valid))) * 0.001, 1.0))
    elif dynamic_range > 1000:
        # symlog keeps zero bins visible and does not invent 1e-300 values
        ax.set_yscale("symlog", linthresh=max(float(np.max(pos)) * 0.001, 1.0))


def draw_th1(h, name, folder):
    nx = h.GetNbinsX()
    if nx <= 0:
        return
    edges = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)], dtype=float)
    y = np.array([h.GetBinContent(i) for i in range(1, nx + 1)], dtype=float)
    errors = np.array([h.GetBinError(i) for i in range(1, nx + 1)], dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 4.5))
    ax.stairs(y, edges, color=COLOR, linewidth=1.5, label="NLO (TGEANT)")
    if nx <= 1200 and np.any(errors > 0):
        y_low = y - errors
        y_high = y + errors
        if np.all(y >= 0):
            y_low = np.maximum(y_low, 0.0)
        # Use a stepped error band; unlike ax.stairs(fill=True), this does
        # not draw spurious rectangles stretching through the zero line.
        ax.fill_between(edges, np.r_[y_low, y_low[-1]],
                        np.r_[y_high, y_high[-1]], step="post",
                        color=COLOR, alpha=0.16, linewidth=0)
    log_or_linear(ax, np.concatenate((y, y-errors, y+errors)))
    title = clean_text(h.GetTitle()) or name.split("/")[-1].replace("_", " ")
    # Some TH1 titles contain inline ;x;y ROOT syntax, not user facing.
    if ";" in title:
        title = title.split(";", 1)[0]
    title, xlabel, ylabel = PLOT_INFO.get(
        name,
        (title, clean_text(h.GetXaxis().GetTitle()) or name.split("/")[-1].replace("_", " "),
         clean_text(h.GetYaxis().GetTitle()) or "Bin content"),
    )
    ax.set_title(title)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xlim(edges[0], edges[-1])
    ax.legend(frameon=False)
    fig.savefig(folder / f"{filename_for(name)}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def draw_th2(h, name, folder):
    nx, ny = h.GetNbinsX(), h.GetNbinsY()
    if nx * ny > 4_000_000:
        print(f"SKIP huge TH2 {name}: {nx}x{ny}")
        return
    x = np.array([h.GetXaxis().GetBinLowEdge(i) for i in range(1, nx + 2)], dtype=float)
    y = np.array([h.GetYaxis().GetBinLowEdge(i) for i in range(1, ny + 2)], dtype=float)
    z = np.array([[h.GetBinContent(ix, iy) for ix in range(1, nx + 1)]
                  for iy in range(1, ny + 1)], dtype=float)
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    # Do not force a log colorbar: allow zeros and signed bins.
    mesh = ax.pcolormesh(x, y, z, shading="flat", cmap="viridis", rasterized=True)
    fig.colorbar(mesh, ax=ax, label="Bin content")
    title = clean_text(h.GetTitle()) or name.split("/")[-1].replace("_", " ")
    if ";" in title:
        title = title.split(";", 1)[0]
    ax.set_title(title)
    ax.set_xlabel(clean_text(h.GetXaxis().GetTitle()))
    ax.set_ylabel(clean_text(h.GetYaxis().GetTitle()))
    fig.savefig(folder / f"{filename_for(name)}.png", dpi=300, bbox_inches="tight")
    plt.close(fig)


def main():
    p = argparse.ArgumentParser()
    p.add_argument("input_root")
    p.add_argument("output_directory")
    args = p.parse_args()
    outdir = Path(args.output_directory)
    outdir.mkdir(parents=True, exist_ok=True)
    import ROOT
    ROOT.gROOT.SetBatch(True)
    f = ROOT.TFile.Open(args.input_root, "READ")
    if not f or f.IsZombie():
        raise RuntimeError(f"Cannot open ROOT: {args.input_root}")
    inventory = ["ROOT object inventory", f"File: {args.input_root}", "",
                 "NOTE: Histogram entries do not necessarily equal simulated event counts.",
                 "NOTE: This inventory does not prove that raw TGEANT data are preserved.", ""]
    nplots = 0
    counter_lines = ["TGEANT analysis totals (all merged batches)",
                     f"ROOT file: {args.input_root}", "",
                     "Entries count parsed primary events, not detector hits.",
                     ""]
    for name, obj in obj_walk(f):
        if obj.InheritsFrom("TH2"):
            inventory.append(f"TH2\t{name}\tentries={obj.GetEntries():.0f}\t"
                             f"content_sum={obj.Integral(0,obj.GetNbinsX()+1,0,obj.GetNbinsY()+1):.8g}")
            draw_th2(obj, name, outdir)
            nplots += 1
        elif obj.InheritsFrom("TH1"):
            inventory.append(f"TH1\t{name}\tentries={obj.GetEntries():.0f}\t"
                             f"content_sum={obj.Integral(0,obj.GetNbinsX()+1):.8g}")
            if name == "h_event_summary":
                for i in range(1, obj.GetNbinsX()+1):
                    k = str(obj.GetXaxis().GetBinLabel(i)) or f"category_{i}"
                    counter_lines.append(f"{k}: {obj.GetBinContent(i):.0f}")
            else:
                draw_th1(obj, name, outdir)
                nplots += 1
        elif obj.InheritsFrom("TTree"):
            inventory.append(f"TREE\t{name}\tentries={obj.GetEntries()}\t"
                             f"branches={','.join(b.GetName() for b in obj.GetListOfBranches())}")
            counter_lines.append(f"{name}_entries: {obj.GetEntries()}")
        else:
            inventory.append(f"{obj.ClassName()}\t{name}")
    inventory += ["", f"Separate plots created: {nplots}"]
    (outdir / "root_object_inventory.txt").write_text("\n".join(inventory) + "\n", encoding="utf-8")
    counter_lines.extend(["", "WARNING: TGEANTTruth does NOT contain detector hits or",
                          "reconstructable raw TGEANT track data. Keep .tgeant.gz until",
                          "after CORAL/PHAST, full validation, and an explicit archive decision."])
    (outdir / "event_statistics.txt").write_text("\n".join(counter_lines) + "\n", encoding="utf-8")
    f.Close()
    print(f"Created {nplots} individual figures under {outdir}")


if __name__ == "__main__":
    main()
