from pymule import *

import matplotlib.pyplot as plt
import numpy as np

from plotting import *
from theo_calc import *

from matplotlib.colors import LogNorm

import argparse
import os
import sys


class Tee:
    """Writes output to both terminal and a file."""

    def __init__(self, filename):
        self.terminal = sys.stdout
        self.file = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


# =========================
# Command line input
# =========================

parser = argparse.ArgumentParser()

parser.add_argument("--homedir", required=True)
parser.add_argument("--outdir", required=True)
parser.add_argument("--outdir-vals", required=True)
parser.add_argument("--run", required=True)
parser.add_argument("--savename", required=True)
parser.add_argument("--nbins", type=int, required=True)

args = parser.parse_args()


homedir = args.homedir
outdir = args.outdir.rstrip("/") + "/"
outdir_vals = args.outdir_vals.rstrip("/") + "/"

run = args.run
savename = args.savename
nbins = args.nbins


savename_base = f"{savename}_{run}"


# =========================
# Redirect stdout
# =========================

log_file = outdir_vals + f"{savename_base}_output.txt"

tee = Tee(log_file)
sys.stdout = tee

print("=======================================================================")
print("        Analysed file: ", savename_base)
print("=======================================================================")


# =========================
# Physics setup
# =========================

run_folder = os.path.join(homedir, run, "out")

setup(folder=run_folder)

lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv
nlo = mergefks(sigma("mp2mpR"),sigma("mp2mpF"),anyxi=sigma("mp2mpA")) * alpha**3 * conv
full = lo + nlo
onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

# =========================
# Extract observables (LAB)
# =========================
lo_th3,   nlo_th3,   full_th3   = lo["th3"],   nlo["th3"],   full["th3"]
lo_Emu,   nlo_Emu,   full_Emu   = lo["Emu"],   nlo["Emu"],   full["Emu"]
lo_th5,   nlo_th5,   full_th5   = lo["th5"],   nlo["th5"],   full["th5"]
lo_Eph,   nlo_Eph,   full_Eph   = lo["Eph"],   nlo["Eph"],   full["Eph"]
lo_phi5,  nlo_phi5,  full_phi5  = lo["phi5"],  nlo["phi5"],  full["phi5"]
lo_costh3,nlo_costh3,full_costh3= lo["costh3"],nlo["costh3"],full["costh3"]
lo_Q2,    nlo_Q2,    full_Q2    = lo["Qsq"],   nlo["Qsq"],   full["Qsq"]

# =========================
# Colors
# =========================
colors = dict(lo="#1f77b4", nlo="#ff7f0e", full="#2ca02c", K="#d62728")

# =========================
# Limits
# =========================
th3vals = finite_bins(full_th3)
th3vals = th3vals[th3vals[:, 1] != 0]
th3_min = np.min(th3vals[:, 0])
th3_max = np.max(th3vals[:, 0])
print("th3_min:    ", th3_min * 1e3, ", th3_max:   ", th3_max * 1e3, "(mrad)")

costh3vals = finite_bins(full_costh3)
costh3vals = costh3vals[costh3vals[:, 1] != 0]
costh3_min = np.min(costh3vals[:, 0])
costh3_max = np.max(costh3vals[:, 0])
print("costh3_min: ", costh3_min, ", costh3_max: ", costh3_max)


# =========================
# Function to make plots & K-factors
# =========================
def make_plots_and_kfactors(*, tag, savename_base,
                             lo_th3, nlo_th3, full_th3,
                             lo_Emu, nlo_Emu, full_Emu,
                             lo_th5, nlo_th5, full_th5,
                             lo_Eph, nlo_Eph, full_Eph,
                             lo_phi5, nlo_phi5, full_phi5,
                             lo_costh3, nlo_costh3, full_costh3,
                             lo_Q2, nlo_Q2, full_Q2,
                             outdir, outdir_vals, colors):
    savename = f"{savename_base}_{tag}"

    # =========================
    # Write value files
    # =========================
    variables = [
        ("th3",   False, False),
        ("Emu",   False, False),
        ("th5",   True,  False),
        ("Eph",   True,  False),
        ("phi5",  True,  False),
        ("costh3",False, True),
        ("Q2",    False, True),
    ]

    data_map = {
        "th3":    {"lo": lo_th3,    "nlo": nlo_th3,    "full": full_th3},
        "Emu":    {"lo": lo_Emu,    "nlo": nlo_Emu,    "full": full_Emu},
        "th5":    {"lo": lo_th5,    "nlo": nlo_th5,    "full": full_th5},
        "Eph":    {"lo": lo_Eph,    "nlo": nlo_Eph,    "full": full_Eph},
        "phi5":   {"lo": lo_phi5,   "nlo": nlo_phi5,   "full": full_phi5},
        "costh3": {"lo": lo_costh3, "nlo": nlo_costh3, "full": full_costh3},
        "Q2":     {"lo": lo_Q2,     "nlo": nlo_Q2,     "full": full_Q2},
    }

    for var, photon_only, lab_only in variables:
        if lab_only and tag != "lab":
            continue
        orders = ["nlo"] if photon_only else ["lo", "nlo", "full"]
        for order in orders:
            arr = data_map.get(var, {}).get(order, None)
            if arr is not None:
                write_file_with_values(outdir_vals + f"{order}_{var}_{savename}.txt", arr,
                                       f"{var}_{tag} bin centers", "value")

    # =========================
    # Combined figure
    # =========================
    nrows = 6
    height_ratios = [3, 1] * (nrows // 2)

    fig, axes = create_figure(nrows=nrows, ncols=2, figsize=(16, 6 + nrows * 2), font_size=12,
                              sharex=False, gridspec_kw={"height_ratios": height_ratios, "hspace": 0.6})

    ax_th3,  ax_Q2   = axes[0]
    ax_K_th3,ax_K_Q2 = axes[1]
    ax_Emu,  ax_Eph  = axes[2]
    ax_K_Emu,ax_K_Eph= axes[3]
    ax_th5,  ax_phi5 = axes[4]
    ax_K_th5,ax_K_phi5 = axes[5]

    ax_K_th3.sharex(ax_th3)
    ax_K_Emu.sharex(ax_Emu)
    ax_K_th5.sharex(ax_th5)
    ax_K_Eph.sharex(ax_Eph)
    ax_K_phi5.sharex(ax_phi5)
    ax_K_Q2.sharex(ax_Q2)

    draw_observable_and_k(ax_th3, ax_K_th3,
        lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3,
        scale_factor=1e-3, x_label_main=r"$\theta_3$ (mrad)", x_label_k=r"$\theta_3$ (mrad)",
        y_label_main=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Muon Scattering Angle ({tag})", force_main_linear=False, colors=colors)

    draw_observable_and_k(ax_Emu, ax_K_Emu,
        lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu,
        scale_factor=1e3, x_label_main=r"$E_\mu$ (GeV)", x_label_k=r"$E_\mu$ (GeV)",
        y_label_main=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        main_title=f"Muon Energy ({tag})", colors=colors)

    draw_observable_and_k(ax_th5, ax_K_th5,
        lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5,
        scale_factor=1e-3, x_label_main=r"$\theta_5$ (mrad)", x_label_k=r"$\theta_5$ (mrad)",
        y_label_main=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Scattering Angle ({tag})", xlim=(-1., 13.), colors=colors)

    draw_observable_and_k(ax_Eph, ax_K_Eph,
        lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph,
        scale_factor=1e3, x_label_main=r"$E_\gamma$ (GeV)", x_label_k=r"$E_\gamma$ (GeV)",
        y_label_main=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        main_title=f"Photon Energy ({tag})", colors=colors)

    draw_observable_and_k(ax_phi5, ax_K_phi5,
        lo_hist=lo_phi5, nlo_hist=nlo_phi5, full_hist=full_phi5,
        scale_factor=1e-3, x_label_main=r"$\phi_5$ (mrad)", x_label_k=r"$\phi_5$ (mrad)",
        y_label_main=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Deflection Angle ({tag})", force_main_linear=False, colors=colors)

    draw_observable_and_k(ax_Q2, ax_K_Q2,
        lo_hist=lo_Q2, nlo_hist=nlo_Q2, full_hist=full_Q2,
        scale_factor=1.e6, x_label_main=r"$Q^2$", x_label_k=r"$Q^2$",
        y_label_main=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
        main_title=f"Momentum Transfer ({tag})", force_main_linear=True, colors=colors)

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)

    # =========================
    # Separate pair plots
    # =========================
    save_single_pair_plot(savename=f"{savename}_th3_pair",
        lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3,
        scale_factor=1.e-3, x_label=r"$\theta_3\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Muon Scattering Angle ({tag})", force_main_linear=False, colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_Emu_pair",
        lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu,
        scale_factor=1.e3, x_label=r"$E_\mu\ (\mathrm{GeV})$",
        y_label=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        main_title=f"Energy of the Scattered Muon ({tag})", colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_th5_pair",
        lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5,
        scale_factor=1.e-3, x_label=r"$\theta_5\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Scattering Angle ({tag})", xlim=(-1., 13.), colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_Eph_pair",
        lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph,
        scale_factor=1.e3, x_label=r"$E_\gamma\ (\mathrm{GeV})$",
        y_label=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        main_title=f"Photon Energy ({tag})", colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_phi5_pair",
        lo_hist=lo_phi5, nlo_hist=nlo_phi5, full_hist=full_phi5,
        scale_factor=1.e-3, x_label=r"$\phi_5\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon X-deflection ({tag})", colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_costh3_pair",
        lo_hist=lo_costh3, nlo_hist=nlo_costh3, full_hist=full_costh3,
        scale_factor=1., x_label=r"$\cos\theta_3$",
        y_label=r"$\frac{d\sigma}{d\cos\theta_3}\ (\mu\mathrm{barn})$",
        main_title=f"Muon Scattering Angle cos ({tag})", xlim=(0.98, 1.),
        force_main_linear=True, colors=colors, outdir=outdir)

    save_single_pair_plot(savename=f"{savename}_Q2_pair",
        lo_hist=lo_Q2, nlo_hist=nlo_Q2, full_hist=full_Q2,
        scale_factor=1.e6, x_label=r"$Q^2\ (\mathrm{GeV}^2)$",
        y_label=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
        main_title=f"$Q^2$ ({tag})", force_main_linear=False, colors=colors, outdir=outdir)


# =========================
# Run for LAB
# =========================
make_plots_and_kfactors(tag="lab", savename_base=savename_base,
                        lo_th3=lo_th3, nlo_th3=nlo_th3, full_th3=full_th3,
                        lo_Emu=lo_Emu, nlo_Emu=nlo_Emu, full_Emu=full_Emu,
                        lo_th5=lo_th5, nlo_th5=nlo_th5, full_th5=full_th5,
                        lo_Eph=lo_Eph, nlo_Eph=nlo_Eph, full_Eph=full_Eph,
                        lo_phi5=lo_phi5, nlo_phi5=nlo_phi5, full_phi5=full_phi5,
                        lo_costh3=lo_costh3, nlo_costh3=nlo_costh3, full_costh3=full_costh3,
                        lo_Q2=lo_Q2, nlo_Q2=nlo_Q2, full_Q2=full_Q2,
                        outdir=outdir, outdir_vals=outdir_vals, colors=colors)


# =========================
# Plot Numeric vs. Analytic
# =========================
plot_costh3_with_analytic(lo_hist=lo_costh3,
                          nlo_hist=nlo_costh3,
                          full_hist=full_costh3,
                          nbins=nbins,
                          th3_min=th3_min, th3_max=th3_max,
                          colors=colors,
                          savename=f"{savename_base}_costh3_analytic",
                          outdir=outdir,
                          outdir_vals=outdir_vals)

plot_Q2_with_analytic(lo_hist=lo_Q2,
                      nlo_hist=nlo_Q2,
                      full_hist=full_Q2,
                      nbins=nbins,
                      ylow_diff=-0.5, yup_diff=1.5,
                      colors=colors,
                      savename=f"{savename_base}_Q2_analytic",
                      outdir=outdir,
                      outdir_vals=outdir_vals)


# =========================
# Calculate total cross section
# =========================
sigma_lo_mb   = lo.value[0]    / 1000
sigma_nlo_mb  = nlo.value[0]   / 1000
sigma_full_mb = full.value[0]  / 1000
Rate     = calculate_rate(sigma_lo_mb)
Rate_nlo  = calculate_rate(sigma_nlo_mb)

print("\n-------------------------------------------------------------")
print("LO cross section:            ", sigma_lo_mb, "mb")
print("LO cross section Paper (Big): 0.255 mb")
print("NLO cross section:           ", sigma_nlo_mb, "mb")
print("Full cross section:          ", sigma_full_mb, "mb")
print("LO Rate:                     ", Rate, "1/s")
print("NLO Rate:                    ", Rate_nlo, "1/s")
print("\n-------------------------------------------------------------")


sys.stdout = tee.terminal
tee.close()
