from pymule import *
import matplotlib.pyplot as plt
import numpy as np

from plotting import *

# =========================
# Paths
# =========================
homedir = "/home/marialei/AMBER_RadCor/"
outdir = homedir + "Figures/"

# =========================
# Input definitions
# =========================
lo_outs = [
    "mp2mp_NLO_22_12",
    "mp2mp_NLO_12_01",
    "mp2mp_NLO_13_01",
    "mp2mp_NLO_19_01",
    "mp2mp_NLO_01_02"
]

nlo_outs = lo_outs
savenames = ["combined", "combined_16Jan", "combined_19Jan", "combined_20_Jan"]

# =========================
# Binning (rad, MeV)
# =========================
nrbins = 500
min_val = np.array([0.3e-3, 95e3, -12e-3, 1e3])
max_val = np.array([2.0e-3, 101e3, 12e-3, 100e3])

binwidths_th3 = np.full(nrbins, (max_val[0] - min_val[0]) / nrbins)
binwidths_Emu = np.full(nrbins, (max_val[1] - min_val[1]) / nrbins)
binwidths_th5 = np.full(nrbins, (max_val[2] - min_val[2]) / nrbins)
binwidths_Eph = np.full(nrbins, (max_val[3] - min_val[3]) / nrbins)

# =========================
# Dataset choice
# =========================
lo_i = 4
nlo_i = 4
savename = savenames[0] + "_" + nlo_outs[nlo_i]

# =========================
# Physics setup
# =========================
setup(folder=homedir + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

setup(folder=homedir + nlo_outs[nlo_i] + "/out")
nlo = (
    mergefks(
        sigma("mp2mpR"),
        sigma("mp2mpF"),
        anyxi=sigma("mp2mpA"),
    )
    * alpha**3
    * conv
)

full = lo + nlo

# =========================
# Extract observables
# =========================
lo_th3, nlo_th3, full_th3 = lo["th3"], nlo["th3"], full["th3"]
lo_Emu, nlo_Emu, full_Emu = lo["Emu"], nlo["Emu"], full["Emu"]
lo_th5, nlo_th5, full_th5 = lo["th5"], nlo["th5"], full["th5"]
lo_Eph, nlo_Eph, full_Eph = lo["Eph"], nlo["Eph"], full["Eph"]

# =========================
# Figure setup (2x4)
# =========================
fig, axes = create_figure(
    nrows=2,
    ncols=4,
    figsize=(24, 8),
    font_size=12,
)

ax_th3, ax_K_th3, ax_Emu, ax_K_Emu = axes[0]
ax_th5, ax_K_th5, ax_Eph, ax_K_Eph = axes[1]

# =========================
# Colors
# =========================
colors = dict(
    lo="#1f77b4",
    nlo="#ff7f0e",
    full="#2ca02c",
    K="#d62728",
)

# =========================
# Observable plots
# =========================
plot_lo_nlo_full(ax_th3, lo_th3, nlo_th3, full_th3, binwidths_th3, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_th3, r"$\theta_3\ (\mathrm{rad})$", "?", "Muon Scattering Angle")
ax_th3.set_xlim(1.3e-3, 1.7e-3)

plot_lo_nlo_full(ax_Emu, lo_Emu, nlo_Emu, full_Emu, binwidths_Emu, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_Emu, r"$E_\mu\ (\mathrm{MeV})$", "?", "Energy of the Scattered Muon")

plot_lo_nlo_full(ax_th5, lo_th5, nlo_th5, full_th5, binwidths_th5, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_th5, r"$\theta_5\ (\mathrm{rad})$", "?", "Photon Scattering Angle")
ax_th5.set_xlim(-2.e-3, 2.e-3)

plot_lo_nlo_full(ax_Eph, lo_Eph, nlo_Eph, full_Eph, binwidths_Eph, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_Eph, r"$E_\gamma\ (\mathrm{MeV})$", "?", "Photon Energy")

ax_th3.set_yscale("linear")
#ax_Emu.set_yscale("linear")
#ax_th5.set_yscale("linear")
#ax_Eph.set_yscale("linear")

# =========================
# K-factors (rebinned)
# =========================
K_th3 = mergebins(divideplots(nlo_th3, full_th3), 5)
K_Emu = mergebins(divideplots(nlo_Emu, full_Emu), 5)
K_th5 = mergebins(divideplots(nlo_th5, full_th5), 5)
K_Eph = mergebins(divideplots(nlo_Eph, full_Eph), 5)

bwK_th3 = np.ones_like(K_th3)
bwK_Emu = np.ones_like(K_Emu)
bwK_th5 = np.ones_like(K_th5)
bwK_Eph = np.ones_like(K_Eph)


plot_K(ax_K_th3, K_th3, bwK_th3, colors["K"], r"$\theta_3\ (\mathrm{rad})$", r"$K(\theta_3)$")
plot_K(ax_K_Emu, K_Emu, bwK_Emu, colors["K"], r"$E_\mu\ (\mathrm{MeV})$", r"$K(E_\mu)$")
plot_K(ax_K_th5, K_th5, bwK_th5, colors["K"], r"$\theta_5\ (\mathrm{rad})$", r"$K(\theta_5)$")
plot_K(ax_K_Eph, K_Eph, bwK_Eph, colors["K"], r"$E_\gamma\ (\mathrm{MeV})$", r"$K(E_\gamma)$")

ax_K_th3.set_xlim(ax_th3.get_xlim())
ax_K_Emu.set_xlim(ax_Emu.get_xlim())
ax_K_th5.set_xlim(ax_th5.get_xlim())
ax_K_Eph.set_xlim(ax_Eph.get_xlim())

ax_K_th3.set_yscale("linear")
#ax_K_Emu.set_yscale("linear")
#ax_K_th5.set_yscale("linear")
#ax_K_Eph.set_yscale("linear")


# =========================
# Layout & save
# =========================
fig.subplots_adjust(hspace=0.35, wspace=0.3)

save_figure(fig, savename, outdir=outdir)
