from pymule import *
import matplotlib.pyplot as plt
import numpy as np

from plotting import *

# =========================
# Paths
# =========================
homedir = "/home/marialei/AMBER_RadCor/"
outdir = homedir + "Figures/"
outdir_vals = homedir + "Vals/"

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


write_file_with_values(outdir_vals + "lo_th3_" + savename + ".txt", lo_th3, "th3 bin centers", "value")
write_file_with_values(outdir_vals + "nlo_th3_" + savename + ".txt", nlo_th3, "th3 bin centers", "value")
write_file_with_values(outdir_vals + "full_th3_" + savename + ".txt", full_th3, "th3 bin centers", "value")
print("[INFO] th3 files written.")
write_file_with_values(outdir_vals + "lo_Emu_" + savename + ".txt", lo_Emu, "Emu bin centers", "value")
write_file_with_values(outdir_vals + "nlo_Emu_" + savename + ".txt", nlo_Emu, "Emu bin centers", "value")
write_file_with_values(outdir_vals + "full_Emu_" + savename + ".txt", full_Emu, "Emu bin centers", "value")
print("[INFO] Emu files written.")
write_file_with_values(outdir_vals + "lo_th5_" + savename + ".txt", lo_th5, "th5 bin centers", "value")
write_file_with_values(outdir_vals + "nlo_th5_" + savename + ".txt", nlo_th5, "th5 bin centers", "value")
write_file_with_values(outdir_vals + "full_th5_" + savename + ".txt", full_th5, "th5 bin centers", "value")
print("[INFO] th5 files written.")
write_file_with_values(outdir_vals + "lo_Eph_" + savename + ".txt", lo_Eph, "Eph bin centers", "value")
write_file_with_values(outdir_vals + "nlo_Eph_" + savename + ".txt", nlo_Eph, "Eph bin centers", "value")
write_file_with_values(outdir_vals + "full_Eph_" + savename + ".txt", full_Eph, "Eph bin centers", "value")
print("[INFO] Eph files written.")


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
lo_th3_scaled = scaleplot(lo_th3, 1.e-3)
nlo_th3_scaled = scaleplot(nlo_th3, 1.e-3)
full_th3_scaled = scaleplot(full_th3, 1.e-3)
plot_lo_nlo_full(ax_th3, lo_th3_scaled, nlo_th3_scaled, full_th3_scaled, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_th3, r"$\theta_3\ (\mathrm{mrad})$", "?", "Muon Scattering Angle")
ax_th3.set_xlim(1.3, 1.7)

lo_Emu_scaled = scaleplot(lo_Emu, 1.e3)
nlo_Emu_scaled = scaleplot(nlo_Emu, 1.e3)
full_Emu_scaled = scaleplot(full_Emu, 1.e3)
plot_lo_nlo_full(ax_Emu, lo_Emu_scaled, nlo_Emu_scaled, full_Emu_scaled, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_Emu, r"$E_\mu\ (\mathrm{GeV})$", "?", "Energy of the Scattered Muon")

lo_th5_scaled = scaleplot(lo_th5, 1.e-3)
nlo_th5_scaled = scaleplot(nlo_th5, 1.e-3)
full_th5_scaled = scaleplot(full_th5, 1.e-3)
plot_lo_nlo_full(ax_th5, lo_th5_scaled, nlo_th5_scaled, full_th5_scaled, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_th5, r"$\theta_5\ (\mathrm{mrad})$", "?", "Photon Scattering Angle")
ax_th5.set_xlim(-2., 2.)

lo_Eph_scaled = scaleplot(lo_Eph, 1.e3)
nlo_Eph_scaled = scaleplot(nlo_Eph, 1.e3)
full_Eph_scaled = scaleplot(full_Eph, 1.e3)
plot_lo_nlo_full(ax_Eph, lo_Eph_scaled, nlo_Eph_scaled, full_Eph_scaled, colors, labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
style_sci_x(ax_Eph, r"$E_\gamma\ (\mathrm{GeV})$", "?", "Photon Energy")

ax_th3.set_yscale("linear")
#ax_Emu.set_yscale("linear")
#ax_th5.set_yscale("linear")
#ax_Eph.set_yscale("linear")


# =========================
# K-factors (rebinned)
# =========================
K_th3 = mergebins(divideplots(nlo_th3_scaled, full_th3_scaled), 5)
K_Emu = mergebins(divideplots(nlo_Emu_scaled, full_Emu_scaled), 5)
K_th5 = mergebins(divideplots(nlo_th5_scaled, full_th5_scaled), 5)
K_Eph = mergebins(divideplots(nlo_Eph_scaled, full_Eph_scaled), 5)

write_file_with_values(outdir_vals + "K_theta_3_" + savename + ".txt", K_th3,"th3 bin center", "K_th3")
write_file_with_values(outdir_vals + "K_E_mu_" + savename + ".txt", K_Emu, "Emu bin center", "K_Emu")
write_file_with_values(outdir_vals + "K_theta_5_" + savename + ".txt", K_th5, "th5 bin center", "K_th5")
write_file_with_values(outdir_vals + "K_E_gamma_" + savename + ".txt", K_Eph, "Eph bin center", "K_Eph")
print("[INFO] K-factor files written.")

plot_K(ax_K_th3, K_th3, colors["K"], r"$\theta_3\ (\mathrm{mrad})$", r"$K(\theta_3)$")
plot_K(ax_K_Emu, K_Emu, colors["K"], r"$E_\mu\ (\mathrm{GeV})$", r"$K(E_\mu)$")
plot_K(ax_K_th5, K_th5, colors["K"], r"$\theta_5\ (\mathrm{mrad})$", r"$K(\theta_5)$")
plot_K(ax_K_Eph, K_Eph, colors["K"], r"$E_\gamma\ (\mathrm{GeV})$", r"$K(E_\gamma)$")

ax_K_th3.set_xlim(ax_th3.get_xlim())
ax_K_Emu.set_xlim(ax_Emu.get_xlim())
ax_K_th5.set_xlim(ax_th5.get_xlim())
ax_K_Eph.set_xlim(ax_Eph.get_xlim())

ax_K_th3.set_yscale("linear")
ax_K_Emu.set_yscale("linear")
ax_K_th5.set_yscale("linear")
#ax_K_Eph.set_yscale("linear")


# =========================
# Layout & save
# =========================
fig.subplots_adjust(hspace=0.35, wspace=0.3)

save_figure(fig, savename, outdir=outdir)
