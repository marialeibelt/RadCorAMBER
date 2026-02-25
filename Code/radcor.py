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
    "mp2mp_NLO_01_02",
    "mp2mp_NLO_24_02"
]
nlo_outs = lo_outs
savenames = ["combined", "combined_16Jan", "combined_19Jan", "combined_20_Jan"]

# =========================
# Dataset choice
# =========================
lo_i = 5
nlo_i = 5
savename_base = savenames[0] + "_" + nlo_outs[nlo_i]

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
# Extract observables (LAB)
# =========================
lo_th3, nlo_th3, full_th3 = lo["th3"], nlo["th3"], full["th3"]
lo_Emu, nlo_Emu, full_Emu = lo["Emu"], nlo["Emu"], full["Emu"]
lo_th5, nlo_th5, full_th5 = lo["th5"], nlo["th5"], full["th5"]
lo_Eph, nlo_Eph, full_Eph = lo["Eph"], nlo["Eph"], full["Eph"]
lo_th5_x, nlo_th5_x, full_th5_x = lo["th5_x"], nlo["th5_x"], full["th5_x"]
lo_th5_y, nlo_th5_y, full_th5_y = lo["th5_y"], nlo["th5_y"], full["th5_y"]

# =========================
# Extract observables (CMS)
# =========================
lo_th3_cms, nlo_th3_cms, full_th3_cms = lo["th3_cm"], nlo["th3_cm"], full["th3_cm"]
lo_Emu_cms, nlo_Emu_cms, full_Emu_cms = lo["Emu_cm"], nlo["Emu_cm"], full["Emu_cm"]
lo_th5_cms, nlo_th5_cms, full_th5_cms = lo["th5_cm"], nlo["th5_cm"], full["th5_cm"]
lo_Eph_cms, nlo_Eph_cms, full_Eph_cms = lo["Eph_cm"], nlo["Eph_cm"], full["Eph_cm"]
lo_th5_x_cms, nlo_th5_x_cms, full_th5_x_cms = lo["th5_x_"], nlo["th5_x_"], full["th5_x_"]
lo_th5_y_cms, nlo_th5_y_cms, full_th5_y_cms = lo["th5_y_"], nlo["th5_y_"], full["th5_y_"]

# =========================
# Colors
# =========================
colors = dict(
    lo="#1f77b4",
    nlo="#ff7f0e",
    full="#2ca02c",
    K="#d62728",
)

def make_plots_and_kfactors(
    *,
    tag,
    savename_base,
    lo_th3, nlo_th3, full_th3,
    lo_Emu, nlo_Emu, full_Emu,
    lo_th5, nlo_th5, full_th5,
    lo_Eph, nlo_Eph, full_Eph,
    lo_th5_x, nlo_th5_x, full_th5_x,
    lo_th5_y, nlo_th5_y, full_th5_y,
    outdir, outdir_vals,
    colors,
):
    savename = f"{savename_base}_{tag}"

    # ---------- write value files ----------
    write_file_with_values(outdir_vals + f"lo_th3_{savename}.txt",   lo_th3,   f"th3_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"nlo_th3_{savename}.txt",  nlo_th3,  f"th3_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"full_th3_{savename}.txt", full_th3, f"th3_{tag} bin centers",  "value")

    write_file_with_values(outdir_vals + f"lo_Emu_{savename}.txt",   lo_Emu,   f"Emu_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"nlo_Emu_{savename}.txt",  nlo_Emu,  f"Emu_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"full_Emu_{savename}.txt", full_Emu, f"Emu_{tag} bin centers",  "value")

    write_file_with_values(outdir_vals + f"lo_th5_{savename}.txt",   lo_th5,   f"th5_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"nlo_th5_{savename}.txt",  nlo_th5,  f"th5_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"full_th5_{savename}.txt", full_th5, f"th5_{tag} bin centers",  "value")

    write_file_with_values(outdir_vals + f"lo_Eph_{savename}.txt",   lo_Eph,   f"Eph_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"nlo_Eph_{savename}.txt",  nlo_Eph,  f"Eph_{tag} bin centers",  "value")
    write_file_with_values(outdir_vals + f"full_Eph_{savename}.txt", full_Eph, f"Eph_{tag} bin centers",  "value")

    write_file_with_values(outdir_vals + f"lo_th5_x_{savename}.txt",   lo_th5_x,   f"th5_x_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"nlo_th5_x_{savename}.txt",  nlo_th5_x,  f"th5_x_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"full_th5_x_{savename}.txt", full_th5_x, f"th5_x_{tag} bin centers", "value")

    write_file_with_values(outdir_vals + f"lo_th5_y_{savename}.txt",   lo_th5_y,   f"th5_y_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"nlo_th5_y_{savename}.txt",  nlo_th5_y,  f"th5_y_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"full_th5_y_{savename}.txt", full_th5_y, f"th5_y_{tag} bin centers", "value")


    fig, axes = create_figure(nrows=3, ncols=4, figsize=(24, 12), font_size=12)

    ax_th3,   ax_K_th3,   ax_Emu,   ax_K_Emu   = axes[0]
    ax_th5,   ax_K_th5,   ax_Eph,   ax_K_Eph   = axes[1]
    ax_th5_x, ax_K_th5_x, ax_th5_y, ax_K_th5_y = axes[2]

    # ---------- observable plots ----------
    lo_th3_s   = scaleplot(lo_th3,   1.e-3)  # -> mrad
    nlo_th3_s  = scaleplot(nlo_th3,  1.e-3)
    full_th3_s = scaleplot(full_th3, 1.e-3)
    plot_lo_nlo_full(ax_th3, lo_th3_s, nlo_th3_s, full_th3_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_th3, r"$\theta_3\ (\mathrm{mrad})$", "?", f"Muon Scattering Angle ({tag})")
    # keep your xlim choice only if it makes sense in cms; otherwise comment it out
    ax_th3.set_xlim(1.3, 1.7)

    lo_Emu_s   = scaleplot(lo_Emu,   1.e3)  # -> GeV
    nlo_Emu_s  = scaleplot(nlo_Emu,  1.e3)
    full_Emu_s = scaleplot(full_Emu, 1.e3)
    plot_lo_nlo_full(ax_Emu, lo_Emu_s, nlo_Emu_s, full_Emu_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_Emu, r"$E_\mu\ (\mathrm{GeV})$", "?", f"Energy of the Scattered Muon ({tag})")

    lo_th5_s   = scaleplot(lo_th5,   1.e-3)  # -> mrad
    nlo_th5_s  = scaleplot(nlo_th5,  1.e-3)
    full_th5_s = scaleplot(full_th5, 1.e-3)
    plot_lo_nlo_full(ax_th5, lo_th5_s, nlo_th5_s, full_th5_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_th5, r"$\theta_5\ (\mathrm{mrad})$", "?", f"Photon Scattering Angle ({tag})")
    ax_th5.set_xlim(-2., 2.)

    lo_Eph_s   = scaleplot(lo_Eph,   1.e3)  # -> GeV
    nlo_Eph_s  = scaleplot(nlo_Eph,  1.e3)
    full_Eph_s = scaleplot(full_Eph, 1.e3)
    plot_lo_nlo_full(ax_Eph, lo_Eph_s, nlo_Eph_s, full_Eph_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_Eph, r"$E_\gamma\ (\mathrm{GeV})$", "?", f"Photon Energy ({tag})")

    ax_th3.set_yscale("linear")

    lo_th5_x_s   = scaleplot(lo_th5_x,   1.e-3)
    nlo_th5_x_s  = scaleplot(nlo_th5_x,  1.e-3)
    full_th5_x_s = scaleplot(full_th5_x, 1.e-3)
    plot_lo_nlo_full(ax_th5_x, lo_th5_x_s, nlo_th5_x_s, full_th5_x_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_th5_x, r"$\theta_{5,x}\ (\mathrm{mrad})$", "?", f"Photon X-deflection ({tag})")
    ax_th5_x.set_xlim(-12., 12.)

    lo_th5_y_s   = scaleplot(lo_th5_y,   1.e-3)
    nlo_th5_y_s  = scaleplot(nlo_th5_y,  1.e-3)
    full_th5_y_s = scaleplot(full_th5_y, 1.e-3)
    plot_lo_nlo_full(ax_th5_y, lo_th5_y_s, nlo_th5_y_s, full_th5_y_s, colors,
                     labels=dict(lo="LO", nlo="NLO", full="LO + NLO"))
    style_sci_x(ax_th5_y, r"$\theta_{5,y}\ (\mathrm{mrad})$", "?", f"Photon Y-deflection ({tag})")
    ax_th5_y.set_xlim(-12., 12.)

    # ---------- K-factors ----------
    K_th3 = mergebins(divideplots(nlo_th3_s,  full_th3_s),  5)
    K_Emu = mergebins(divideplots(nlo_Emu_s,  full_Emu_s),  5)
    K_th5 = mergebins(divideplots(nlo_th5_s,  full_th5_s),  5)
    K_Eph = mergebins(divideplots(nlo_Eph_s,  full_Eph_s),  5)
    K_th5_x = mergebins(divideplots(nlo_th5_x_s, full_th5_x_s), 5)
    K_th5_y = mergebins(divideplots(nlo_th5_y_s, full_th5_y_s), 5)

    write_file_with_values(outdir_vals + f"K_theta_3_{savename}.txt", K_th3, f"th3_{tag} bin center", "K_th3")
    write_file_with_values(outdir_vals + f"K_E_mu_{savename}.txt",    K_Emu, f"Emu_{tag} bin center", "K_Emu")
    write_file_with_values(outdir_vals + f"K_theta_5_{savename}.txt", K_th5, f"th5_{tag} bin center", "K_th5")
    write_file_with_values(outdir_vals + f"K_E_gamma_{savename}.txt", K_Eph, f"Eph_{tag} bin center", "K_Eph")
    write_file_with_values(outdir_vals + f"K_theta_5x_{savename}.txt", K_th5_x, f"th5_x_{tag} bin center", "K_th5_x")
    write_file_with_values(outdir_vals + f"K_theta_5y_{savename}.txt", K_th5_y, f"th5_y_{tag} bin center", "K_th5_y")
    print(f"[INFO] K-factor files written for {tag}.")

    plot_K(ax_K_th3, K_th3, colors["K"], r"$\theta_3\ (\mathrm{mrad})$", r"$K(\theta_3)$")
    plot_K(ax_K_Emu, K_Emu, colors["K"], r"$E_\mu\ (\mathrm{GeV})$",     r"$K(E_\mu)$")
    plot_K(ax_K_th5, K_th5, colors["K"], r"$\theta_5\ (\mathrm{mrad})$", r"$K(\theta_5)$")
    plot_K(ax_K_Eph, K_Eph, colors["K"], r"$E_\gamma\ (\mathrm{GeV})$",  r"$K(E_\gamma)$")
    plot_K(ax_K_th5_x, K_th5_x, colors["K"], r"$\theta_{5,x}\ (\mathrm{mrad})$", r"$K(\theta_{5,x})$")
    plot_K(ax_K_th5_y, K_th5_y, colors["K"], r"$\theta_{5,y}\ (\mathrm{mrad})$", r"$K(\theta_{5,y})$")

    ax_K_th3.set_xlim(ax_th3.get_xlim())
    ax_K_Emu.set_xlim(ax_Emu.get_xlim())
    ax_K_th5.set_xlim(ax_th5.get_xlim())
    ax_K_Eph.set_xlim(ax_Eph.get_xlim())
    ax_K_th5_x.set_xlim(ax_th5_x.get_xlim())
    ax_K_th5_y.set_xlim(ax_th5_y.get_xlim())

    ax_K_th3.set_yscale("linear")
    ax_K_Emu.set_yscale("linear")
    ax_K_th5.set_yscale("linear")
    ax_K_th5_x.set_yscale("linear")
    ax_K_th5_y.set_yscale("linear")

    fig.subplots_adjust(hspace=0.40, wspace=0.30)
    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)



make_plots_and_kfactors(
    tag="lab",
    savename_base=savename_base,
    lo_th3=lo_th3, nlo_th3=nlo_th3, full_th3=full_th3,
    lo_Emu=lo_Emu, nlo_Emu=nlo_Emu, full_Emu=full_Emu,
    lo_th5=lo_th5, nlo_th5=nlo_th5, full_th5=full_th5,
    lo_Eph=lo_Eph, nlo_Eph=nlo_Eph, full_Eph=full_Eph,
    lo_th5_x=lo_th5_x, nlo_th5_x=nlo_th5_x, full_th5_x=full_th5_x,
    lo_th5_y=lo_th5_y, nlo_th5_y=nlo_th5_y, full_th5_y=full_th5_y,
    outdir=outdir, outdir_vals=outdir_vals,
    colors=colors,
)

make_plots_and_kfactors(
    tag="cms",
    savename_base=savename_base,
    lo_th3=lo_th3_cms, nlo_th3=nlo_th3_cms, full_th3=full_th3_cms,
    lo_Emu=lo_Emu_cms, nlo_Emu=nlo_Emu_cms, full_Emu=full_Emu_cms,
    lo_th5=lo_th5_cms, nlo_th5=nlo_th5_cms, full_th5=full_th5_cms,
    lo_Eph=lo_Eph_cms, nlo_Eph=nlo_Eph_cms, full_Eph=full_Eph_cms,
    lo_th5_x=lo_th5_x_cms, nlo_th5_x=nlo_th5_x_cms, full_th5_x=full_th5_x_cms,
    lo_th5_y=lo_th5_y_cms, nlo_th5_y=nlo_th5_y_cms, full_th5_y=full_th5_y_cms,
    outdir=outdir, outdir_vals=outdir_vals,
    colors=colors,
)