from pymule import *
import matplotlib.pyplot as plt
import numpy as np

from plotting import *

# =========================
# Paths
# =========================
#homedir = "/home/marialei/AMBER_RadCor/"                                           #Laptop
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"    #Office
outdir = homedir + "Figures/"
outdir_vals = homedir + "Vals/"

# =========================
# Input definitions
# =========================
lo_outs = [
    "mp2mp_NLO_19_01",
    "mp2mp_NLO_01_02",
    "mp2mp_NLO_24_02",
    "mp2mp_NLO_15_03"
]
nlo_outs = lo_outs
savenames = ["combined", "15_03", "17_03"]

# =========================
# Dataset choice
# =========================
lo_i = 3
nlo_i = 3
savename_base = savenames[2] + "_" + nlo_outs[nlo_i]

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
lo_phi5_x, nlo_phi5_x, full_phi5_x = lo["phi5_x"], nlo["phi5_x"], full["phi5_x"]
lo_phi5_y, nlo_phi5_y, full_phi5_y = lo["phi5_y"], nlo["phi5_y"], full["phi5_y"]

# =========================
# Extract observables (CMS)
# =========================
lo_th3_cms, nlo_th3_cms, full_th3_cms = lo["th3_cms"], nlo["th3_cms"], full["th3_cms"]
lo_Emu_cms, nlo_Emu_cms, full_Emu_cms = lo["Emu_cms"], nlo["Emu_cms"], full["Emu_cms"]
lo_th5_cms, nlo_th5_cms, full_th5_cms = lo["th5_cms"], nlo["th5_cms"], full["th5_cms"]
lo_Eph_cms, nlo_Eph_cms, full_Eph_cms = lo["Eph_cms"], nlo["Eph_cms"], full["Eph_cms"]
lo_phi5_x_cms, nlo_phi5_x_cms, full_phi5_x_cms = lo["phi5_x_cms"], nlo["phi5_x_cms"], full["phi5_x_cms"]
lo_phi5_y_cms, nlo_phi5_y_cms, full_phi5_y_cms = lo["phi5_y_cms"], nlo["phi5_y_cms"], full["phi5_y_cms"]

# =========================
# Colors
# =========================
colors = dict(
    lo="#1f77b4",
    nlo="#ff7f0e",
    full="#2ca02c",
    K="#d62728",
)


def draw_observable_and_k(
    ax_main,
    ax_k,
    *,
    lo_hist,
    nlo_hist,
    full_hist,
    scale_factor,
    x_label_main,
    x_label_k,
    y_label_main,
    main_title,
    xlim=None,
    main_yscale="log",
    force_main_linear=False,
    colors=None,
    hide_main_xticks=True,
):
    lo_s   = scaleplot(lo_hist, scale_factor)
    nlo_s  = scaleplot(nlo_hist, scale_factor)
    full_s = scaleplot(full_hist, scale_factor)

    plot_lo_nlo_full(
        ax_main,
        lo_s,
        nlo_s,
        full_s,
        colors,
        labels=dict(lo="LO", nlo="NLO", full="LO + NLO"),
    )

    yscale_to_use = main_yscale
    if force_main_linear:
        yscale_to_use = "linear"
    else:
        if np.all(full_s[:, 1] <= 0):
            yscale_to_use = "linear"

    style_sci_x(
        ax_main,
        x_label_main,
        y_label_main,
        main_title,
        yscale=yscale_to_use,
        sharex=False,
    )

    if xlim is not None:
        ax_main.set_xlim(*xlim)

    K = mergebins(divideplots(nlo_s, full_s), 5)

    plot_K(
        ax_k,
        K,
        colors["K"],
        x_label_k,
    )

    if xlim is not None:
        ax_k.set_xlim(*xlim)
    else:
        ax_k.set_xlim(ax_main.get_xlim())

    if hide_main_xticks:
        ax_main.tick_params(axis="x", labelbottom=False)
        ax_main.set_xlabel(None)

    return lo_s, nlo_s, full_s, K


def save_single_pair_plot(
    *,
    savename,
    lo_hist,
    nlo_hist,
    full_hist,
    scale_factor,
    x_label,
    y_label,
    main_title,
    xlim=None,
    main_yscale="log",
    force_main_linear=False,
    colors=None,
    outdir=None,
):
    fig, axes = create_figure(
        nrows=2,
        ncols=1,
        figsize=(7, 6),
        font_size=12,
        sharex=True,
        gridspec_kw={
            "height_ratios": [3, 1],
            "hspace": 0,
        },
    )

    ax_main = axes[0, 0]
    ax_k    = axes[1, 0]

    draw_observable_and_k(
        ax_main,
        ax_k,
        lo_hist=lo_hist,
        nlo_hist=nlo_hist,
        full_hist=full_hist,
        scale_factor=scale_factor,
        x_label_main=None,
        x_label_k=x_label,
        y_label_main= y_label,
        main_title=main_title,
        xlim=xlim,
        main_yscale=main_yscale,
        force_main_linear=force_main_linear,
        colors=colors,
        hide_main_xticks=True,
    )

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)


def make_plots_and_kfactors(
    *,
    tag,
    savename_base,
    lo_th3, nlo_th3, full_th3,
    lo_Emu, nlo_Emu, full_Emu,
    lo_th5, nlo_th5, full_th5,
    lo_Eph, nlo_Eph, full_Eph,
    lo_phi5_x, nlo_phi5_x, full_phi5_x,
    lo_phi5_y, nlo_phi5_y, full_phi5_y,
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

    write_file_with_values(outdir_vals + f"lo_phi5_x_{savename}.txt",   lo_phi5_x,   f"phi5_x_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"nlo_phi5_x_{savename}.txt",  nlo_phi5_x,  f"phi5_x_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"full_phi5_x_{savename}.txt", full_phi5_x, f"phi5_x_{tag} bin centers", "value")

    write_file_with_values(outdir_vals + f"lo_phi5_y_{savename}.txt",   lo_phi5_y,   f"phi5_y_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"nlo_phi5_y_{savename}.txt",  nlo_phi5_y,  f"phi5_y_{tag} bin centers", "value")
    write_file_with_values(outdir_vals + f"full_phi5_y_{savename}.txt", full_phi5_y, f"phi5_y_{tag} bin centers", "value")

    # ---------- combined figure ----------
    fig, axes = create_figure(
        nrows=6,
        ncols=2,
        figsize=(16, 16),
        font_size=12,
        sharex='col',
        gridspec_kw={
            "height_ratios": [3, 1, 3, 1, 3, 1],
            "hspace": 0.5,
        },
    )

    ax_th3,   ax_Emu       = axes[0]
    ax_K_th3, ax_K_Emu     = axes[1]

    ax_th5,   ax_Eph       = axes[2]
    ax_K_th5, ax_K_Eph     = axes[3]

    ax_phi5_x, ax_phi5_y     = axes[4]
    ax_K_phi5_x, ax_K_phi5_y = axes[5]

    _, _, _, K_th3 = draw_observable_and_k(
        ax_th3, ax_K_th3,
        lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3,
        scale_factor=1.e-3,
        x_label_main=r"$\theta_3\ (\mathrm{mrad})$",
        x_label_k=r"$\theta_3\ (\mathrm{mrad})$",
        y_label_main=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Muon Scattering Angle ({tag})",
        xlim=(1.3, 1.7),
        main_yscale="log",
        force_main_linear=True,
        colors=colors,
        hide_main_xticks=False,
    )

    _, _, _, K_Emu = draw_observable_and_k(
        ax_Emu, ax_K_Emu,
        lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu,
        scale_factor=1.e3,
        x_label_main=r"$E_\mu\ (\mathrm{GeV})$",
        x_label_k=r"$E_\mu\ (\mathrm{GeV})$",
        y_label_main=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",  
        main_title=f"Energy of the Scattered Muon ({tag})",
        colors=colors,
        hide_main_xticks=False,
    )

    _, _, _, K_th5 = draw_observable_and_k(
        ax_th5, ax_K_th5,
        lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5,
        scale_factor=1.e-3,
        x_label_main=r"$\theta_5\ (\mathrm{mrad})$",
        x_label_k=r"$\theta_5\ (\mathrm{mrad})$",
        y_label_main=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Scattering Angle ({tag})",
        xlim=(-2., 2.),
        colors=colors,
        hide_main_xticks=False,
    )

    _, _, _, K_Eph = draw_observable_and_k(
        ax_Eph, ax_K_Eph,
        lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph,
        scale_factor=1.e3,
        x_label_main=r"$E_\gamma\ (\mathrm{GeV})$",
        x_label_k=r"$E_\gamma\ (\mathrm{GeV})$",
        y_label_main=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",  
        main_title=f"Photon Energy ({tag})",
        colors=colors,
        hide_main_xticks=False,
    )

    _, _, _, K_phi5_x = draw_observable_and_k(
        ax_phi5_x, ax_K_phi5_x,
        lo_hist=lo_phi5_x, nlo_hist=nlo_phi5_x, full_hist=full_phi5_x,
        scale_factor=1.e-3,
        x_label_main=r"$\phi_{5,x}\ (\mathrm{mrad})$",
        x_label_k=r"$\phi_{5,x}\ (\mathrm{mrad})$",
        y_label_main=r"$\frac{d\sigma}{d\phi_{5,x}}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon X-deflection ({tag})",
        xlim=(-12., 12.),
        colors=colors,
        hide_main_xticks=False,
    )

    _, _, _, K_phi5_y = draw_observable_and_k(
        ax_phi5_y, ax_K_phi5_y,
        lo_hist=lo_phi5_y, nlo_hist=nlo_phi5_y, full_hist=full_phi5_y,
        scale_factor=1.e-3,
        x_label_main=r"$\phi_{5,y}\ (\mathrm{mrad})$",
        x_label_k=r"$\phi_{5,y}\ (\mathrm{mrad})$",
        y_label_main=r"$\frac{d\sigma}{d\phi_{5,y}}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Y-deflection ({tag})",
        xlim=(-12., 12.),
        colors=colors,
        hide_main_xticks=False,
    )

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)

    # ---------- write K-value files ----------
    write_file_with_values(outdir_vals + f"K_theta_3_{savename}.txt",  K_th3,   f"th3_{tag} bin center",   "K_th3")
    write_file_with_values(outdir_vals + f"K_E_mu_{savename}.txt",     K_Emu,   f"Emu_{tag} bin center",   "K_Emu")
    write_file_with_values(outdir_vals + f"K_theta_5_{savename}.txt",  K_th5,   f"th5_{tag} bin center",   "K_th5")
    write_file_with_values(outdir_vals + f"K_E_gamma_{savename}.txt",  K_Eph,   f"Eph_{tag} bin center",   "K_Eph")
    write_file_with_values(outdir_vals + f"K_phi_5x_{savename}.txt", K_phi5_x, f"phi5_x_{tag} bin center", "K_phi5_x")
    write_file_with_values(outdir_vals + f"K_phi_5y_{savename}.txt", K_phi5_y, f"phi5_y_{tag} bin center", "K_phi5_y")

    print(f"[INFO] K-factor files written for {tag}.")

    # ---------- separate pair plots ----------
    save_single_pair_plot(
        savename=f"{savename}_th3_pair",
        lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3,
        scale_factor=1.e-3,
        x_label=r"$\theta_3\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Muon Scattering Angle ({tag})",
        xlim=(1.3, 1.7),
        force_main_linear=True,
        colors=colors,
        outdir=outdir,
    )

    save_single_pair_plot(
        savename=f"{savename}_Emu_pair",
        lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu,
        scale_factor=1.e3,
        x_label=r"$E_\mu\ (\mathrm{GeV})$",
        y_label=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",   
        main_title=f"Energy of the Scattered Muon ({tag})",
        colors=colors,
        outdir=outdir,
    )

    save_single_pair_plot(
        savename=f"{savename}_th5_pair",
        lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5,
        scale_factor=1.e-3,
        x_label=r"$\theta_5\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Scattering Angle ({tag})",
        xlim=(-2., 2.),
        colors=colors,
        outdir=outdir,
    )

    save_single_pair_plot(
        savename=f"{savename}_Eph_pair",
        lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph,
        scale_factor=1.e3,
        x_label=r"$E_\gamma\ (\mathrm{GeV})$",
        y_label=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",   
        main_title=f"Photon Energy ({tag})",
        colors=colors,
        outdir=outdir,
    )

    save_single_pair_plot(
        savename=f"{savename}_phi5_x_pair",
        lo_hist=lo_phi5_x, nlo_hist=nlo_phi5_x, full_hist=full_phi5_x,
        scale_factor=1.e-3,
        x_label=r"$\phi_{5,x}\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\phi_{5,x}}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon X-deflection ({tag})",
        xlim=(-12., 12.),
        colors=colors,
        outdir=outdir,
    )

    save_single_pair_plot(
        savename=f"{savename}_phi5_y_pair",
        lo_hist=lo_phi5_y, nlo_hist=nlo_phi5_y, full_hist=full_phi5_y,
        scale_factor=1.e-3,
        x_label=r"$\phi_{5,y}\ (\mathrm{mrad})$",
        y_label=r"$\frac{d\sigma}{d\phi_{5,y}}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        main_title=f"Photon Y-deflection ({tag})",
        xlim=(-12., 12.),
        colors=colors,
        outdir=outdir,
    )


make_plots_and_kfactors(
    tag="lab",
    savename_base=savename_base,
    lo_th3=lo_th3, nlo_th3=nlo_th3, full_th3=full_th3,
    lo_Emu=lo_Emu, nlo_Emu=nlo_Emu, full_Emu=full_Emu,
    lo_th5=lo_th5, nlo_th5=nlo_th5, full_th5=full_th5,
    lo_Eph=lo_Eph, nlo_Eph=nlo_Eph, full_Eph=full_Eph,
    lo_phi5_x=lo_phi5_x, nlo_phi5_x=nlo_phi5_x, full_phi5_x=full_phi5_x,
    lo_phi5_y=lo_phi5_y, nlo_phi5_y=nlo_phi5_y, full_phi5_y=full_phi5_y,
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
    lo_phi5_x=lo_phi5_x_cms, nlo_phi5_x=nlo_phi5_x_cms, full_phi5_x=full_phi5_x_cms,
    lo_phi5_y=lo_phi5_y_cms, nlo_phi5_y=nlo_phi5_y_cms, full_phi5_y=full_phi5_y_cms,
    outdir=outdir, outdir_vals=outdir_vals,
    colors=colors,
)