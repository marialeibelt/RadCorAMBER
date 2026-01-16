from pymule import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib as mpl

mpl.rcParams["text.usetex"] = False

from plotting import *


# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to where you are (Office/Home)
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# homedir = "/home/marialei/"
# outdir = homedir + "AMBER_RadCor/Figures/"
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"
outdir = homedir + "Figures/"


# =========================
# Input definitions
# =========================
lo_outs = [
    "mp2mp_NLO_22_12",
    "mp2mp_NLO_12_01",
    "mp2mp_NLO_13_01",
]

nlo_outs = lo_outs
savenames = ["combined", "combined_16Jan"]


# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to what you want to analyze
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
lo_i = 2
nlo_i = 2
savename = savenames[1] + "_" + nlo_outs[nlo_i]


# =========================
# Physics setup
# =========================
# LO setup
setup(folder=homedir + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

# NLO setup
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

# full observable
full = lo + nlo


# =========================
# Figure setup (2x3 grid)
# =========================
fig, axes = create_figure(
    nrows=2,
    ncols=3,
    figsize=(16, 9),
    font_size=12,
)

ax_th3 = axes[0, 0]
ax_Emu = axes[0, 1]
ax_K   = axes[0, 2]
ax_th5 = axes[1, 0]
ax_Eph = axes[1, 1]

# unused panels
axes[1, 2].axis("off")


# =========================
# Colors for consistent plots
# =========================
colors = dict(
    lo="#1f77b4",     # blue
    nlo="#ff7f0e",    # orange
    full="#2ca02c",   # green
    K="#d62728",      # red
    th5="#ff00c8",    # pink
    Eph="#ff00c8"     # pink
)


# =========================
# th3 
# =========================
th3_lo = lo["th3"]
th3_nlo = nlo["th3"]
th3_full = addplots(th3_lo, th3_nlo)

plot_lo_nlo_full(
    ax_th3,
    th3_lo,
    th3_nlo,
    th3_full,
    colors,
    labels=dict(lo="LO", nlo="NLO", full="LO + NLO"),
)

style_sci_x(
    ax_th3,
    xlabel=r"$\theta_3\ (\mathrm{rad})$",
    ylabel="Counts",
    title="Scattering Angle",
    sharex=True,
)

ax_th3.set_xlim(1.3e-3, 1.7e-3)


# =========================
# K-Factor
# =========================
thK = mergebins(divideplots(th3_nlo, th3_full), 5)

plot_errorband(ax_K, thK, colors["K"])
ax_K.plot([], [], color=colors["K"], label="K-factor")
ax_K.legend()

style_sci_x(
    ax_K,
    xlabel=r"$\theta_3\ (\mathrm{rad})$",
    ylabel=r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$",
    title="K Factor",
    yscale="linear",
)

ax_K.set_xlim(ax_th3.get_xlim())


# =========================
# Emu
# =========================
plot_lo_nlo_full(
    ax_Emu,
    lo["Emu"],
    nlo["Emu"],
    full["Emu"],
    colors,
    labels=dict(lo="LO", nlo="NLO correction", full="LO + NLO"),
)

style_sci_x(
    ax_Emu,
    xlabel=r"$E_\mu\ (\mathrm{MeV})$",
    ylabel="Counts",
    title="Energy of the Scattered Muon",
)

# =========================
# th5
# =========================
plot_lo_nlo_full(
    ax_th5,
    lo["th5"],
    nlo["th5"],
    full["th5"],
    colors,
    labels=dict(lo="LO", nlo="NLO correction", full="LO + NLO"),
)

style_sci_x(
    ax_th5,
    xlabel=r"$\theta_5\ (\mathrm{rad})$",
    ylabel="Counts",
    title="Scattering Angle of the Photon",
)

# =========================
# Eph
# =========================
plot_lo_nlo_full(
    ax_Eph,
    lo["Eph"],
    nlo["Eph"],
    full["Eph"],
    colors,
    labels=dict(lo="LO", nlo="NLO correction", full="LO + NLO"),
)

style_sci_x(
    ax_Eph,
    xlabel=r"$E_\gamma\ (\mathrm{MeV})$",
    ylabel="Counts",
    title="Energy of the Photon",
)


# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to where you are (Office/Home)
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# mulify(fig, delx=0.0, dely=0.0)


# =========================
# Layout & Save
# =========================
fig.subplots_adjust(hspace=0.35, wspace=0.3)

save_figure(
    fig,
    savename,
    outdir=outdir,
)
