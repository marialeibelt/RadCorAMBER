from pymule import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt
import matplotlib as mpl

from plotting import *

# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to where you are (Office/Home)
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
#mpl.rcParams["text.usetex"] = False

# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to where you are (Office/Home)
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
homedir = "/home/marialei/AMBER_RadCor/"
outdir = homedir + "Figures/"
#homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"
#outdir = homedir + "Figures/"


# =========================
# Input definitions
# =========================
lo_outs = [
    "mp2mp_NLO_22_12",
    "mp2mp_NLO_12_01",
    "mp2mp_NLO_13_01",
    "mp2mp_NLO_19_01"
]

nlo_outs = lo_outs
savenames = ["combined", "combined_16Jan","combined_19Jan"]

# =========================
# Bin-Einstellungen aus Fortran Userfile
# =========================
nrbins = 500
min_val = np.array([0.3e-3, 95e3, -12e-3, 1e3])  # th3, Emu, th5, Eph
max_val = np.array([2.0e-3, 101e3, 12e-3, 100e3])

# Binbreiten berechnen
binwidths_th3  = np.full(nrbins, (max_val[0] - min_val[0]) / nrbins) * 1e3    # mrad
binwidths_Emu  = np.full(nrbins, (max_val[1] - min_val[1]) / nrbins) / 1e3    # GeV
binwidths_th5  = np.full(nrbins, (max_val[2] - min_val[2]) / nrbins) * 1e3    # mrad
binwidths_Eph  = np.full(nrbins, (max_val[3] - min_val[3]) / nrbins) / 1e3    # GeV


# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to what you want to analyze
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
lo_i = 3
nlo_i = 3
savename = savenames[2] + "_" + nlo_outs[nlo_i]


# =========================
# Physics setup
# =========================
# LO setup
setup(folder=homedir + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv
#print("conv: ",conv)

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
th3_lo_mrad   = lo["th3"] * 1e3
th3_nlo_mrad  = nlo["th3"] * 1e3
th3_full_mrad = full["th3"] * 1e3

plot_lo_nlo_full(
    ax_th3,
    th3_lo_mrad,
    th3_nlo_mrad,
    th3_full_mrad,
    binwidths_th3,
    colors,
    labels=dict(lo="LO", nlo="NLO", full="LO + NLO"),
)

style_sci_x(
    ax_th3,
    xlabel=r"$\theta_3\ (\mathrm{mrad})$",
    ylabel=r"?",
    title="Scattering Angle",
)

ax_th3.set_xlim(1.3, 1.7)



# =========================
# K-Factor
# =========================
# K = NLO / (LO + NLO), rebinned
thK = mergebins(divideplots(th3_nlo_mrad, th3_full_mrad), 5)

# Dummy binwidths (K is dimensionless)
binwidths_K = np.ones(len(thK))

plot_lo_nlo_full(
    ax_K,
    lo=thK,        # we pass the same array three times
    nlo=thK,       # to reuse the unified plotting function
    full=thK,
    colors=dict(lo=colors["K"], nlo=colors["K"], full=colors["K"]),
    labels=dict(full="K-factor"),
    binwidths=binwidths_K,
)

style_sci_x(
    ax_K,
    xlabel=r"$\theta_3\ (\mathrm{mrad})$",
    ylabel=r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$",
    title="K Factor",
    yscale="linear",
)

ax_K.set_xlim(ax_th3.get_xlim())



# =========================
# Emu
# =========================
lo_Emu_GeV   = lo["Emu"] / 1e3
nlo_Emu_GeV  = nlo["Emu"] / 1e3
full_Emu_GeV = full["Emu"] / 1e3

plot_lo_nlo_full(
    ax_Emu,
    lo_Emu_GeV,
    nlo_Emu_GeV,
    full_Emu_GeV,
    binwidths_Emu,
    colors,
    labels=dict(lo="LO", nlo="NLO", full="LO + NLO")
)

style_sci_x(
    ax_Emu,
    xlabel=r"$E_\mu\ (\mathrm{GeV})$",
    ylabel=r"?",
    title="Energy of the Scattered Muon",
)
ax_Emu.set_yscale("linear")
ax_Emu.set_ylim(-1e-8,1e-8)

# =========================
# th5
# =========================
th5_lo_mrad   = lo["th5"] * 1e3
th5_nlo_mrad  = nlo["th5"] * 1e3
th5_full_mrad = full["th5"] * 1e3

plot_lo_nlo_full(
    ax_th5,
    th5_lo_mrad,
    th5_nlo_mrad,
    th5_full_mrad,
    binwidths_th5,
    colors,
    labels=dict(lo="LO", nlo="NLO", full="LO + NLO")
)

style_sci_x(
    ax_th5,
    xlabel=r"$\theta_5\ (\mathrm{mrad})$",
    ylabel=r"?",
    title="Scattering Angle of the Photon",
)
ax_th5.set_yscale("linear")

# =========================
# Eph
# =========================
lo_Eph_GeV   = lo["Eph"] / 1e3
nlo_Eph_GeV  = nlo["Eph"]/ 1e3
full_Eph_GeV = full["Eph"]/ 1e3

plot_lo_nlo_full(
    ax_Eph,
    lo_Eph_GeV,
    nlo_Eph_GeV,
    full_Eph_GeV,
    binwidths_Eph,
    colors,
    labels=dict(lo="LO", nlo="NLO", full="LO + NLO")
)

style_sci_x(
    ax_Eph,
    xlabel=r"$E_\gamma\ (\mathrm{GeV})$",
    ylabel=r"?",
    title="Energy of the Photon",
)
ax_Eph.set_xlim(0.,35.)
ax_Eph.set_yscale("linear")

# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
# CHANGE according to where you are (Office/Home)
# !/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!/!
#mulify(fig, delx=0.0, dely=0.0)


# =========================
# Layout & Save
# =========================
fig.subplots_adjust(hspace=0.35, wspace=0.3)

save_figure(
    fig,
    savename,
    outdir=outdir,
)
