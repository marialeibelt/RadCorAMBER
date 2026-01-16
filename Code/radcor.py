# =========================
# Imports
# =========================
from pymule import *
from matplotlib.ticker import ScalarFormatter
import matplotlib.pyplot as plt

from plotting import (
    create_figure,
    style_axis,
    save_figure
)

# =========================
# Input definitions
# =========================
lo_outs = [
    "mcmule-release",
    "mp2mp_NLO_22_12",
    "mp2mp_NLO_12_01",
    "mp2mp_NLO_13_01"
]

nlo_outs = [
    "mp2mp_testNLO",
    "mp2mp_NLO_22_12",
    "mp2mp_NLO_12_01",
    "mp2mp_NLO_13_01"
]

savenames = ["combined"]

lo_i = 3
nlo_i = 3
savename = savenames[0] + "_" + nlo_outs[nlo_i]

# =========================
# Physics setup
# =========================
# LO setup
setup(folder="/home/marialei/" + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

# NLO setup
setup(folder="/home/marialei/" + nlo_outs[nlo_i] + "/out")
nlo = (
    mergefks(
        sigma("mp2mpR"),
        sigma("mp2mpF"),
        anyxi=sigma("mp2mpA")
    )
    * alpha**3
    * conv
)

# full observable
full = lo + nlo

# =========================
# Figure setup (2x2 grid)
# =========================
fig, axes = create_figure(
    nrows=2,
    ncols=2,
    figsize=(16, 9),
    font_size=12
)

ax_th3 = axes[0, 0]   # top left: th3 spectrum
ax_Emu = axes[0, 1]   # top right: Emu spectrum
ax_K   = axes[1, 0]   # bottom left: K-factor

# bottom right empty
axes[1, 1].axis("off")

# =========================
# Colors for consistent plots
# =========================
color_lo = "#1f77b4"     # blue
color_nlo = "#ff7f0e"    # orange
color_full = "#2ca02c"   # green
color_K = "#d62728"      # red

# =========================
# Observable 1: th3 (top left)
# =========================
th3_lo = lo["th3"]
th3_nlo = nlo["th3"]
th3_added = addplots(th3_lo, th3_nlo)

plt.sca(ax_th3)
eb_lo = errorband(th3_lo)
eb_nlo = errorband(th3_nlo)
eb_added = errorband(th3_added)

for line in eb_lo:
    line.set_color(color_lo)
for line in eb_nlo:
    line.set_color(color_nlo)
for line in eb_added:
    line.set_color(color_full)

# legend with proper colors
ax_th3.plot([], [], color=color_lo, label="LO")
ax_th3.plot([], [], color=color_nlo, label="NLO")
ax_th3.plot([], [], color=color_full, label="LO + NLO")
ax_th3.legend()

# style axes
style_axis(
    ax_th3,
    xlabel=r"$\theta_3\ (\mathrm{rad})$",
    ylabel="Counts",
    title=r"Scattering Angle",
    yscale="log"
)

ax_th3.set_xlim(1.3e-3, 1.7e-3)
ax_th3.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax_th3.ticklabel_format(style="sci", axis="x", scilimits=(-3, 3))
ax_th3.tick_params(labelbottom=False)  # shared x with K-panel

# =========================
# K-Factor for th3 (bottom left)
# =========================
thK = mergebins(divideplots(th3_nlo, th3_added), 5)

plt.sca(ax_K)
eb_K = errorband(thK)
for line in eb_K:
    line.set_color(color_K)

ax_K.plot([], [], color=color_K, label="K-factor")
ax_K.legend()

style_axis(
    ax_K,
    xlabel=r"$\theta_3\ (\mathrm{rad})$",
    ylabel=r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$",
    title=r"K Factor",
    yscale="linear",
    legend=False
)

ax_K.set_xlim(ax_th3.get_xlim())
ax_K.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax_K.ticklabel_format(style="sci", axis="x", scilimits=(-3, 3))

# =========================
# Observable 2: Emu (top right)
# =========================
Emu_lo = lo["Emu"]
Emu_nlo = nlo["Emu"]
Emu_full = full["Emu"]

plt.sca(ax_Emu)
eb_Emu_lo = errorband(Emu_lo)
eb_Emu_nlo = errorband(Emu_nlo)
eb_Emu_full = errorband(Emu_full)

for line in eb_Emu_lo:
    line.set_color(color_lo)
for line in eb_Emu_nlo:
    line.set_color(color_nlo)
for line in eb_Emu_full:
    line.set_color(color_full)

ax_Emu.plot([], [], color=color_lo, label="LO")
ax_Emu.plot([], [], color=color_nlo, label="NLO correction")
ax_Emu.plot([], [], color=color_full, label="LO + NLO")
ax_Emu.legend()

style_axis(
    ax_Emu,
    xlabel=r"$E_\mu\ (\mathrm{MeV})$",
    ylabel="Counts",
    title="Energy of the Scattered Muon",
    yscale="log"
)

ax_Emu.xaxis.set_major_formatter(ScalarFormatter(useMathText=True))
ax_Emu.ticklabel_format(style="sci", axis="x", scilimits=(-3, 3))

# =========================
# Layout & Save
# =========================
mulify(fig, delx=0.0, dely=0.0)

# increase vertical spacing to avoid overlap
fig.subplots_adjust(hspace=0.4, wspace=0.3)

save_figure(
    fig,
    savename,
    outdir="/home/marialei/AMBER_RadCor/Figures/"
)
