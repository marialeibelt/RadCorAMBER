from pymule import *
import matplotlib.pyplot as plt
import numpy as np
from plotting import *
import sys

# =========================
# Paths
# =========================
#homedir = "/home/marialei/AMBER_RadCor/"  # Laptop
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"  # Office
outdir  = homedir + "Figures/"

# =========================
# Runs to compare (label, folder)
# =========================
runs = [
    ("thmu small", "23_06_200MeV_Q2big_xi01_thmusmall"),
    ("thmu norm",  "23_06_200MeV_Q2big_xi01"),
    ("thmu big",   "23_06_200MeV_Q2big_xi01_thmubig"),
]

savename = "compare_th5_thmu"

# =========================
# Colors
# =========================
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

# =========================
# Plot
# =========================
plt.rcParams.update({
    "font.size": 12,
    "font.family": "serif",
    "text.usetex": True,
})

fig, ax = plt.subplots(figsize=(7, 5))

for (label, folder), color in zip(runs, colors):
    setup(folder=homedir + folder + "/out")
    nlo   = mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv
    onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

    th5 = onlyR["th5"]
    th5_finite = finite_bins(th5)
    th5_finite = th5_finite[th5_finite[:, 1] != 0]

    x = th5_finite[:, 0] * 1e-3   # mrad -> rad, x-axis in mrad via scale below
    y = th5_finite[:, 1] * 1e-3   # µbarn/mrad
    e = th5_finite[:, 2] * 1e-3

    x_mrad = th5_finite[:, 0]     # keep in mrad for x-axis

    ax.plot(x_mrad, y, color=color, label=label)
    ax.fill_between(x_mrad, y - e, y + e, color=color, alpha=0.2)

ax.set_xlabel(r"$\theta_5\ (\mathrm{mrad})$")
ax.set_ylabel(r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$")
ax.set_title(r"Photon angle $\theta_5$ for different $\theta_\mu$ cuts")
ax.set_yscale("log")
ax.set_xlim(left=0.)
ax.legend(framealpha=0)
ax.grid(True, alpha=0.4, linestyle="dotted")

fig.tight_layout()

from pathlib import Path
Path(outdir).mkdir(parents=True, exist_ok=True)
fig.savefig(outdir + savename + ".png", dpi=300, bbox_inches="tight")
plt.close(fig)

print(f"Saved: {outdir}{savename}.png")
