from pymule import *
import matplotlib.pyplot as plt
import numpy as np
from plotting import *
import sys
from pathlib import Path

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
print(runs)
print(len(runs))

savename = "compare_th5_thmu"

# =========================
# Colors
# =========================
colors = ["#1f77b4", "#ff7f0e", "#2ca02c"]

# =========================
# Plot
# =========================

fig, ax = plt.subplots(figsize=(7, 5))

print("Start loop")
for (label, folder), color in zip(runs, colors):
    setup(folder=homedir + folder + "/out")
    nlo   = mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv
    onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

    th5 = onlyR["th5"]
    th5_finite = finite_bins(th5)

    print("shape =", th5_finite.shape)
    if len(th5_finite) == 0:
        print("th5_finite is empty!")
        continue

    x = th5_finite[:, 0] * 1e3   # mrad
    y = th5_finite[:, 1] * 1e-3   # µbarn/mrad
    e = th5_finite[:, 2] * 1e-3

    x_mrad = th5_finite[:, 0] * 1e3

    print(np.min(th5_finite[:, 0]), np.max(th5_finite[:, 0]))

    ax.plot(x_mrad, y, color=color, label=label)
print("End loop")
ax.set_xlabel(r"$\theta_5\ (\mathrm{mrad})$")
ax.set_ylabel(r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$")
ax.set_title(r"Photon angle $\theta_5$ for different $\theta_\mu$ cuts")
ax.set_yscale("log")
ax.set_xlim(-1, 14.)
ax.legend(framealpha=0)
ax.grid(True, alpha=0.4, linestyle="dotted")

fig.tight_layout()
fig.savefig(outdir + savename + ".png", dpi=300, bbox_inches="tight")

plt.close(fig)
print(f"Saved: {outdir}{savename}.png")
