from pymule import *
from theo_calc import *
from plotting import *
import sys

# =========================
# Paths
# =========================
homedir = "/home/marialei/AMBER_RadCor/"  # Laptop
#homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"  # Office
outdir_vals = homedir + "Vals/"

# =========================
# Dataset definitions
# =========================
lo_outs = ["mp2mp_NLO_12_05_BIG",
           "mp2mp_NLO_12_05_SMALL",
           "mp2mp_NLO_12_05_TH100MeV_BIG",
           "mp2mp_NLO_13_05_TH100MeV_SMALL",
           "mp2mp_NLO_13_05_TH500MeV_BIG",
           "mp2mp_NLO_13_05_TH500MeV_SMALL"]

labels = [
    ("200 MeV", "BIG",   0),
    ("200 MeV", "SMALL", 1),
    ("100 MeV", "BIG",   2),
    ("100 MeV", "SMALL", 3),
    ("500 MeV", "BIG",   4),
    ("500 MeV", "SMALL", 5),
]

outfile = outdir_vals + "crosssections_summary.txt"

# =========================
# Loop & collect
# =========================
results = []

for th_label, size_label, idx in labels:
    setup(folder=homedir + lo_outs[idx] + "/out")
    lo    = mergefks(sigma("mp2mp0")) * alpha**2 * conv
    nlo   = mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv
    full  = lo + nlo
    onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

    sigma_lo_mb   = lo.value    / 1000
    sigma_nlo_mb  = nlo.value   / 1000
    sigma_full_mb = full.value  / 1000
    sigma_Rph_mb  = onlyR.value / 1000

    Rate_lo   = calculate_rate(sigma_lo_mb)
    Rate_nlo  = calculate_rate(sigma_nlo_mb)
    Rate_full = calculate_rate(sigma_full_mb)
    Rate_Rph  = calculate_rate(sigma_Rph_mb)

    results.append((th_label, size_label,
                    sigma_lo_mb, sigma_nlo_mb, sigma_full_mb, sigma_Rph_mb,
                    Rate_lo, Rate_nlo, Rate_full, Rate_Rph))

# =========================
# Write summary file
# =========================
with open(outfile, "w") as f:
    header = (f"{'Threshold':<12} {'Range':<8} "
              f"{'sigma_LO (mb)':<20} {'sigma_NLO (mb)':<20} {'sigma_full (mb)':<20} {'sigma_Rph (mb)':<20} "
              f"{'Rate_LO (1/s)':<20} {'Rate_NLO (1/s)':<20} {'Rate_full (1/s)':<20} {'Rate_Rph (1/s)'}\n")
    f.write(header)
    f.write("-" * len(header) + "\n")
    for th, sz, slo, snlo, sfull, sRph, Rlo, Rnlo, Rfull, RRph in results:
        f.write(f"{th:<12} {sz:<8} "
                f"{str(slo):<20} {str(snlo):<20} {str(sfull):<20} {str(sRph):<20} "
                f"{str(Rlo):<20} {str(Rnlo):<20} {str(Rfull):<20} {str(RRph)}\n")

print(f"Saved to {outfile}")