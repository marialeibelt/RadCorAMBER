from pymule import *
import numpy as np
from plotting import *
from theo_calc import *

# =========================
# Paths
# =========================
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"
outdir_vals = homedir + "Vals/"

# =========================
# Input definitions
# =========================
lo_outs = [
    "20_05_Eph100MeV_SMALL", "20_05_Eph200MeV_SMALL", "20_05_Eph500MeV_SMALL",
    "20_05_Eph100MeV_BIG",   "20_05_Eph200MeV_BIG",   "20_05_Eph500MeV_BIG",
    "20_05_100MeV_large_BIG","20_05_200MeV_large_BIG", "20_05_500MeV_large_BIG",
    "09_06_thmunorm_noEphcut_noxycut","09_06_thmu4xsmaller_noEphcut_noxycut","09_06_thm4xbigger_noEphcut_noxycut",
]
nlo_outs = lo_outs

savename_base = "09_06"
table_file_cs = outdir_vals + f"{savename_base}_CS.txt"
table_file_r  = outdir_vals + f"{savename_base}_Rate.txt"

# =========================
# Column headers
# =========================
headers_cs = [
    "file",
    "Rph CS [mb]",  "Rph CS err",
    "LO CS [mb]",   "LO CS err",
    "NLO CS [mb]",  "NLO CS err",
    "Full CS [mb]", "Full CS err",
]

headers_r = [
    "file",
    "Rph R [1/s]",  "Rph R err",
    "LO R [1/s]",   "LO R err",
    "NLO R [1/s]",  "NLO R err",
    "Full R [1/s]", "Full R err",
]

# =========================
# Fill rows
# =========================
col_width = 25
rows_cs = []
rows_r  = []

for lo_out, nlo_out in zip(lo_outs, nlo_outs):
    setup(folder=homedir + lo_out + "/out")
    lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

    setup(folder=homedir + nlo_out + "/out")
    nlo   = mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv
    full  = lo + nlo
    onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

    CS_Rph_val,  CS_Rph_err  = onlyR.value /1000
    CS_LO_val,   CS_LO_err   = lo.value /1000
    CS_NLO_val,  CS_NLO_err  = nlo.value /1000
    CS_Full_val, CS_Full_err = full.value /1000

    R_Rph_val,  R_Rph_err  = calculate_rate(onlyR.value /1000)
    R_LO_val,   R_LO_err   = calculate_rate(lo.value /1000)
    R_NLO_val,  R_NLO_err  = calculate_rate(nlo.value /1000)
    R_Full_val, R_Full_err = calculate_rate(full.value /1000)

    rows_cs.append([
        lo_out,
        CS_Rph_val,  CS_Rph_err,
        CS_LO_val,   CS_LO_err,
        CS_NLO_val,  CS_NLO_err,
        CS_Full_val, CS_Full_err,
    ])

    rows_r.append([
        lo_out,
        R_Rph_val,  R_Rph_err,
        R_LO_val,   R_LO_err,
        R_NLO_val,  R_NLO_err,
        R_Full_val, R_Full_err,
    ])

# =========================
# Write tables
# =========================
def write_table(filename, headers, rows, comment=""):
    with open(filename, "w") as f:
        if comment:
            f.write(f"# {comment}\n\n")
        f.write("".join(h.ljust(col_width) for h in headers) + "\n")
        f.write("-" * (col_width * len(headers)) + "\n")
        for row in rows:
            line = row[0].ljust(col_width)
            line += "".join(f"{v:.6e}".ljust(col_width) for v in row[1:])
            f.write(line + "\n")
    print(f"[INFO] Table written to {filename}")

write_table(table_file_cs, headers_cs, rows_cs, comment="0.001 < Q2 (GeV2/c2) < 0.04 — Cross Sections")
write_table(table_file_r,  headers_r,  rows_r,  comment="0.001 < Q2 (GeV2/c2) < 0.04 — Rates")