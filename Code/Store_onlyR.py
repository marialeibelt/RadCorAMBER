from pymule import *
import numpy as np
from plotting import *
from theo_calc import *

# =========================
# Paths
# =========================
#homedir = "/home/marialei/AMBER_RadCor/" # Laptop
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"
outdir_vals = homedir + "Vals/"

# =========================
# Input definitions
# =========================
outs = [
    "26_06_200MeV_Q2big_xi001_onlyR","26_06_200MeV_Q2big_xi0001_onlyR",
    "01_07_200MeV_Q2big_xi01_onlyR_thmubig","01_07_200MeV_Q2big_xi01_onlyR_thmusmall",
    "01_07_100MeV_Q2big_xi01_onlyR","25_06_200MeV_Q2big_xi01_onlyR","01_07_500MeV_Q2big_xi01_onlyR","08_07_200MeV_Q2big_xi01_onlyR_Mehran_precthmurange","09_07_200MeV_xi01_onlyR_Mehran",
    "14_07_02to70GeV_xi01_onlyR_Mehran" 
]

savename_base = "14_07"
table_file_cs_R = outdir_vals + f"{savename_base}_CS_R.txt"

# =========================
# Column headers
# =========================
headers = [
    "file",
    "CS [mb]",  "CS err",
    "R [1/s]",  "R err",
]

# =========================
# Fill rows
# =========================
first_col_width = 45
col_width = 25
rows = []
for out in outs:
    setup(folder=homedir + out + "/out")
    onlyR = mergefks(sigma("mp2mpR")) * alpha**3 * conv

    CS_Rph_val,  CS_Rph_err  = onlyR.value /1000

    R_Rph_val,  R_Rph_err  = calculate_rate(onlyR.value /1000)

    rows.append([
        out,
        CS_Rph_val,  CS_Rph_err,
        R_Rph_val,  R_Rph_err,
    ])

# =========================
# Write tables
# =========================
def write_table(filename, headers, rows, comment=""):
    with open(filename, "w") as f:
        if comment:
            f.write(f"# {comment}\n\n")
        f.write(headers[0].ljust(first_col_width))
        f.write("".join(h.ljust(col_width) for h in headers[1:]))
        f.write("\n")
        for row in rows:
            line = row[0].ljust(first_col_width)
            line += "".join(f"{v:.6e}".ljust(col_width) for v in row[1:])
            f.write(line + "\n")
    print(f"[INFO] Table written to {filename}")

write_table(table_file_cs_R, headers, rows, comment="0.001 < Q2 (GeV2/c2) < 0.04")