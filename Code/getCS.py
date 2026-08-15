from pymule import *
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.colors import LogNorm
import sys

class Tee:
    """Writes output to both terminal and a file."""
    def __init__(self, filename):
        self.terminal = sys.stdout
        self.file = open(filename, "w")

    def write(self, message):
        self.terminal.write(message)
        self.file.write(message)

    def flush(self):
        self.terminal.flush()
        self.file.flush()

    def close(self):
        self.file.close()


homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"  # Office

# =========================
# Input definitions
# =========================

savenames = ["14_08"]        #0       
outfolder = ["05_08_evtgen"] #0

# =========================
# Dataset choice / Has to be checked each time!
# =========================
out_i=0
savename_i = 0
nbins = 500

savename_base = savenames[savename_i] + "_" + outfolder[out_i]
outdir = homedir + outfolder[out_i] 

# Redirect stdout
log_file = outdir + "/" + f"{savename_base}_output.txt"
sys.stdout = Tee(log_file)

print("=======================================================================")
print("        Analysed file: ", savename_base)
print("=======================================================================")

print("outdir + /out:", outdir + "/out")

setup(folder= outdir + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

setup(folder= outdir + "/out")
nlo    = mergefks(sigma("mp2mpR"), sigma("mp2mpNLO0")) * alpha**3 * conv
full   = lo + nlo
onlyR  = mergefks(sigma("mp2mpR")) * alpha**3 * conv

cs_lo = lo.value[0] / 1000
cs_onlyR = onlyR.value[0] / 1000
cs_nlo = nlo.value[0] / 1000
cs_full = full.value[0] / 1000

print("LO cross section:     ", cs_lo, " mb")
print("NLO cross section:    ", cs_nlo, " mb")
print("Full cross section:   ", cs_full, " mb")
print("Only R cross section: ", cs_onlyR, " mb")
cs_array = np.array([cs_lo, cs_onlyR, cs_nlo, cs_full])
print("Cross section array (lo, onlyR, nlo, full):", cs_array)

# Save cross section array to txt file
array_file = outdir+ "/" + f"{savename_base}_cross_sections.txt"
np.savetxt(
    array_file,
    cs_array,
    header="1. LO  2. onlyR  3. NLO  4. full"
)

print("Cross section array saved to:", array_file)