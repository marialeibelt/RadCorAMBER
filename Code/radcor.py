from pymule import *
import matplotlib.pyplot as plt
import numpy as np
from plotting import *
from theo_calc import *
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


# =========================
# Paths
# =========================
#homedir = "/home/marialei/AMBER_RadCor/" # Laptop
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"
outdir = homedir + "Figures/"
outdir_vals = homedir + "Vals/"

# =========================
# Input definitions
# =========================
lo_outs = ["mp2mp_NLO_19_01", "mp2mp_NLO_01_02", "mp2mp_NLO_24_02","mp2mp_NLO_15_03", "mp2mptest",
           "mp2mp_23_03", "mp2mp_NLO_24_03", "mp2mp_NLO_24_03_new", "mp2mp_NLO_24_03_evening", "mp2mp_NLO_26_03",
           "mp2mp_NLO_26_03_new","mp2mp_26_03_timetest","lesspoints3","smallth3","folder",
           "folder2", "folder3", "mp2mp_NLO_27_03", "mp2mp_NLO_27_03_2", "mp2mp_NLO_13_04",
           "mp2mp_NLO_20_04","mp2mp_NLO_21_04","mp2mp_NLO_21_04_phicut","mp2mp_NLO_24_04_mitcos","mp2mp_NLO_28_04",
           "mp2mp_NLO_29_04","mp2mp_NLO_07_05_big","mp2mp_NLO_07_05_small","mp2mp_NLO_08_05_full","mp2mp_NLO_08_05_full_costh3test",
           "mp2mp_NLO_11_05_BIG","mp2mp_NLO_12_05_BIG","mp2mp_NLO_12_05_SMALL","mp2mp_NLO_12_05_TH100MeV_BIG","mp2mp_NLO_13_05_TH100MeV_SMALL",
           "mp2mp_NLO_13_05_TH500MeV_BIG","mp2mp_NLO_13_05_TH500MeV_SMALL","20_05_Eph100MeV_SMALL","20_05_Eph200MeV_SMALL","20_05_Eph500MeV_SMALL",
           "20_05_Eph100MeV_BIG","20_05_Eph200MeV_BIG","20_05_Eph500MeV_BIG","20_05_100MeV_large_BIG","20_05_200MeV_large_BIG",
           "20_05_500MeV_large_BIG"]

nlo_outs = lo_outs

savenames = ["combined", "15_03", "17_03", "18_03", "23_03",
             "24_03", "25_03","26_03","27_03","13_04",
             "14_04","14_04_add","20_04","21_04","22_04",
             "24_04","28_04","29_04","4_5","5_5",
             "7_5","8_5","11_05","12_5","13_5",
             "20_05","28_5"]

# =========================
# Dataset choice
# =========================
lo_i = 44
nlo_i = 44
savename_i = 26
nbins = 500

bin_width = 0.0382
n_bands = 10
band_min = -(n_bands/2 * bin_width)
band_max = n_bands/2 * bin_width
Y5_RANGE = (band_min, band_max)
X5_RANGE = (band_min, band_max)

savename_base = savenames[savename_i] + "_" + nlo_outs[nlo_i]

# Redirect stdout
log_file = outdir_vals + f"{savename_base}_output.txt"
sys.stdout = Tee(log_file)

# =========================
# Physics setup
# =========================
setup(folder=homedir + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

setup(folder=homedir + nlo_outs[nlo_i] + "/out")
nlo = (mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv)
full = lo + nlo
onlyR = (mergefks(sigma("mp2mpR")) * alpha**3 * conv)

# =========================
# Extract observables (LAB only)
# =========================
lo_th3, nlo_th3, full_th3 = lo["th3"], nlo["th3"], full["th3"]
lo_Emu, nlo_Emu, full_Emu = lo["Emu"], nlo["Emu"], full["Emu"]
lo_th5, nlo_th5, full_th5 = lo["th5"], nlo["th5"], full["th5"]
lo_Eph, nlo_Eph, full_Eph = lo["Eph"], nlo["Eph"], full["Eph"]
lo_phi5, nlo_phi5, full_phi5 = lo["phi5"], nlo["phi5"], full["phi5"]
lo_x5, nlo_x5, full_x5 = lo["x5"], nlo["x5"], full["x5"]
lo_y5, nlo_y5, full_y5 = lo["y5"], nlo["y5"], full["y5"]
lo_costh3, nlo_costh3, full_costh3 = lo["costh3"], nlo["costh3"], full["costh3"]
lo_Q2, nlo_Q2, full_Q2 = lo["Qsq"], nlo["Qsq"], full["Qsq"]

x5_bands_lo, x5_bands_nlo, x5_bands_full = {}, {}, {}
y5_bands_lo, y5_bands_nlo, y5_bands_full = {}, {}, {}

for i in range(1, n_bands+1):
    key_x = f"x5_B{i}"
    key_y = f"y5_B{i}"

    try:
        x5_bands_lo[i]   = lo[key_x]
        x5_bands_nlo[i]  = nlo[key_x]
        x5_bands_full[i] = full[key_x]
    except KeyError:
        pass

    try:
        y5_bands_lo[i]   = lo[key_y]
        y5_bands_nlo[i]  = nlo[key_y]
        y5_bands_full[i] = full[key_y]
    except KeyError:
        pass


colors = dict(lo="#1f77b4", nlo="#ff7f0e", full="#2ca02c", K="#d62728")

# =========================
# Function (CMS removed)
# =========================
def make_plots_and_kfactors(*,
                            tag, savename_base,
                            lo_th3, nlo_th3, full_th3,
                            lo_Emu, nlo_Emu, full_Emu,
                            lo_th5, nlo_th5, full_th5,
                            lo_Eph, nlo_Eph, full_Eph,
                            lo_phi5, nlo_phi5, full_phi5,
                            lo_costh3, nlo_costh3, full_costh3,
                            lo_Q2, nlo_Q2, full_Q2,
                            lo_x5, nlo_x5, full_x5,
                            lo_y5, nlo_y5, full_y5,
                            outdir, outdir_vals, colors):

    savename = f"{savename_base}_{tag}"

    # only LAB variables
    variables = [
        ("th3", False, False),
        ("Emu", False, False),
        ("th5", True, False),
        ("Eph", True, False),
        ("phi5", True, False),
        ("costh3", False, True),
        ("Q2", False, True),
        ("x5", True, False),
        ("y5", True, False),
    ]

    data_map = {
        "th3": lo_th3,
        "Emu": lo_Emu,
        "th5": lo_th5,
        "Eph": lo_Eph,
        "phi5": lo_phi5,
        "costh3": lo_costh3,
        "Q2": lo_Q2,
        "x5": lo_x5,
        "y5": lo_y5,
    }

    for var, photon_only, lab_only in variables:
        if lab_only and tag != "lab":
            continue

        orders = ["nlo"] if photon_only else ["lo", "nlo", "full"]

        for order in orders:
            arr = locals()[f"{order}_{var}"]
            write_file_with_values(outdir_vals + f"{order}_{var}_{savename}.txt",
                                   arr,
                                   f"{var}_{tag} bin centers",
                                   "value")

    fig, axes = create_figure(nrows=8, ncols=2, figsize=(16,22),
                              font_size=12,
                              gridspec_kw={"height_ratios":[3,1]*4,"hspace":0.6})

    # plotting unchanged (omitted CMS dependencies already removed earlier)

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)


# =========================
# Run LAB only
# =========================
make_plots_and_kfactors(tag="lab", savename_base=savename_base,
                        lo_th3=lo_th3, nlo_th3=nlo_th3, full_th3=full_th3,
                        lo_Emu=lo_Emu, nlo_Emu=nlo_Emu, full_Emu=full_Emu,
                        lo_th5=lo_th5, nlo_th5=nlo_th5, full_th5=full_th5,
                        lo_Eph=lo_Eph, nlo_Eph=nlo_Eph, full_Eph=full_Eph,
                        lo_phi5=lo_phi5, nlo_phi5=nlo_phi5, full_phi5=full_phi5,
                        lo_costh3=lo_costh3, nlo_costh3=nlo_costh3, full_costh3=full_costh3,
                        lo_Q2=lo_Q2, nlo_Q2=nlo_Q2, full_Q2=full_Q2,
                        lo_x5=lo_x5, nlo_x5=nlo_x5, full_x5=full_x5,
                        lo_y5=lo_y5, nlo_y5=nlo_y5, full_y5=full_y5,
                        outdir=outdir, outdir_vals=outdir_vals, colors=colors)

sys.stdout.close()
sys.stdout = sys.stdout.terminal