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
homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"  # Office
outdir = homedir + "Figures/"
outdir_vals = homedir + "Vals/"

# =========================
# Input definitions
# =========================
lo_outs = ["mp2mp_NLO_19_01", "mp2mp_NLO_01_02", "mp2mp_NLO_24_02","mp2mp_NLO_15_03", "mp2mptest",                                              #0-4
           "mp2mp_23_03", "mp2mp_NLO_24_03", "mp2mp_NLO_24_03_new", "mp2mp_NLO_24_03_evening", "mp2mp_NLO_26_03",                               #5-9
           "mp2mp_NLO_26_03_new","mp2mp_26_03_timetest","lesspoints3","smallth3","folder",                                                      #10-14
           "folder2", "folder3", "mp2mp_NLO_27_03", "mp2mp_NLO_27_03_2", "mp2mp_NLO_13_04",                                                     #15-19
           "mp2mp_NLO_20_04","mp2mp_NLO_21_04","mp2mp_NLO_21_04_phicut","mp2mp_NLO_24_04_mitcos","mp2mp_NLO_28_04",                             #20-24
           "mp2mp_NLO_29_04","mp2mp_NLO_07_05_big","mp2mp_NLO_07_05_small","mp2mp_NLO_08_05_full","mp2mp_NLO_08_05_full_costh3test",                  #25-29
           "mp2mp_NLO_11_05_BIG","mp2mp_NLO_12_05_BIG","mp2mp_NLO_12_05_SMALL","mp2mp_NLO_12_05_TH100MeV_BIG","mp2mp_NLO_13_05_TH100MeV_SMALL",       #30-34
           "mp2mp_NLO_13_05_TH500MeV_BIG","mp2mp_NLO_13_05_TH500MeV_SMALL","20_05_Eph100MeV_SMALL","20_05_Eph200MeV_SMALL","20_05_Eph500MeV_SMALL",   #35-39
           "20_05_Eph100MeV_BIG","20_05_Eph200MeV_BIG","20_05_Eph500MeV_BIG","20_05_100MeV_large_BIG","20_05_200MeV_large_BIG",                 #40-44
           "20_05_500MeV_large_BIG"] #45
nlo_outs = lo_outs
savenames = ["combined", "15_03", "17_03", "18_03", "23_03",    #0-4
             "24_03", "25_03","26_03","27_03","13_04",          #5-9
             "14_04","14_04_add","20_04","21_04","22_04",       #10-14
             "24_04","28_04","29_04","4_5","5_5",               #15-19
             "7_5","8_5","11_05","12_5","13_5",                 #20-24  
             "20_05","28_5"]                                    #25-26

# =========================
# Dataset choice/ Has to be checked each time!
# =========================
lo_i = 44
nlo_i = 44
savename_i = 26
nbins = 500

bin_width = 0.0382 #ECal2 with 10x cells with 38.2 mm x 38.2 mm ->active area x&y: [-19.1;19.1]
n_bands = 10
band_min = -(n_bands/2 * bin_width)
band_max = n_bands/2 * bin_width
Y5_RANGE = (band_min, band_max)
X5_RANGE = (band_min, band_max)


savename_base = savenames[savename_i] + "_" + nlo_outs[nlo_i]

# Redirect stdout
log_file = outdir_vals + f"{savename_base}_output.txt"
sys.stdout = Tee(log_file)
print("=========================","\nBIG -- 0.001  <  Q2 (GeV2/c2) < 0.04","\n=========================")
#print("=========================","\nSMALL -- 0.0005  <  Q2 (GeV2/c2) < 0.001","\n=========================")

# =========================
# Physics setup
# =========================
#print("PATH =", homedir + lo_outs[lo_i] + "/out")
#exit()
setup(folder=homedir + lo_outs[lo_i] + "/out")
lo = mergefks(sigma("mp2mp0")) * alpha**2 * conv

setup(folder=homedir + nlo_outs[nlo_i] + "/out")
nlo = (mergefks(sigma("mp2mpR"), sigma("mp2mpF"), anyxi=sigma("mp2mpA")) * alpha**3 * conv)
full = lo + nlo
onlyR = (mergefks(sigma("mp2mpR")) * alpha**3 * conv)

# =========================
# Extract observables (LAB)
# =========================
lo_th3, nlo_th3, full_th3 = lo["th3"], nlo["th3"], full["th3"]
lo_Emu, nlo_Emu, full_Emu = lo["Emu"], nlo["Emu"], full["Emu"]
lo_th5, nlo_th5, full_th5 = lo["th5"], nlo["th5"], full["th5"]
lo_Eph, nlo_Eph, full_Eph = lo["Eph"], nlo["Eph"], full["Eph"]
lo_phi5, nlo_phi5, full_phi5 = lo["phi5"], nlo["phi5"], full["phi5"]
lo_x5, nlo_x5, full_x5 = lo["x5"], nlo["x5"], full["x5"]
lo_y5, nlo_y5, full_y5 = lo["y5"], nlo["y5"], full["y5"]
lo_ql51, nlo_ql51, full_ql51 = lo["ql5(1)"], nlo["ql5(1)"], full["ql5(1)"]
lo_ql52, nlo_ql52, full_ql52 = lo["ql5(2)"], nlo["ql5(2)"], full["ql5(2)"]
lo_costh3, nlo_costh3, full_costh3 = lo["costh3"], nlo["costh3"], full["costh3"]
lo_Q2, nlo_Q2, full_Q2 = lo["Qsq"], nlo["Qsq"], full["Qsq"]

x5_bands_lo, x5_bands_nlo, x5_bands_full = {}, {}, {}
y5_bands_lo, y5_bands_nlo, y5_bands_full = {}, {}, {}


# =========================
# Fill bands robustly with try/except to avoid KeyError
# =========================
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

# =========================
# Extract observables (CMS)
# =========================
lo_th3_cms, nlo_th3_cms, full_th3_cms = lo["th3_cms"], nlo["th3_cms"], full["th3_cms"]
lo_Emu_cms, nlo_Emu_cms, full_Emu_cms = lo["Emu_cms"], nlo["Emu_cms"], full["Emu_cms"]
lo_th5_cms, nlo_th5_cms, full_th5_cms = lo["th5_cms"], nlo["th5_cms"], full["th5_cms"]
lo_Eph_cms, nlo_Eph_cms, full_Eph_cms = lo["Eph_cms"], nlo["Eph_cms"], full["Eph_cms"]
lo_phi5_cms, nlo_phi5_cms, full_phi5_cms = lo["phi5_cms"], nlo["phi5_cms"], full["phi5_cms"]


# =========================
# Colors
# =========================
colors = dict(lo="#1f77b4", nlo="#ff7f0e", full="#2ca02c", K="#d62728")
#colors = dict(lo="#2ca02c", nlo="#4fa3ff", full="#ff69b4", K="#d62728")


# =========================
# Limits
# =========================
th3vals = finite_bins(full_th3)
th3_min = np.min(th3vals[:, 0])
th3_max = np.max(th3vals[:, 0])
print("th3_min:    ",th3_min*1e3,", th3_max:    ",th3_max*1e3, " (in mrad)")
costh3vals = finite_bins(full_costh3)
costh3_min = np.min(costh3vals[:, 0])
costh3_max = np.max(costh3vals[:, 0])
"""
print(full_th3[0, 0])
print(full_costh3[0, 0])
print("")
print(full_th3[1, 0])
print(full_costh3[1, 0])
print("")
print(full_th3[2, 0])
print(full_costh3[2, 0])
#exit()
"""

print("costh3_min: ",costh3_min,", costh3_max: ",costh3_max)


# =========================
# Function to make plots & K-factors
# =========================    
def make_plots_and_kfactors( *, tag, savename_base, 
                            lo_th3, nlo_th3, full_th3, 
                            lo_Emu, nlo_Emu, full_Emu, 
                            lo_th5, nlo_th5, full_th5, 
                            lo_Eph, nlo_Eph, full_Eph, 
                            lo_phi5, nlo_phi5, full_phi5, 
                            lo_costh3, nlo_costh3, full_costh3,
                            lo_Q2, nlo_Q2, full_Q2,
                            nlo_ql51=None, nlo_ql52=None,
                            lo_x5, nlo_x5, full_x5,
                            lo_y5, nlo_y5, full_y5,
                            outdir, outdir_vals, colors, ): 
    savename = f"{savename_base}_{tag}" 

    # =========================
    # write value files for all variables
    # =========================   
    # (var, has_cms, photon_only, lab_only)

    variables = [
        ("th3",  True,  False, False),
        ("Emu",  True,  False, False),
        ("th5",  True,  True,  False),
        ("Eph",  True,  True,  False),
        ("phi5", False, True,  False),
        ("costh3", False, False,  True),
        ("Q2", False, False,  True),
        ("x5",   False, True,  False),
        ("y5",   False, True,  False),
        ("ql51", False, True,  True),
        ("ql52", False, True,  True),
    ]

    data_map = {
        "th3":  {"lo": lo_th3,  "nlo": nlo_th3,  "full": full_th3,
                "lo_cms": lo_th3_cms, "nlo_cms": nlo_th3_cms, "full_cms": full_th3_cms},
        "Emu":  {"lo": lo_Emu,  "nlo": nlo_Emu,  "full": full_Emu,
                "lo_cms": lo_Emu_cms, "nlo_cms": nlo_Emu_cms, "full_cms": full_Emu_cms},
        "th5":  {"lo": lo_th5,  "nlo": nlo_th5,  "full": full_th5,
                "lo_cms": lo_th5_cms, "nlo_cms": nlo_th5_cms, "full_cms": full_th5_cms},
        "Eph":  {"lo": lo_Eph,  "nlo": nlo_Eph,  "full": full_Eph,
                "lo_cms": lo_Eph_cms, "nlo_cms": nlo_Eph_cms, "full_cms": full_Eph_cms},
        "phi5": {"lo": lo_phi5, "nlo": nlo_phi5, "full": full_phi5},
        "costh3": {"lo": lo_costh3, "nlo": nlo_costh3, "full": full_costh3},
        "Q2": {"lo": lo_Q2, "nlo": nlo_Q2, "full": full_Q2},
        "x5":   {"lo": lo_x5,   "nlo": nlo_x5,   "full": full_x5},
        "y5":   {"lo": lo_y5,   "nlo": nlo_y5,   "full": full_y5},
        "ql51": {"nlo": nlo_ql51},
        "ql52": {"nlo": nlo_ql52},
    }

    for var, has_cms, photon_only, lab_only in variables:
        if lab_only and tag != "lab":
            continue

        orders = ["nlo"] if photon_only else ["lo", "nlo", "full"]

        # lab frame
        for order in orders:
            arr = data_map.get(var, {}).get(order, None)
            if arr is not None:
                write_file_with_values(outdir_vals + f"{order}_{var}_{savename}.txt",arr,f"{var}_{tag} bin centers","value")

        # cms frame
        if has_cms:
            for order in orders:
                arr_cms = data_map.get(var, {}).get(f"{order}_cms", None)
                if arr_cms is not None:
                    write_file_with_values(outdir_vals + f"{order}_{var}_cms_{savename}.txt",arr_cms,f"{var}_cms_{tag} bin centers","value")

                    
    # =========================
    # Create combined figure
    # =========================   
    fig, axes = create_figure(nrows=8, ncols=2, figsize=(16,22), font_size=12,sharex=False, gridspec_kw={"height_ratios":[3,1]*4,"hspace":0.6})

    ax_th3, ax_Q2 = axes[0]
    ax_K_th3, ax_K_Q2 = axes[1]
    ax_Emu, ax_Eph = axes[2]
    ax_K_Emu, ax_K_Eph = axes[3]
    ax_th5, ax_phi5  = axes[4]
    ax_K_th5, ax_K_phi5 = axes[5]
    ax_x5, ax_y5 = axes[6]
    ax_K_x5, ax_K_y5 = axes[7]

    ax_K_th3.sharex(ax_th3)
    ax_K_Emu.sharex(ax_Emu)
    ax_K_th5.sharex(ax_th5)
    ax_K_Eph.sharex(ax_Eph)
    ax_K_phi5.sharex(ax_phi5)
    ax_K_Q2.sharex(ax_Q2)
    ax_K_x5.sharex(ax_x5)
    ax_K_y5.sharex(ax_y5)

    # ---------- LAB/CMS plots ----------
    _, _, _, K_th3 = draw_observable_and_k(ax_th3, ax_K_th3,
                                           lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3,
                                           scale_factor=1e-3, x_label_main=r"$\theta_3$ (mrad)",
                                           x_label_k=r"$\theta_3$ (mrad)", y_label_main=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
                                           main_title=f"Muon Scattering Angle ({tag})", force_main_linear=False, colors=colors)

    _, _, _, K_Emu = draw_observable_and_k(ax_Emu, ax_K_Emu,
                                           lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu,
                                           scale_factor=1e3, x_label_main=r"$E_\mu$ (GeV)",
                                           x_label_k=r"$E_\mu$ (GeV)", y_label_main=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
                                           main_title=f"Muon Energy ({tag})", colors=colors)

    _, _, _, K_th5 = draw_observable_and_k(ax_th5, ax_K_th5,
                                           lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5,
                                           scale_factor=1e-3, x_label_main=r"$\theta_5$ (mrad)",
                                           x_label_k=r"$\theta_5$ (mrad)", y_label_main=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
                                           main_title=f"Photon Scattering Angle ({tag})", xlim=(-1.,13.), colors=colors) #, force_main_linear=True,

    _, _, _, K_Eph = draw_observable_and_k(ax_Eph, ax_K_Eph,
                                           lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph,
                                           scale_factor=1e3, x_label_main=r"$E_\gamma$ (GeV)",
                                           x_label_k=r"$E_\gamma$ (GeV)", y_label_main=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
                                           main_title=f"Photon Energy ({tag})", colors=colors)

    _, _, _, K_phi5 = draw_observable_and_k(ax_phi5, ax_K_phi5,
                                            lo_hist=lo_phi5, nlo_hist=nlo_phi5, full_hist=full_phi5,
                                            scale_factor=1e-3, x_label_main=r"$\phi_5$ (mrad)",
                                            x_label_k=r"$\phi_5$ (mrad)", y_label_main=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
                                            main_title=f"Photon Deflection Angle ({tag})", force_main_linear=False, colors=colors)
    _, _, _, K_Q2 = draw_observable_and_k(ax_Q2, ax_K_Q2,
                                            lo_hist=lo_Q2, nlo_hist=nlo_Q2, full_hist=full_Q2,
                                            scale_factor=1.e6,
                                            x_label_main=r"$Q^2$",
                                            x_label_k=r"$Q^2$",
                                            y_label_main=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
                                            main_title=f"Momentum Transfer({tag})",
                                            force_main_linear=True,
                                            colors=colors)

    _, _, _, K_x5 = draw_observable_and_k(ax_x5, ax_K_x5,
                                            lo_hist=lo_x5, nlo_hist=nlo_x5, full_hist=full_x5,
                                            scale_factor=1., x_label_main=r"$x_5$ (mrad)",
                                            x_label_k=r"$x_5$ (m)", y_label_main=r"$\frac{d\sigma}{dx_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
                                            main_title=f"Photon x({tag})", xlim=(-0.2,0.2), colors=colors)

    _, _, _, K_y5 = draw_observable_and_k(ax_y5, ax_K_y5,
                                            lo_hist=lo_y5, nlo_hist=nlo_y5, full_hist=full_y5,
                                            scale_factor=1., x_label_main=r"$y_5$ (mrad)",
                                            x_label_k=r"$y_5$ (m)", y_label_main=r"$\frac{d\sigma}{dy_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
                                            main_title=f"Photon y ({tag})", xlim=(-0.2,0.2), colors=colors)


    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)

    
    # x5 in y5-slices
    plot_bands(x5_bands_nlo,
            xlabel=r"$x_5\ (\mathrm{m})$",
            ylabel=r"$\frac{d\sigma}{dx_5} (\mu\mathrm{barn}/\mathrm{m})$",
            title=f"x5 distribution in y5-slices ({tag})",
            savename=f"{savename}_x5_allbands",
            outdir=outdir,
            colors=colors,
            slice_name="y5",
            slice_range=Y5_RANGE,
            yscale="log")

    # y5 in x5-slices
    plot_bands(y5_bands_nlo,
            xlabel=r"$y_5\ (\mathrm{m})$",
            ylabel=r"$\frac{d\sigma}{dy_5} (\mu\mathrm{barn}/\mathrm{m})$",
            title=f"y5 distribution in x5-slices ({tag})",
            savename=f"{savename}_y5_allbands",
            outdir=outdir,
            colors=colors,
            slice_name="x5",
            slice_range=X5_RANGE,
            yscale="log")

    # =========================
    # 2D plot: x5 vs y5
    # =========================
    keys = sorted(x5_bands_nlo.keys())

    rows = []
    for i in keys:
        band = np.array(x5_bands_nlo[i])
        if band is not None and len(band) > 0:
            vals = band[:, 1]
            n_rebin = len(vals) // n_bands
            #vals_rebinned = vals[:n_bands]
            #rows.append(vals_rebinned * bin_width)  

            #vals_rebinned = vals[:n_bands * n_rebin].reshape(n_bands, n_rebin).sum(axis=1)
            vals_rebinned = vals[:n_bands * n_rebin].reshape(n_bands, n_rebin).mean(axis=1) #dsigma/dx *Dx !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!new check
            rows.append(vals_rebinned * bin_width)

    Z = np.array(rows)  # (10, 10)
    sigma_photons_2D = np.sum(Z)
    print("\nsigma_photons from 2D distribution: ", sigma_photons_2D*1e-3, " mb")
    Rate_ECAL = calculate_rate(sigma_photons_2D*1e-3)
    print("Rate_ECAL from 2D distribution:     ", Rate_ECAL, " 1/s")

    x_centers = np.linspace(-0.191 + bin_width/2, 0.191 - bin_width/2, 10)

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(Z,extent=[band_min, band_max, band_min, band_max],origin="lower", aspect="auto", cmap="viridis",norm=LogNorm())
    grid_ticks = np.arange(band_min, band_max + bin_width, bin_width)
    grid_ticks = np.round(grid_ticks, 6)
    for t in grid_ticks:
        ax.axvline(t, linestyle="--", linewidth=0.4, alpha=0.5, color="white")
        ax.axhline(t, linestyle="--", linewidth=0.4, alpha=0.5, color="white")
    ax.set_xticks(grid_ticks, minor=True)
    ax.set_yticks(grid_ticks, minor=True)
    ax.grid(which="minor", linestyle="--", linewidth=0.4, alpha=0.5)

    cbar = plt.colorbar(im, ax=ax)
    cbar.set_label(r"$\Delta\sigma\ \text{per ECAL cell}\ (\mu\mathrm{barn})$")
    ax.set_xlabel(r"$x_5\ (\mathrm{m})$")
    ax.set_ylabel(r"$y_5\ (\mathrm{m})$")
    ax.set_title(f"2D ECAL cell distribution ({tag})")
    ax.set_xlim(-0.2, 0.2)
    ax.set_ylim(-0.2, 0.2)

    save_figure(fig, f"{savename}_x5y5_2D", outdir=outdir)
    plt.close(fig)
    '''
    # =========================
    # write K-value files
    # =========================
    write_file_with_values(outdir_vals + f"K_theta_3_{savename}.txt", K_th3, f"th3_{tag} bin center", "K_th3") 
    write_file_with_values(outdir_vals + f"K_E_mu_{savename}.txt", K_Emu, f"Emu_{tag} bin center", "K_Emu") 
    write_file_with_values(outdir_vals + f"K_theta_5_{savename}.txt", K_th5, f"th5_{tag} bin center", "K_th5") 
    write_file_with_values(outdir_vals + f"K_E_gamma_{savename}.txt", K_Eph, f"Eph_{tag} bin center", "K_Eph")
    if K_phi5 is not None:
        write_file_with_values(outdir_vals + f"K_phi5_{savename}.txt", K_phi5, f"phi5_{tag} bin center", "K_phi5")
    if K_Q2 is not None:
        write_file_with_values(outdir_vals + f"K_Q2_{savename}.txt", K_Q2, f"costh3_{tag} bin center", "K_Q2")
    if K_x5 is not None:
        write_file_with_values(outdir_vals + f"K_x5_{savename}.txt", K_x5, f"x5_{tag} bin center", "K_x5")
    if K_y5 is not None:
        write_file_with_values(outdir_vals + f"K_y5_{savename}.txt", K_y5, f"y5_{tag} bin center", "K_y5")
    '''
    # =========================
    # separate pair plots
    # =========================
    save_single_pair_plot( savename=f"{savename}_th3_pair", 
                          lo_hist=lo_th3, nlo_hist=nlo_th3, full_hist=full_th3, 
                          scale_factor=1.e-3, x_label=r"$\theta_3\ (\mathrm{mrad})$", y_label=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$", 
                          main_title=f"Muon Scattering Angle ({tag})", force_main_linear=False, colors=colors, outdir=outdir, ) 
    save_single_pair_plot( savename=f"{savename}_Emu_pair", 
                          lo_hist=lo_Emu, nlo_hist=nlo_Emu, full_hist=full_Emu, 
                          scale_factor=1.e3, x_label=r"$E_\mu\ (\mathrm{GeV})$", y_label=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$", 
                          main_title=f"Energy of the Scattered Muon ({tag})", colors=colors, outdir=outdir, ) 
    save_single_pair_plot( savename=f"{savename}_th5_pair", 
                          lo_hist=lo_th5, nlo_hist=nlo_th5, full_hist=full_th5, 
                          scale_factor=1.e-3, x_label=r"$\theta_5\ (\mathrm{mrad})$", y_label=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$", 
                          main_title=f"Photon Scattering Angle ({tag})", xlim=(-1., 13.), colors=colors, outdir=outdir, ) 
    save_single_pair_plot( savename=f"{savename}_Eph_pair", 
                          lo_hist=lo_Eph, nlo_hist=nlo_Eph, full_hist=full_Eph, 
                          scale_factor=1.e3, x_label=r"$E_\gamma\ (\mathrm{GeV})$", y_label=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$", 
                          main_title=f"Photon Energy ({tag})", colors=colors, outdir=outdir, ) 
    save_single_pair_plot( savename=f"{savename}_phi5_pair", 
                          lo_hist=lo_phi5, nlo_hist=nlo_phi5, full_hist=full_phi5, 
                          scale_factor=1.e-3, x_label=r"$\phi_5\ (\mathrm{mrad})$", y_label=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$", 
                          main_title=f"Photon X-deflection ({tag})", colors=colors, outdir=outdir, ) 
    save_single_pair_plot(savename=f"{savename}_costh3_pair",
                          lo_hist=lo_costh3, nlo_hist=nlo_costh3, full_hist=full_costh3,
                          scale_factor=1., x_label=r"$\cos\theta_3$", y_label=r"$\frac{d\sigma}{d\cos\theta_3}\ (\mu\mathrm{barn})$",
                          main_title=f"Muon Scattering Angle cos({tag})",xlim=(0.98,1.),force_main_linear=True,colors=colors, outdir=outdir,)
    save_single_pair_plot(savename=f"{savename}_Q2_pair",
                          lo_hist=lo_Q2, nlo_hist=nlo_Q2, full_hist=full_Q2,
                          scale_factor=1.e6, x_label=r"$Q^2\ (\mathrm{GeV}^2)$", y_label=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
                          main_title=f"$Q^2$({tag})",force_main_linear=False,colors=colors, outdir=outdir,)
    save_single_pair_plot( savename=f"{savename}_x5_pair", 
                          lo_hist=lo_x5, nlo_hist=nlo_x5, full_hist=full_x5, 
                          scale_factor=1., x_label=r"$x_5\ (\mathrm{m})$", y_label=r"$\frac{d\sigma}{dx_5}\ (\mu\mathrm{barn}/\mathrm{m})$", 
                          main_title=f"Photon X Hit ({tag})", colors=colors, outdir=outdir, )
    save_single_pair_plot( savename=f"{savename}_y5_pair", 
                          lo_hist=lo_y5, nlo_hist=nlo_y5, full_hist=full_y5, 
                          scale_factor=1., x_label=r"$y_5\ (\mathrm{m})$", y_label=r"$\frac{d\sigma}{dy_5}\ (\mu\mathrm{barn}/\mathrm{m})$", 
                          main_title=f"Photon Y Hit ({tag})", colors=colors, outdir=outdir, )


# =========================
# Run for LAB and CMS
# =========================
make_plots_and_kfactors(tag="lab", savename_base=savename_base,
                        lo_th3=lo_th3, nlo_th3=nlo_th3, full_th3=full_th3,
                        lo_Emu=lo_Emu, nlo_Emu=nlo_Emu, full_Emu=full_Emu,
                        lo_th5=lo_th5, nlo_th5=nlo_th5, full_th5=full_th5,
                        lo_Eph=lo_Eph, nlo_Eph=nlo_Eph, full_Eph=full_Eph,
                        lo_phi5=lo_phi5, nlo_phi5=nlo_phi5, full_phi5=full_phi5,
                        lo_costh3=lo_costh3, nlo_costh3=nlo_costh3, full_costh3=full_costh3,
                        lo_Q2=lo_Q2, nlo_Q2=nlo_Q2, full_Q2=full_Q2,
                        nlo_ql51=nlo_ql51, nlo_ql52=nlo_ql52,
                        lo_x5=lo_x5, nlo_x5=nlo_x5, full_x5=full_x5,
                        lo_y5=lo_y5, nlo_y5=nlo_y5, full_y5=full_y5,
                        outdir=outdir, outdir_vals=outdir_vals, colors=colors)

# make_plots_and_kfactors(tag="cms", savename_base=savename_base,
#                         lo_th3=lo_th3_cms, nlo_th3=nlo_th3_cms, full_th3=full_th3_cms,
#                         lo_Emu=lo_Emu_cms, nlo_Emu=nlo_Emu_cms, full_Emu=full_Emu_cms,
#                         lo_th5=lo_th5_cms, nlo_th5=nlo_th5_cms, full_th5=full_th5_cms,
#                         lo_Eph=lo_Eph_cms, nlo_Eph=nlo_Eph_cms, full_Eph=full_Eph_cms,
#                         nlo_ql51=None, nlo_ql52=None,
#                         lo_phi5=None, nlo_phi5=None, full_phi5=None,
#                         lo_x5=None, nlo_x5=None, full_x5=None,
#                         lo_y5=None, nlo_y5=None, full_y5=None,
#                         outdir=outdir, outdir_vals=outdir_vals, colors=colors)

# =========================
# Plot Numeric vs. Analytic
# =========================
plot_costh3_with_analytic(lo_hist=lo_costh3,
                          nlo_hist=nlo_costh3,
                          full_hist=full_costh3,
                          nbins=nbins,
                          th3_min=th3_min, th3_max=th3_max,
                          colors=colors,
                          savename=f"{savename_base}_costh3_analytic",
                          outdir=outdir,
                          outdir_vals=outdir_vals)

plot_Q2_with_analytic(lo_hist=lo_Q2,
                          nlo_hist=nlo_Q2,
                          full_hist=full_Q2,
                          nbins=nbins,
                          ylow_diff=-0.5,yup_diff=1.5,
                          colors=colors,
                          savename=f"{savename_base}_Q2_analytic",
                          outdir=outdir,
                          outdir_vals=outdir_vals)


# =========================
# Calculate total cross section
# =========================
lo_Q2_finite   = finite_bins(lo_Q2)
nlo_Q2_finite  = finite_bins(nlo_Q2)

sigma_lo = lo.value
sigma_nlo = nlo.value
sigma_full = full.value
sigma_Rph = onlyR.value 
sigma_lo_mb = sigma_lo /1000
sigma_nlo_mb = sigma_nlo /1000
sigma_full_mb = sigma_full /1000
sigma_Rph_mb = sigma_Rph /1000

# =========================
# Calculate Rate
# =========================
Rate = calculate_rate(sigma_lo_mb)
Rate_Rph = calculate_rate(sigma_Rph_mb)


print("\n-------------------------","\nRESULTS","\n-------------------------")
print("LO cross section:                         ", sigma_lo_mb, "mb")
print("LO cross section Paper Big Range:           0.255 mb")

print("\nLO Rate:               ", Rate,"1/s")
#print("Paper:                   86.6 1/s")
#print("Rate Paper Small Range:  89.0 1/s")

print("\n------------------------- NLO STUFF -------------------------")
print("Cross section real photon: ", sigma_Rph_mb, "mb")
print("Rate real photon:          ", Rate_Rph, "1/s")
#print("sigma_R/sigma_nlo:         ", sigma_Rph_mb/sigma_nlo_mb)

#print("\nLO cross section:   ",sigma_lo_mb, "mb")
print("NLO cross section:  ",sigma_nlo_mb, "mb")
#print("sigma_nlo/sigma_lo: ", sigma_nlo_mb/sigma_lo_mb)

print("\nsigma_full: ",sigma_full_mb,"mb")


# =========================
# Fraction of photons with E_gamma > 2 GeV
# =========================

E_cut = 2000  # MeV
onlyR_Eph = onlyR["Eph"]

Eph_finite = finite_bins(onlyR_Eph)

bin_centers = Eph_finite[:,0]
bin_values  = Eph_finite[:,1]

# average bin width
bin_width = np.mean(np.diff(bin_centers))

# total photon cross section
sigma_total = np.sum(bin_values * bin_width)

# only photons above cut
mask = bin_centers > E_cut
sigma_above = np.sum(bin_values[mask] * bin_width)

fraction = sigma_above / sigma_total
percent = fraction * 100

print("\n-------------------------")
print("onlyR.value: ",onlyR.value)
print("sigma_total: ",sigma_total)
print("sigma_above: ",sigma_above)
print(f"\nFraction of photons with Eγ > {E_cut/1000} GeV:")
print(percent, "%")



sys.stdout.close()
sys.stdout = sys.stdout.terminal  # restore normal stdout