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
homedir = "/home/marialei/AMBER_RadCor/" # Laptop
#homedir = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/"  # Office
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
           "20_05_500MeV_large_BIG","28_05_200MeV_Q2rel_fixedcuts","02_06_0MeV_Q2big_1601"] #45-47
nlo_outs = lo_outs
savenames = ["combined", "15_03", "17_03", "18_03", "23_03",    #0-4
             "24_03", "25_03","26_03","27_03","13_04",          #5-9
             "14_04","14_04_add","20_04","21_04","22_04",       #10-14
             "24_04","28_04","29_04","4_5","5_5",               #15-19
             "7_5","8_5","11_05","12_5","13_5",                 #20-24  
             "20_05","28_5","29_5","01_06","02_06"]             #25-29

# =========================
# Dataset choice/ Has to be checked each time!
# =========================
lo_i = 47
nlo_i = 47
savename_i = 29
nbins = 500

bin_width = 0.0382 #ECal2 with 10x cells with 38.2 mm x 38.2 mm ->active area x&y: [-.191;19.1]
n_bands = 10
band_min = -(n_bands/2 * bin_width)
band_max = n_bands/2 * bin_width
Y5_RANGE = (band_min, band_max)
X5_RANGE = (band_min, band_max)


savename_base = savenames[savename_i] + "_" + nlo_outs[nlo_i]

# Redirect stdout
log_file = outdir_vals + f"{savename_base}_output.txt"
sys.stdout = Tee(log_file)
print(savename_base)
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
lo_costh3, nlo_costh3, full_costh3 = lo["costh3"], nlo["costh3"], full["costh3"]
lo_Q2, nlo_Q2, full_Q2 = lo["Qsq"], nlo["Qsq"], full["Qsq"]


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
                            outdir, outdir_vals, colors, ): 
    savename = f"{savename_base}_{tag}" 

    # =========================
    # write value files for all variables
    # =========================   
    # (var, photon_only, lab_only)

    variables = [
        ("th3",  False, False),
        ("Emu",  False, False),
        ("th5",  True,  False),
        ("Eph",  True,  False),
        ("phi5", True,  False),
        ("costh3", False,  True),
        ("Q2", False,  True),
    ]

    data_map = {
        "th3":  {"lo": lo_th3,  "nlo": nlo_th3,  "full": full_th3},
        "Emu":  {"lo": lo_Emu,  "nlo": nlo_Emu,  "full": full_Emu},
        "th5":  {"lo": lo_th5,  "nlo": nlo_th5,  "full": full_th5},
        "Eph":  {"lo": lo_Eph,  "nlo": nlo_Eph,  "full": full_Eph},
        "phi5": {"lo": lo_phi5, "nlo": nlo_phi5, "full": full_phi5},
        "costh3": {"lo": lo_costh3, "nlo": nlo_costh3, "full": full_costh3},
        "Q2": {"lo": lo_Q2, "nlo": nlo_Q2, "full": full_Q2},
    }

    for var, photon_only, lab_only in variables:
        if lab_only and tag != "lab":
            continue

        orders = ["nlo"] if photon_only else ["lo", "nlo", "full"]

        # lab frame
        for order in orders:
            arr = data_map.get(var, {}).get(order, None)
            if arr is not None:
                write_file_with_values(outdir_vals + f"{order}_{var}_{savename}.txt",arr,f"{var}_{tag} bin centers","value")

                    
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

    ax_K_th3.sharex(ax_th3)
    ax_K_Emu.sharex(ax_Emu)
    ax_K_th5.sharex(ax_th5)
    ax_K_Eph.sharex(ax_Eph)
    ax_K_phi5.sharex(ax_phi5)
    ax_K_Q2.sharex(ax_Q2)

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

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)


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
                        outdir=outdir, outdir_vals=outdir_vals, colors=colors)

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