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
outs = ["25_06_200MeV_Q2big_xi01_onlyR","26_06_200MeV_Q2big_xi1_onlyR","26_06_200MeV_Q2big_xi001_onlyR","26_06_200MeV_Q2big_xi0001_onlyR","02_07_200MeV_Q2big_xi01_onlyR_Mehran_noEmucut_correct", #0-4
        "01_07_100MeV_Q2big_xi01_onlyR","01_07_500MeV_Q2big_xi01_onlyR","08_07_200MeV_Q2big_xi01_onlyR_Mehran_precthmurange","09_07_200MeV_xi01_onlyR_Mehran","14_07_02to70GeV_xi01_onlyR_Mehran" ] #5-9

savenames = ["26_06", "29_06","01_07","03_07","07_07","09_07",           #0-5
             "14_07"] #6

# =========================
# Dataset choice / Has to be checked each time!
# =========================
outs_i = 9
savename_i = 6
nbins = 500
# HAS_BANDS = True
HAS_BANDS = False
HAS_XY = True
thmu = "NOR"
#thmu = "SMA"
#thmu = "BIG"

if thmu == "NOR":
    thmu_low = 0.3e-3
    thmu_up  = 2.e-3
if thmu == "SMA":
    thmu_low = 0.075e-3
    thmu_up  = 0.5e-3
if thmu == "BIG":
    thmu_low = 1.2e-3
    thmu_up  = 8.e-3


bin_width = 0.0382  # ECal2 with 10x cells with 38.2 mm x 38.2 mm -> active area x&y: [-.191;.191]
n_bands = 10
band_min = -(n_bands / 2 * bin_width)
band_max =   n_bands / 2 * bin_width
Y5_RANGE = (band_min, band_max)
X5_RANGE = (band_min, band_max)

savename_base = savenames[savename_i] + "_" + outs[outs_i]

# Redirect stdout
log_file = outdir_vals + f"{savename_base}_output.txt"
sys.stdout = Tee(log_file)
print("=======================================================================")
print("        Analysed file: ", savename_base)
print("=======================================================================")

# =========================
# Physics setup
# =========================

setup(folder=homedir + outs[outs_i] + "/out")
onlyR  = mergefks(sigma("mp2mpR")) * alpha**3 * conv #alpha=0.0072973525692838015, conv = 389379365.556916 pb


# =========================
# Extract observables (LAB)
# =========================
onlyR_th3 = onlyR["th3"]
onlyR_Emu = onlyR["Emu"]
onlyR_th5 = onlyR["th5"]
onlyR_Eph = onlyR["Eph"]
onlyR_phi5 = onlyR["phi5"]
onlyR_costh3 = onlyR["costh3"]
onlyR_Q2 = onlyR["Qsq"]

if HAS_XY:
    onlyR_x5  = onlyR["x5"]
    onlyR_y5 = onlyR["y5"]

if HAS_BANDS:
    x5_bands_onlyR, y5_bands_onlyR = {}, {}
    for i in range(1, n_bands + 1):
        try:
            x5_bands_onlyR[i] = onlyR[f"x5_B{i}"]
        except KeyError:
            pass
        try:
            y5_bands_onlyR[i] = onlyR[f"y5_B{i}"]
        except KeyError:
            pass

# =========================
# Colors
# =========================
colors = dict(onlyR="#1f77b4", K="#d62728")

# =========================
# Limits
# =========================
th3vals = finite_bins(onlyR_th3)
th3vals = th3vals[th3vals[:, 1] != 0]
th3_min = np.min(th3vals[:, 0])
th3_max = np.max(th3vals[:, 0])
print("th3_min:    ", th3_min * 1e3, ", th3_max:   ", th3_max * 1e3, "(mrad)")

costh3vals = finite_bins(onlyR_costh3)
costh3vals = costh3vals[costh3vals[:, 1] != 0]
costh3_min = np.min(costh3vals[:, 0])
costh3_max = np.max(costh3vals[:, 0])
print("costh3_min: ", costh3_min, ", costh3_max: ", costh3_max)


# =========================
# Function to make plots & K-factors
# =========================
def make_plots(*, tag, savename_base,
                             onlyR_th3,
                             onlyR_Emu,
                             onlyR_th5,
                             onlyR_Eph,
                             onlyR_phi5,
                             onlyR_costh3,
                             onlyR_Q2,
                             onlyR_x5=None,
                             onlyR_y5=None,
                             outdir, colors):
    savename = f"{savename_base}_{tag}"

    # =========================
    # Combined figure
    # =========================
    has_xy = (onlyR_x5 is not None and onlyR_y5 is not None)

    nrows = 4 if has_xy else 3

    fig, axes = create_figure(
        nrows=nrows,
        ncols=2,
        figsize=(16, 3 * nrows + 2),
        font_size=12,
        sharex=False,
        gridspec_kw={"hspace": 0.6}
    )

    ax_th3, ax_Q2   = axes[0]
    ax_Emu, ax_Eph  = axes[1]
    ax_th5, ax_phi5 = axes[2]

    if has_xy:
        ax_x5, ax_y5 = axes[3]

    draw_observable_onlyR(
        ax_th3,
        onlyR_hist=onlyR_th3,
        scale_factor=1e-3,
        x_label=r"$\theta_3$ (mrad)",
        y_label=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        title=f"Muon Scattering Angle ({tag})",
        xlim=(0, th3_max * 1e3 + 0.5),  # Add a small margin to the right
        force_linear=False,
        colors=colors,
    )

    draw_observable_onlyR(
        ax_Emu,
        onlyR_hist=onlyR_Emu,
        scale_factor=1e3, 
        x_label=r"$E_\mu$ (GeV)",
        y_label=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        title=f"Muon Energy ({tag})", 
        colors=colors)

    draw_observable_onlyR(ax_th5,
        onlyR_hist=onlyR_th5,
        scale_factor=1e-3, x_label=r"$\theta_5$ (mrad)", y_label=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        title=f"Photon Scattering Angle ({tag})", xlim=(-1., 13.), colors=colors)

    draw_observable_onlyR(ax_Eph,
        onlyR_hist=onlyR_Eph,
        scale_factor=1e3, x_label=r"$E_\gamma$ (GeV)", y_label=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
        title=f"Photon Energy ({tag})", colors=colors)

    draw_observable_onlyR(ax_phi5,  
        onlyR_hist=onlyR_phi5,
        scale_factor=1e-3, x_label=r"$\phi_5$ (mrad)", y_label=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
        title=f"Photon Deflection Angle ({tag})", force_linear=False, colors=colors)

    draw_observable_onlyR(ax_Q2,
        onlyR_hist=onlyR_Q2,
        scale_factor=1.e6, x_label=r"$Q^2$", y_label=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
        title=f"Momentum Transfer ({tag})", force_linear=True, colors=colors)

    if has_xy:
        draw_observable_onlyR(ax_x5,
            onlyR_hist=onlyR_x5,
            scale_factor=1., x_label=r"$x_5$ (m)", y_label=r"$\frac{d\sigma}{dx_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
            title=f"Photon x ({tag})", xlim=(-0.2, 0.2), colors=colors)

        draw_observable_onlyR(ax_y5,
            onlyR_hist=onlyR_y5,
            scale_factor=1., x_label=r"$y_5$ (m)", y_label=r"$\frac{d\sigma}{dy_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
            title=f"Photon y ({tag})", xlim=(-0.2, 0.2), colors=colors)

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)

    # =========================
    # Separate pair plots
    # =========================
    # save_single_plot_onlyR(savename=f"{savename}_th3",
    #     onlyR_hist=onlyR_th3,
    #     scale_factor=1.e-3, x_label=r"$\theta_3\ (\mathrm{mrad})$",
    #     y_label=r"$\frac{d\sigma}{d\theta_3}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
    #     title=f"Muon Scattering Angle ({tag})",
    #     xlim=(0,thmu_up*1e3 + 0.5)  ,force_linear=False, colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_Emu",
    #     onlyR_hist=onlyR_Emu,
    #     scale_factor=1.e3, x_label=r"$E_\mu\ (\mathrm{GeV})$",
    #     y_label=r"$\frac{d\sigma}{dE_\mu}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
    #     title=f"Energy of the Scattered Muon ({tag})", colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_th5",
    #     onlyR_hist=onlyR_th5,
    #     scale_factor=1.e-3, x_label=r"$\theta_5\ (\mathrm{mrad})$",
    #     y_label=r"$\frac{d\sigma}{d\theta_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
    #     title=f"Photon Scattering Angle ({tag})", xlim=(-1., 13.), colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_Eph",
    #     onlyR_hist=onlyR_Eph,
    #     scale_factor=1.e3, x_label=r"$E_\gamma\ (\mathrm{GeV})$",
    #     y_label=r"$\frac{d\sigma}{dE_\gamma}\ (\mu\mathrm{barn}/\mathrm{GeV})$",
    #     title=f"Photon Energy ({tag})", colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_phi5",
    #     onlyR_hist=onlyR_phi5,
    #     scale_factor=1.e-3, x_label=r"$\phi_5\ (\mathrm{mrad})$",
    #     y_label=r"$\frac{d\sigma}{d\phi_5}\ (\mu\mathrm{barn}/\mathrm{mrad})$",
    #     title=f"Photon X-deflection ({tag})", colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_costh3",
    #     onlyR_hist=onlyR_costh3,
    #     scale_factor=1., x_label=r"$\cos\theta_3$",
    #     y_label=r"$\frac{d\sigma}{d\cos\theta_3}\ (\mu\mathrm{barn})$",
    #     title=f"Muon Scattering Angle cos ({tag})", xlim=(0.99, 1.),
    #     force_linear=True, colors=colors, outdir=outdir)

    # save_single_plot_onlyR(savename=f"{savename}_Q2",
    #     onlyR_hist=onlyR_Q2,
    #     scale_factor=1.e6, x_label=r"$Q^2\ (\mathrm{GeV}^2)$",
    #     y_label=r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn})$",
    #     title=f"$Q^2$ ({tag})", force_linear=False, colors=colors, outdir=outdir)

    # if has_xy:
    #     save_single_plot_onlyR(savename=f"{savename}_x5",
    #         onlyR_hist=onlyR_x5,
    #         scale_factor=1., x_label=r"$x_5\ (\mathrm{m})$",
    #         y_label=r"$\frac{d\sigma}{dx_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
    #         title=f"Photon X Hit ({tag})", colors=colors, outdir=outdir)

    #     save_single_plot_onlyR(savename=f"{savename}_y5",
    #         onlyR_hist=onlyR_y5,
    #         scale_factor=1., x_label=r"$y_5\ (\mathrm{m})$",
    #         y_label=r"$\frac{d\sigma}{dy_5}\ (\mu\mathrm{barn}/\mathrm{m})$",
    #         title=f"Photon Y Hit ({tag})", colors=colors, outdir=outdir)


def make_band_plots(*, tag, savename_base,
                    x5_bands_onlyR, y5_bands_onlyR,
                    n_bands, bin_width,
                    X5_RANGE, Y5_RANGE,
                    outdir, colors):
    """Band-Plots und 2D ECAL-Verteilung. Nur aufrufen wenn der Run x5/y5-Bands enthält."""
    savename = f"{savename_base}_{tag}"

    # x5 in y5-Slices
    plot_bands(x5_bands_onlyR,
               xlabel=r"$x_5\ (\mathrm{m})$",
               ylabel=r"$\frac{d\sigma}{dx_5} (\mu\mathrm{barn}/\mathrm{m})$",
               title=f"x5 distribution in y5-slices ({tag})",
               savename=f"{savename}_x5_allbands",
               outdir=outdir, colors=colors,
               slice_name="y5", slice_range=Y5_RANGE, yscale="log")

    # y5 in x5-Slices
    plot_bands(y5_bands_onlyR,
               xlabel=r"$y_5\ (\mathrm{m})$",
               ylabel=r"$\frac{d\sigma}{dy_5} (\mu\mathrm{barn}/\mathrm{m})$",
               title=f"y5 distribution in x5-slices ({tag})",
               savename=f"{savename}_y5_allbands",
               outdir=outdir, colors=colors,
               slice_name="x5", slice_range=X5_RANGE, yscale="log")

    # 2D ECAL-Verteilung
    keys = sorted(x5_bands_onlyR.keys())
    rows = []
    for i in keys:
        band = np.array(x5_bands_onlyR[i])
        if band is not None and len(band) > 0:
            vals = band[:, 1]
            n_rebin = len(vals) // n_bands
            vals_rebinned = vals[:n_bands * n_rebin].reshape(n_bands, n_rebin).mean(axis=1)
            rows.append(vals_rebinned * bin_width)

    Z = np.array(rows)
    sigma_photons_2D = np.sum(Z)
    print("\nPhoton cross section from 2D distribution: ", sigma_photons_2D * 1e-3, " mb")
    print("Photon rate from 2D distribution:          ", calculate_rate(sigma_photons_2D * 1e-3), " 1/s")

    band_min = -(n_bands / 2 * bin_width)
    band_max =   n_bands / 2 * bin_width

    fig, ax = plt.subplots(figsize=(6, 5))
    im = ax.imshow(Z, extent=[band_min, band_max, band_min, band_max],
                   origin="lower", aspect="auto", cmap="viridis", norm=LogNorm())
    grid_ticks = np.round(np.arange(band_min, band_max + bin_width, bin_width), 6)
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


# =========================
# Run for LAB
# =========================
make_plots(tag="lab", savename_base=savename_base,
                        onlyR_th3=onlyR_th3,
                        onlyR_Emu=onlyR_Emu,
                        onlyR_th5=onlyR_th5,
                        onlyR_Eph=onlyR_Eph,
                        onlyR_phi5=onlyR_phi5,
                        onlyR_costh3=onlyR_costh3,
                        onlyR_Q2=onlyR_Q2,
                        onlyR_x5=onlyR_x5 if HAS_XY else None,
                        onlyR_y5=onlyR_y5 if HAS_XY else None,
                        outdir=outdir, colors=colors)

if HAS_BANDS:
    make_band_plots(tag="lab", savename_base=savename_base,
                    x5_bands_onlyR=x5_bands_onlyR, y5_bands_onlyR=y5_bands_onlyR,
                    n_bands=n_bands, bin_width=bin_width,
                    X5_RANGE=X5_RANGE, Y5_RANGE=Y5_RANGE,
                    outdir=outdir, colors=colors)



sigma_Rph   = onlyR.value
sigma_Rph_mb  = onlyR.value[0] / 1000

# =========================
# Calculate Rate
# =========================
Rate_Rph = calculate_rate(sigma_Rph_mb)

print("\n-------------------------", "\nRESULTS", "\n-------------------------")
print("Real photon cross section:   ", sigma_Rph_mb, "mb")
print("\nRate real photon:            ", Rate_Rph, "1/s")
print("\n-------------------------------------------------------------")

# =========================
# Fraction of photons with E_gamma > E_cut
# =========================
E_cut = 2000  # MeV
onlyR_Eph = onlyR["Eph"]
Eph_finite = finite_bins(onlyR_Eph)

bin_centers = Eph_finite[:, 0]
bin_values  = Eph_finite[:, 1]
bin_width_Eph = np.mean(np.diff(bin_centers))  # local var, doesn't overwrite ECAL bin_width

mask        = bin_centers > E_cut
sigma_above = np.sum(bin_values[mask] * bin_width_Eph)
sigma_above_mb = sigma_above / 1000
onlyR_mb       = onlyR.value[0] / 1000

percent  = sigma_above_mb / onlyR_mb * 100   # = sigma_above / onlyR.value[0] * 100, mb kürzen sich weg

if percent > 100 or percent < 0:
    print(f"Photon fraction with Eγ > {E_cut/1000} GeV:       unreliable (FKS bin artefacts)")
    print(f"Photon cross section with Eγ > {E_cut/1000} GeV:  unreliable (FKS bin artefacts)")
else:
    print(f"Photon cross section with Eγ > {E_cut/1000} GeV:  {sigma_above_mb:.6e} mb")
    print(f"Photon fraction with Eγ > {E_cut/1000} GeV:       {percent:.2f} %")

sys.stdout.close()
sys.stdout = sys.stdout.terminal  # restore normal stdout
