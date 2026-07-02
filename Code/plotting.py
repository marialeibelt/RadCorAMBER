import matplotlib.pyplot as plt
import numpy as np
from pathlib import Path
from matplotlib.ticker import ScalarFormatter
from pymule import *
from theo_calc import *

# =========================
# Basics
# =========================
def finite_bins(hist):
    """Entfernt Bins mit nicht-finiten x-Werten (McMule underflow/overflow)."""
    mask = np.isfinite(hist[:, 0])
    return hist[mask]


def write_file_with_values(filename, parameter_array, xlabel, parameter_label):
    with open(filename, "w") as f:
        f.write(f"{xlabel}    {parameter_label}    Error \n")
        for i, row in enumerate(parameter_array):
            bin_center = row[0]
            value      = row[1]
            error      = row[2]
            f.write(f"{bin_center: .6e}  {value: .6e}  {error: .6e}\n")


def integrate_hist(hist, xmin=None, xmax=None):

    hist = np.array(hist)
    hist = finite_bins(hist)

    centers = hist[:, 0]
    values  = hist[:, 1]

    widths = np.diff(centers)
    widths = np.append(widths, widths[-1])

    mask = np.ones_like(centers, dtype=bool)

    if xmin is not None:
        mask &= (centers >= xmin)

    if xmax is not None:
        mask &= (centers <= xmax)

    sigma = np.sum(values[mask] * widths[mask])

    return sigma

def calculate_rate(sigma_mb):
    ltarget = 80 #cm
    Hpres = 20 #bar
    Npvol = 2.687*1e19 #Protons/cm^3
    Ibeam = 2*1e6 #1/s

    Np = Npvol * 2 * (Hpres/1.013) * ltarget # #Protons/cm^2 = target thickness
    rate = sigma_mb * 1e-27 * Np * Ibeam
    return rate #1/s


# =========================
# Figure + saving
# =========================
def create_figure(nrows=1, ncols=1, figsize=(8,5), font_size=7, sharex=False, gridspec_kw=None):
    plt.rcParams.update({
        "font.size": font_size,
        "font.family": "serif",
        "text.usetex": True,
    })

    fig, axes = plt.subplots(nrows, ncols, figsize=figsize, sharex=sharex, gridspec_kw=gridspec_kw)

    if nrows == 1 and ncols == 1:
        axes = np.array([[axes]])
    elif nrows == 1:
        axes = np.array([axes])
    elif ncols == 1:
        axes = np.array([[ax] for ax in axes])

    return fig, axes


def save_figure(fig, name, outdir, dpi=700):
    outdir = Path(outdir)
    outdir.mkdir(parents=True, exist_ok=True)
    #fig.savefig(outdir / f"{name}.pdf", bbox_inches="tight")
    fig.savefig(outdir / f"{name}.png", dpi=dpi, bbox_inches="tight")


# =========================
# Styling
# =========================
def style_axis(ax, xlabel=None, ylabel=None, title=None, xscale="linear", yscale="linear", legend=False):
    if xlabel:
        ax.set_xlabel(xlabel)
    if ylabel:
        ax.set_ylabel(ylabel)
    if title:
        ax.set_title(title)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)
    ax.grid(True, alpha=0.7, ls="dotted")

    if legend:
        ax.legend(framealpha=0)


def style_sci_x(ax, xlabel, ylabel, title, yscale="log"):
    style_axis(ax, xlabel=xlabel, ylabel=ylabel, title=title, yscale=yscale, legend=False)


def add_secondary_xaxis(ax, transform_func, xlabel=None, scale="linear"):
    ax_sec = ax.twiny()
    ax_sec.set_xlim(transform_func(np.array(ax.get_xlim())))
    ax_sec.set_xscale(scale)
    if xlabel:
        ax_sec.set_xlabel(xlabel, fontsize=9)
    return ax_sec


# =========================
# Plotting
# =========================
def plot_errorband(ax, hist, color):
    plt.sca(ax)
    for line in errorband(hist):
        line.set_color(color)
        line.set_alpha(0.6) #Macht Linien bisschen durchsichtig
        #line.set_linewidth(1.2)


def plot_lo_nlo_full(ax, lo, nlo, full, colors, labels=None):
    if lo is not None:
        plot_errorband(ax, lo, colors["lo"])
    if nlo is not None:
        plot_errorband(ax, nlo, colors["nlo"])
    if full is not None:
        plot_errorband(ax, full, colors["full"])

    if labels is not None:
        for key, label in labels.items():
            ax.plot([], [], color=colors[key], label=label)
        ax.legend(framealpha=0)


def plot_bands(bands_dict, *,
               xlabel, ylabel, title, savename, outdir, colors,
               slice_name="y5", slice_range=(-0.09, 0.09),
               yscale="linear", xscale="linear"):

    fig, ax = plt.subplots(figsize=(6, 5))

    if not bands_dict:
        ax.set_xlabel(xlabel)
        ax.set_ylabel(ylabel)
        ax.set_title(title)
        ax.set_xscale(xscale)
        ax.set_yscale(yscale)
        save_figure(fig, savename, outdir=outdir)
        plt.close(fig)
        return

    band_indices = sorted(bands_dict.keys())
    n_bands = len(band_indices)
    min_val, max_val = slice_range
    edges = np.linspace(min_val, max_val, n_bands + 1)
    cmap = plt.cm.viridis

    for idx, i in enumerate(band_indices):
        band = bands_dict[i]
        if band is None or len(band) == 0:
            continue

        band = np.array(band, copy=True)

        mask = np.isfinite(band[:, 0]) & np.isfinite(band[:, 1])
        band = band[mask]
        if len(band) == 0:
            continue

        color = cmap(idx / max(n_bands - 1, 1))
        low, high = edges[idx], edges[idx + 1]
        label = f"{low:.3f} < {slice_name} < {high:.3f}"

        ax.plot(band[:, 0], band[:, 1], color=color, label=label)

    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_title(title)

    ax.set_xscale(xscale)
    ax.set_yscale(yscale)

    ax.legend(fontsize=8, ncol=2, framealpha=0)

    save_figure(fig, savename, outdir=outdir)
    plt.close(fig)


def plot_K(ax, K, colors_K, xlabel):
    plt.sca(ax)
    for line in errorband(K):
        line.set_color(colors_K)
    ax.plot([], [], color=colors_K, label="K-factor")
    style_sci_x(ax, xlabel, r"$K = \mathrm{NLO}/(\mathrm{LO+NLO})$", None, yscale="linear")
    ax.legend(framealpha=0)


def draw_observable_and_k(ax_main, ax_k, *, lo_hist, nlo_hist, full_hist,
                          scale_factor, x_label_main, x_label_k, y_label_main,
                          main_title, xlim=None, ylim=None, main_yscale="log", force_main_linear=False,
                          colors=None, hide_main_xticks=True):

    # Wenn keine Daten übergeben wurden, alles ausblenden und zurückgeben
    if nlo_hist is None:
        ax_main.set_visible(False)
        ax_k.set_visible(False)
        return None, None, None, None

    lo_s   = scaleplot(lo_hist, scale_factor) if lo_hist is not None else None
    nlo_s  = scaleplot(nlo_hist, scale_factor)
    full_s = scaleplot(full_hist, scale_factor) if full_hist is not None else None

    # Plotten
    plot_lo_nlo_full(ax_main, lo_s, nlo_s, full_s, colors,
                    labels=dict(
                        **({"lo": "LO"} if lo_s is not None else {}),
                        nlo="NLO",
                        **({"full": "LO + NLO"} if full_s is not None else {})
                    ))

    yscale_to_use = main_yscale
    if force_main_linear:
        yscale_to_use = "linear"
    else:
        if full_s is not None and np.all(full_s[:, 1] <= 0):
            yscale_to_use = "linear"

    style_sci_x(ax_main, x_label_main, y_label_main, main_title, yscale=yscale_to_use)
    if xlim is not None:
        ax_main.set_xlim(*xlim)

    if ylim is not None:
        ax_k.set_ylim(*ylim)

    if full_s is not None:
        K = mergebins(divideplots(nlo_s, full_s), 5)
        plot_K(ax_k, K, colors["K"], x_label_k)
    else:
        K = None
        ax_k.set_visible(False)

    if xlim is not None:
        ax_k.set_xlim(*xlim)
    else:
        ax_k.set_xlim(ax_main.get_xlim())

    if hide_main_xticks:
        ax_main.tick_params(axis="x", labelbottom=False)
        ax_main.set_xlabel(None)

    ax_k.tick_params(axis="x", labelbottom=True)
    return lo_s, nlo_s, full_s, K    

def draw_observable_onlyR(ax, *, onlyR_hist,
                          scale_factor, x_label, y_label,
                          title, xlim=None, ylim=None,
                          yscale="log", force_linear=False,
                          colors=None):

    # Keine Daten → Achse ausblenden
    if onlyR_hist is None:
        ax.set_visible(False)
        return None

    # Histogramm skalieren
    onlyR_s = scaleplot(onlyR_hist, scale_factor)

    # Farbe
    plot_color = colors.get("onlyR", colors.get("full", "black")) if colors else "black"

    # Plot
    plot_errorband(ax, onlyR_s, plot_color)

    # y-Skalierung
    yscale_to_use = "linear" if force_linear else yscale
    if not force_linear and np.all(onlyR_s[:, 1] <= 0):
        yscale_to_use = "linear"

    style_sci_x(ax, x_label, y_label, title, yscale=yscale_to_use)

    if xlim is not None:
        ax.set_xlim(*xlim)

    if ylim is not None:
        ax.set_ylim(*ylim)

    return onlyR_s


def save_single_pair_plot( *, savename, lo_hist, nlo_hist, full_hist, 
                          scale_factor, x_label, y_label, main_title, 
                          xlim=None, ylim=None, main_yscale="log", force_main_linear=False, 
                          colors=None, outdir=None, ): 
    fig, axes = create_figure( nrows=2, ncols=1, figsize=(7, 6), 
                              font_size=12, sharex=True, gridspec_kw={ "height_ratios": [3, 1], "hspace": 0., }, 
                              ) 
    ax_main = axes[0, 0] 
    ax_k = axes[1, 0] 
    draw_observable_and_k( ax_main, ax_k, lo_hist=lo_hist, nlo_hist=nlo_hist, full_hist=full_hist,
                           scale_factor=scale_factor, x_label_main=None, x_label_k=x_label, y_label_main= y_label, main_title=main_title, 
                           xlim=xlim,ylim=ylim, main_yscale=main_yscale, force_main_linear=force_main_linear, colors=colors, ) 
    mulify(fig, delx=4.5, dely=1.)
    save_figure(fig, savename, outdir=outdir) 
    plt.close(fig) 

def save_single_plot_onlyR( *, savename, onlyR_hist, 
                          scale_factor, x_label, y_label, title, 
                          xlim=None, ylim=None, yscale="log", force_linear=False, 
                          colors=None, outdir=None, ): 
    # Da kein K-Faktor mehr geplottet wird, reicht 1 Zeile (nrows=1)
    fig, axes = create_figure( nrows=1, ncols=1, figsize=(7, 4.5), 
                              font_size=12,
                              ) 
    # Je nachdem, wie create_figure Achsen zurückgibt (2D-Array oder einzelne Achse bei nrows=1):
    # Wenn create_figure bei 1x1 ein 2D-Array liefert: axes[0, 0]. Wenn es ein einzelnes Ax-Objekt liefert: axes.
    # Hier wird angenommen, dass es wie bei subplots ein einzelnes Objekt ist, falls nicht, wieder axes[0, 0] nutzen.
    ax = axes[0, 0] if hasattr(axes, "shape") else axes
    
    # Aufruf der neuen Version von draw_observable_and_k
    draw_observable_onlyR( ax, onlyR_hist=onlyR_hist,
                           scale_factor=scale_factor, x_label=x_label, y_label=y_label, title=title, 
                           xlim=xlim, ylim=ylim, yscale=yscale, force_linear=force_linear, colors=colors)
    
    mulify(fig, delx=4.5, dely=1.)
    save_figure(fig, savename, outdir=outdir) 
    plt.close(fig)

def plot_costh3_with_analytic(lo_hist, nlo_hist, full_hist, nbins, th3_min,th3_max,colors, savename, outdir, outdir_vals):

    fig_main, axes = create_figure(nrows=2, ncols=1, figsize=(7,6),gridspec_kw={"height_ratios":[3,1], "hspace":0})
    ax_main = axes[0,0]
    ax_k    = axes[1,0]

    # -------------------------
    # num scaling + plotting
    # -------------------------
    lo_s   = scaleplot(lo_hist, 1.0)
    nlo_s  = scaleplot(nlo_hist, 1.0)
    full_s = scaleplot(full_hist, 1.0)

    plot_lo_nlo_full(ax_main, lo_s, nlo_s, full_s, colors,labels={"lo":"LO", "nlo":"NLO", "full":"LO+NLO"})

    # -------------------------
    # analytic curve
    # -------------------------
    theta_grid = np.linspace(th3_min, th3_max, nbins)
    costh_grid = np.cos(theta_grid)

    dsig_grid = np.array([dsigma_dcosth(t) for t in theta_grid])

    ax_main.plot(costh_grid, dsig_grid,color="black", linestyle="--", label="analytic")
    ax_main.legend()
    #ax_main.set_xlim(0.9999979,0.9999991)
    #ax_main.set_ylim(1e4,1.3*1e7)

    style_sci_x(ax_main,r"$\cos\theta_3$",r"$\frac{d\sigma}{d\cos\theta_3}\ (\mu\mathrm{barn})$","Muon Scattering Angle (lab)",yscale="linear")


    # -------------------------
    # Bin-by-bin Vergleich Numerik vs. Analytik
    # -------------------------
    num = finite_bins(scaleplot(lo_hist, 1.0))

    x_num = num[:, 0]  # das sind cos(theta)-Werte
    y_num = num[:, 1]
    err_num = num[:, 2]

    theta_at_bins = np.arccos(x_num)
    theo_at_bins = np.array([dsigma_dcosth(t) for t in theta_at_bins])

    diff = y_num - theo_at_bins
    rel_diff = 100 * diff / np.where(np.abs(theo_at_bins) > 1e-20, theo_at_bins, np.nan)

    out = np.column_stack([x_num, y_num, err_num, theo_at_bins, diff, rel_diff])
    np.savetxt(outdir_vals + "costh3.csv",out,delimiter=",",header="cos_theta,num,num_err,ana,diff,rel_diff_percent",comments="")

    fig_diff, ax = plt.subplots()
    ax.errorbar(x_num, rel_diff, yerr=100*err_num/theo_at_bins, fmt='o', markersize=2)
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(r"$\cos\theta_3$")
    ax.set_ylabel(r"$(num - ana)/\mathrm{ana}\ (\%)$")
    ax.set_title("Relative difference numeric vs. analytic (LO)")
    save_figure(fig_diff, f"{savename}_rel_diff", outdir=outdir)
    plt.close(fig_diff)


    K = mergebins(divideplots(nlo_s, full_s), 5)
    plot_K(ax_k, K, colors["K"], r"$\cos\theta_3$")

    ax_k.set_xlim(ax_main.get_xlim())

    mulify(fig_main, delx=4.5, dely=1.)
    save_figure(fig_main, savename, outdir=outdir)
    plt.close(fig_main)


# =========================
# Plot cos(th3) analytical vs. numerical
# =========================    
def plot_Q2_with_analytic(lo_hist, nlo_hist, full_hist,nbins,ylow_diff,yup_diff, colors, savename, outdir, outdir_vals):
    fig_main, axes = create_figure(nrows=2, ncols=1, figsize=(7,6), gridspec_kw={"height_ratios":[3,1], "hspace":0})
    ax_main = axes[0,0]
    ax_k    = axes[1,0]

    # -------------------------
    # num scaling + plotting
    # µbarn/MeV² -> µbarn/GeV² (factor 1e6)
    # -------------------------
    lo_s   = scaleplot(lo_hist,   1.e6)
    nlo_s  = scaleplot(nlo_hist,  1.e6)
    full_s = scaleplot(full_hist, 1.e6)

    plot_lo_nlo_full(ax_main, lo_s, nlo_s, full_s, colors, labels={"lo":"LO", "nlo":"NLO", "full":"LO+NLO"})

    # -------------------------
    # Q² range from histogram (bin centers still in MeV², convert to GeV²)
    # -------------------------
    num = finite_bins(scaleplot(lo_hist, 1.))
    x_num     = num[:, 0]  *1e-6  # MeV² -> GeV²
    y_num     = num[:, 1]  *1e6  # 1/MeV² -> 1/GeV²
    err_num   = num[:, 2]  *1e6  # 1/MeV² -> 1/GeV²

    # -------------------------
    # analytic curve over correct Q² range in GeV²
    # -------------------------
    q2_grid  = np.linspace(x_num[0], x_num[-1], nbins)
    dsig_grid = np.array([float(dsigma_dQ2(t)) for t in q2_grid])  # µbarn/GeV²

    ax_main.plot(q2_grid, dsig_grid, color="black", linestyle="--", label="analytic")
    ax_main.set_ylim(1e-3,1*1e6)
    ax_main.legend()

    style_sci_x(ax_main,
                r"$Q^2\ (\mathrm{GeV}^2)$",
                r"$\frac{d\sigma}{dQ^2}\ (\mu\mathrm{barn}/\mathrm{GeV}^2)$",
                r"$Q^2$ distribution (lab)",
                yscale="log")

    # -------------------------
    # Bin-by-bin comparison numeric vs. analytic
    # -------------------------
    theo_at_bins = np.array([dsigma_dQ2(t) for t in x_num])     
    #print("y_num: ",y_num)  
    #print("theo_at_bins: ",theo_at_bins) 
    
    #print("scaled y_num[0]:", y_num[0])
    #print("analytic[0]:", theo_at_bins[0])

    diff     = y_num - theo_at_bins
    #print("diff: ",diff)
    rel_diff = 100. * diff / np.where(np.abs(theo_at_bins) > 1e-20, theo_at_bins, np.nan)
    #print("rel_diff: ",rel_diff)
    out = np.column_stack([x_num, y_num, err_num, theo_at_bins, diff, rel_diff])
    np.savetxt(outdir_vals + "Q2.csv", out, delimiter=",",header="Q2_GeV2,num,num_err,ana,diff,rel_diff_percent", comments="")


    fig_diff, ax = plt.subplots()
    ax.errorbar(x_num, rel_diff, yerr=100.*err_num/theo_at_bins,fmt='o', markersize=2)
    ax.axhline(0, color='k', linestyle='--')
    ax.set_xlabel(r"$Q^2\ (\mathrm{GeV}^2)$")
    ax.set_ylabel(r"$(num - ana)/ana\ (\%)$")
    ax.set_title("Relative difference numeric vs. analytic (LO)")
    ax.set_ylim(ylow_diff,yup_diff)
    save_figure(fig_diff, f"{savename}_rel_diff", outdir=outdir)
    plt.close(fig_diff)


    K = mergebins(divideplots(nlo_s, full_s), 5)
    plot_K(ax_k, K, colors["K"], r"$Q^2\ (\mathrm{GeV}^2)$")
    
    ax_k.set_xlim(ax_main.get_xlim())
    mulify(fig_main, delx=4.5, dely=1.) 
    save_figure(fig_main, savename, outdir=outdir)
    plt.close(fig_main)

