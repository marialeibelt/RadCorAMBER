#!/usr/bin/env python3
"""Standalone, fast PHAST-reconstruction plotting from the merged ROOT file.

Uses existing *reconstructed* TH1 histograms from UserEvent300.
Uses saved equal-width RECO TH1 histograms except the muon energy,
Q2 and both muon angles. These are re-histogrammed from the reconstructed
TTree using full, data-driven ranges because the original TH1s clip at
60-105 GeV, Q2=0.05 GeV^2 and angle=5 mrad.
This is plotting only; it does not rerun reconstruction.
No MC truth, no new reconstruction, no event reweighting.
Photon energy denotes neutral ECAL2 cluster *candidates*, not verified gammas.

The photon-energy linear/log-x plots use exactly the SAME uniform-width
underlying histogram. Log-x uses positive bin centers because the first
ROOT bin begins at zero; no rebinning and no artificially fine bins.
Distributions (including photon energy, E_mu, angles and Q2) use log-y
by default; the cutflow stays linear. Nonpositive bins are hidden on log-y.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import LogLocator, LogFormatterMathtext
import numpy as np

COLOR = '#2ca02c'
BAND = '#2ca02c'
PLOTTER_VERSION = 'logy_uniform_full_muon_distributions_plus_theta5_q2_zoom'
TREE_PATH = 'UserEvent300/PRM'
CUTFLOW_PATH = 'UserEvent300/00_EventSelection/cutflow'
DIR_RECO = 'UserEvent300/01_CORAL_Reco'
DIR_GAMMA = 'UserEvent300/05_Photon_CORAL_Reco'

# These four plots are rebuilt from event-level reco values so that the
# original ROOT histogram overflow/underflow cannot silently cut the tails.
REBUILT_TREE_HISTS = {
    f'{DIR_RECO}/E_mu_reco',
    f'{DIR_RECO}/Q2_reco',
    f'{DIR_RECO}/theta_mu_reco_z',
    f'{DIR_RECO}/theta_mu_reco_beam',
    f'{DIR_RECO}/theta_mu_reco_z_theta5_zoom',
    f'{DIR_RECO}/Q2_reco_theta5_zoom',
}

plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.titlesize': 18,
    'axes.labelsize': 14,
    'xtick.labelsize': 12,
    'ytick.labelsize': 12,
    'axes.grid': True,
    'grid.alpha': 0.28,
    'grid.linestyle': ':',
    'grid.linewidth': 0.8,
    'savefig.dpi': 300,
    'text.usetex': False,
})

# All bins are equidistant; E_mu, Q2 and muon angles are rebuilt from the reco tree.
# (filename, ROOT histogram path, title, x axis, y axis, log x)
PLOTS = [
    ('photon_candidate_energy_linear', f'{DIR_GAMMA}/E_gamma_reco_candidates',
     'Neutral ECAL2 Cluster Energy', r'$E_{\mathrm{cluster}}\ (\mathrm{GeV})$', 'Clusters / bin', False),
    ('photon_candidate_energy_logx', f'{DIR_GAMMA}/E_gamma_reco_candidates',
     'Neutral ECAL2 Cluster Energy', r'$E_{\mathrm{cluster}}\ (\mathrm{GeV})$', 'Clusters / bin', True),
    ('photon_candidate_angle', f'{DIR_GAMMA}/theta_gamma_reco_candidates',
     'Neutral ECAL2 Cluster Polar Angle', r'$\theta_{\mathrm{cluster}}\ (\mathrm{mrad})$', 'Clusters / bin', False),
    ('photon_candidate_multiplicity', f'{DIR_GAMMA}/n_gamma_reco_candidates',
     'Neutral ECAL2 Cluster Multiplicity', r'$N_{\mathrm{neutral\ clusters}}$', 'Events / bin', False),
    ('photon_candidate_multiplicity_Egt200MeV', f'{DIR_GAMMA}/n_gamma_reco_candidates_Egt200MeV',
     'Neutral ECAL2 Clusters Above 200 MeV', r'$N_{\mathrm{clusters},\ E>200\,\mathrm{MeV}}$', 'Events / bin', False),
    ('muon_energy', f'{DIR_RECO}/E_mu_reco',
     'Scattered Muon Energy', r'$E_\mu^{\mathrm{reco}}\ (\mathrm{GeV})$', 'Events / bin', False),
    ('muon_angle_lab', f'{DIR_RECO}/theta_mu_reco_z',
     'Muon Scattering Angle', r'$\theta_\mu^{\mathrm{lab,reco}}\ (\mathrm{mrad})$', 'Events / bin', False),
    ('muon_angle_to_beam', f'{DIR_RECO}/theta_mu_reco_beam',
     'Muon Angle Relative to Beam', r'$\theta_{\mu,\mathrm{beam}}^{\mathrm{reco}}\ (\mathrm{mrad})$', 'Events / bin', False),
    ('momentum_transfer', f'{DIR_RECO}/Q2_reco',
     'Momentum Transfer', r'$Q_{\mathrm{reco}}^2\ (\mathrm{GeV}^2)$', 'Events / bin', False),
    ('muon_angle_lab_zoom', f'{DIR_RECO}/theta_mu_reco_z_theta5_zoom',
     'Muon Scattering Angle (Zoom: 0–5 mrad)', r'$\theta_\mu^{\mathrm{lab,reco}}\ (\mathrm{mrad})$', 'Events / bin', False),
    ('momentum_transfer_zoom', f'{DIR_RECO}/Q2_reco_theta5_zoom',
     'Momentum Transfer (Zoom, theta < 5 mrad)', r'$Q_{\mathrm{reco}}^2\ (\mathrm{GeV}^2)$', 'Events / bin', False),
    ('reconstructed_vertex_z', f'{DIR_RECO}/vertex_z_reco',
     'Reconstructed Primary Vertex', r'$z_{\mathrm{vtx}}^{\mathrm{reco}}\ (\mathrm{cm})$', 'Events / bin', False),
    ('muon_track_chi2ndf', f'{DIR_RECO}/track_chi2ndf',
     'Reconstructed Muon Track Quality', r'$\chi^2/\mathrm{ndf}$', 'Events / bin', False),
    ('primary_vertex_chi2', f'{DIR_RECO}/vertex_chi2',
     'Primary Vertex Fit Quality', r'$\chi^2_{\mathrm{vtx}}$', 'Events / bin', False),
]


def read_hist(root_file, path):
    return read_hist_object(root_file.Get(path), path)


def read_hist_object(obj, path):
    if not obj or not obj.InheritsFrom('TH1') or obj.InheritsFrom('TH2'):
        return None
    n = int(obj.GetNbinsX())
    edges = np.asarray([float(obj.GetXaxis().GetBinLowEdge(i)) for i in range(1, n+2)], dtype=float)
    widths = np.diff(edges)
    # All physical plots must retain equidistant bins. Never silently use
    # variable-width bins or transform the stored distribution.
    if not np.all(np.isfinite(widths)) or not np.allclose(widths, widths[0], rtol=1e-8, atol=1e-13):
        raise RuntimeError(f'ROOT histogram {path} does not have equidistant bins.')
    counts = np.asarray([float(obj.GetBinContent(i)) for i in range(1, n+1)], dtype=float)
    errors = np.asarray([float(obj.GetBinError(i)) for i in range(1, n+1)], dtype=float)
    return {
        'edges': edges, 'counts': counts, 'errors': errors,
        'underflow': float(obj.GetBinContent(0)),
        'overflow': float(obj.GetBinContent(n+1)),
        'entries': float(obj.GetEntries()),
        'obj': obj,
    }



def make_full_muon_histograms(tree):
    """Re-histogram *all* reconstructed muon energies, Q2 and angles.

    Every x bin has the same width in its own histogram. Ranges extend
    beyond the observed maxima, never truncating at the old fixed TH1
    limits. Zero is the physical lower edge. No event-selection cuts
    beyond has_muon_reco and finite/nonnegative values are imposed.
    """
    import math
    import ROOT

    variables = (
        ('E_mu_reco', 180, 105.0, 5.0, 'GeV'),
        ('Q2_reco', 250, 0.18, 0.01, 'GeV^2'),
        ('theta_mu_reco_z_mrad', 250, 5.0, 0.5, 'mrad'),
        ('theta_mu_reco_beam_mrad', 250, 5.0, 0.5, 'mrad'),
    )
    needed = {'has_muon_reco'} | {item[0] for item in variables}
    branches = {str(b.GetName()) for b in tree.GetListOfBranches()}
    absent = sorted(needed - branches)
    if absent:
        raise RuntimeError('Missing reconstructed tree branches: ' + ', '.join(absent))

    df = ROOT.RDataFrame(tree).Filter('has_muon_reco != 0')
    dfs = {}
    extrema = {}
    for branch, _, _, _, _ in variables:
        valid = df.Filter(f'std::isfinite({branch}) && {branch} >= 0.0')
        dfs[branch] = valid
        extrema[branch] = (valid.Count(), valid.Max(branch))

    # Run all maxima/counts in one pass over the merged reco TTree.
    ROOT.RDF.RunGraphs([pointer for pair in extrema.values() for pointer in pair])

    hist_pointers = {}
    for branch, nbins, default_upper, rounding, unit in variables:
        count_ptr, max_ptr = extrema[branch]
        nvalid = int(count_ptr.GetValue())
        xmax = float(max_ptr.GetValue()) if nvalid else 0.0
        if nvalid and (not math.isfinite(xmax) or xmax < 0):
            raise RuntimeError(f'Invalid maximum for {branch}: {xmax}')
        # +2% padding and round UP to a convenient upper limit.
        upper = max(default_upper, math.ceil(xmax * 1.02 / rounding) * rounding)
        # Guard against exact equality to an upper bin edge, without using
        # the current sample's maximum as the histogram's actual last edge.
        if nvalid and xmax >= upper:
            upper += rounding
        print(f'FULL RANGE: {branch} = [0, {upper:.9g}] {unit}; '
              f'valid entries={nvalid:,}; observed max={xmax:.9g} {unit}')
        name = f'{branch}_full_reco_plot'
        hist_pointers[branch] = (dfs[branch].Histo1D((name, name, nbins, 0.0, upper), branch), nvalid)

    # Additional plots only: same angular selection for both zoom histograms.
    # No change to the four full-range histograms above.
    zoom = df.Filter(
        'std::isfinite(theta_mu_reco_z_mrad) && '
        'theta_mu_reco_z_mrad >= 0.0 && theta_mu_reco_z_mrad < 5.0 && '
        'std::isfinite(Q2_reco) && Q2_reco >= 0.0'
    )
    n_zoom = zoom.Count()
    angle_zoom = zoom.Histo1D(
        ('theta_mu_reco_z_theta5_zoom', 'Reco theta 0-5 mrad', 250, 0.0, 5.0),
        'theta_mu_reco_z_mrad'
    )
    q2_zoom = zoom.Histo1D(
        ('Q2_reco_theta5_zoom', 'Reco Q2 0-0.18 GeV2, theta <5 mrad', 250, 0.0, 0.18),
        'Q2_reco'
    )

    # Second pass: fill all full-range and zoom histograms together.
    ROOT.RDF.RunGraphs([pair[0] for pair in hist_pointers.values()]
                       + [n_zoom, angle_zoom, q2_zoom])

    keys = {
        'E_mu_reco': f'{DIR_RECO}/E_mu_reco',
        'Q2_reco': f'{DIR_RECO}/Q2_reco',
        'theta_mu_reco_z_mrad': f'{DIR_RECO}/theta_mu_reco_z',
        'theta_mu_reco_beam_mrad': f'{DIR_RECO}/theta_mu_reco_beam',
    }
    output = {}
    for branch, (hist_ptr, nvalid) in hist_pointers.items():
        hist_obj = hist_ptr.GetValue().Clone(branch + '_full_plot_detached')
        hist_obj.SetDirectory(0)
        h = read_hist_object(hist_obj, keys[branch])
        if h is None:
            raise RuntimeError('Cannot rebuild histogram for ' + branch)
        included = float(np.sum(h['counts']))
        # A floating-point ROOT histogram can contain weighted fills; the
        # reconstructed muon histograms here are unweighted (1 per event).
        if abs(included - nvalid) > 0.5 or h['underflow'] or h['overflow']:
            raise RuntimeError(f'Histogram dropped entries for {branch}: '
                               f'expected={nvalid}, included={included:g}, '
                               f'underflow={h["underflow"]:g}, '
                               f'overflow={h["overflow"]:g}')
        output[keys[branch]] = h
    # The zoomed angle histogram contains the full common angular selection.
    # The Q2 plot shows only its 0-0.18 GeV^2 portion; Q2 outliers remain
    # available in the original momentum_transfer.png full-range plot.
    n_selected = int(n_zoom.GetValue())
    for tag, ptr, path in (
        ('angle', angle_zoom, f'{DIR_RECO}/theta_mu_reco_z_theta5_zoom'),
        ('Q2', q2_zoom, f'{DIR_RECO}/Q2_reco_theta5_zoom'),
    ):
        obj = ptr.GetValue().Clone(tag + '_theta5_zoom_detached')
        obj.SetDirectory(0)
        h = read_hist_object(obj, path)
        if h is None:
            raise RuntimeError('Could not make zoom histogram ' + path)
        recorded = float(np.sum(h['counts']) + h['underflow'] + h['overflow'])
        if abs(recorded - n_selected) > 0.5:
            raise RuntimeError(f'Zoom {tag}: selected={n_selected}, histogram={recorded:g}')
        output[path] = h
    print(f'ADDITIONAL ZOOMS: theta_mu < 5 mrad; selected={n_selected:,}; '
          f'Q2 in [0, 0.18) GeV^2={int(round(np.sum(output[f"{DIR_RECO}/Q2_reco_theta5_zoom"]["counts"]))):,}')
    return output


def make_plot(hist, output, title, x_label, y_label, logx=False, logy=True):
    """Draw saved equal-width ROOT histogram, with explicit y-axis scaling.

    In log-y mode only strictly positive contents are displayed. No zeros
    are replaced by artificial event counts. A histogram without positive
    entries still receives a logarithmic y-axis (but no plotted distribution).
    """
    edges = hist['edges']
    counts = np.asarray(hist['counts'], dtype=float)
    errors = np.asarray(hist['errors'], dtype=float)
    centers = (edges[1:] + edges[:-1]) / 2.0
    positive = counts[np.isfinite(counts) & (counts > 0)]
    has_positive = bool(positive.size)

    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    if logy:
        # Set log scale REGARDLESS of whether the histogram is empty.
        ax.set_yscale('log', nonpositive='mask')
        ax.yaxis.set_major_locator(LogLocator(base=10.0, numticks=12))
        ax.yaxis.set_major_formatter(LogFormatterMathtext(base=10.0))
        ax.yaxis.set_minor_locator(LogLocator(base=10.0, subs=np.arange(2, 10), numticks=50))
        if has_positive:
            floor = max(float(np.min(positive)) * .45, np.finfo(float).tiny)
            # Include uncertainty bands without inventing positive entries.
            upper_values = counts[np.isfinite(counts) & (counts > 0)] + errors[np.isfinite(counts) & (counts > 0)]
            ymax = float(np.max(upper_values))
            ceiling = max(ymax * 1.3, floor * 10.0)
            ax.set_ylim(floor, ceiling)
        else:
            floor = 0.1
            ax.set_ylim(0.1, 10.0)
        plotted = np.where(counts > 0, counts, np.nan)
        band_lower = np.where(counts > 0, np.maximum(counts - errors, floor), np.nan)
        band_upper = np.where(counts > 0, counts + errors, np.nan)
    else:
        plotted = counts
        band_lower = np.maximum(counts - errors, 0.0)
        band_upper = counts + errors
        ymax = float(np.max(np.maximum(counts+errors, counts))) if counts.size else 0.0
        ax.set_ylim(0, max(1.1, ymax*1.09))

    # Drawing method changes with x-axis only. Both photon-energy PNGs use
    # the *same* underlying equal-width ROOT bins.
    if logx:
        eligible = centers > 0
        x = centers[eligible]
        if len(x) and has_positive:
            ax.plot(x, plotted[eligible], color=COLOR, lw=2.1,
                    label='NLO', zorder=3)
            ax.fill_between(x, band_lower[eligible], band_upper[eligible],
                            color=BAND, alpha=0.17, zorder=2)
        if len(x):
            ax.set_xlim(x[0]*0.8, max(edges[-1], x[-1]))
        ax.set_xscale('log')
    else:
        if has_positive:
            ax.stairs(plotted, edges, color=COLOR, linewidth=2.0,
                      label='NLO', zorder=3)
            ax.fill_between(centers, band_lower, band_upper,
                            color=BAND, alpha=0.17, step='mid', zorder=2)
        ax.set_xlim(float(edges[0]), float(edges[-1]))

    if has_positive:
        ax.legend(loc='upper right', frameon=False, fontsize=13)
    else:
        message = ('No positive reconstructed entries' if logy
                   else 'No reconstructed candidates / no entries')
        ax.text(0.5, 0.5, message, horizontalalignment='center',
                verticalalignment='center', transform=ax.transAxes,
                fontsize=13)

    ax.set_title(title, pad=16)
    ax.set_xlabel(x_label)
    ax.set_ylabel(y_label)
    ax.set_axisbelow(True)
    for spine in ax.spines.values():
        spine.set_linewidth(0.9)
    # Fail fast if a future code change accidentally resets the y scale.
    expected_scale = 'log' if logy else 'linear'
    if ax.get_yscale() != expected_scale:
        raise RuntimeError(f'Unexpected y-axis scale for {output}: '
                           f'expected {expected_scale}, got {ax.get_yscale()}')
    fig.savefig(output, dpi=300,
                metadata={'Description': f'CORAL/PHAST {PLOTTER_VERSION}; '
                                         f'xscale={ax.get_xscale()}; '
                                         f'yscale={ax.get_yscale()}; '
                                         'uniform ROOT histogram bin widths'})
    plt.close(fig)


def make_cutflow(hist, output, logy=False):
    obj = hist['obj']
    labels = [str(obj.GetXaxis().GetBinLabel(i)) or f'Stage {i}'
              for i in range(1, obj.GetNbinsX()+1)]
    counts = hist['counts']
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    x = np.arange(len(counts))
    ax.bar(x, counts, width=0.62, color=COLOR, alpha=0.80)
    ax.set_xticks(x, labels, rotation=15, ha='right')
    ax.set_ylabel('Events')
    ax.set_title('CORAL Reconstruction and Event Selection', pad=16)
    if logy and np.any(counts > 0):
        pos = counts[counts > 0]
        ax.set_yscale('log', nonpositive='mask')
        ax.set_ylim(float(np.min(pos))*0.5, max(float(np.max(counts))*1.25, 1))
    else:
        ax.set_ylim(0, max(float(np.max(counts))*1.12, 1.0))
    ax.set_axisbelow(True)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def whole_or_float(x):
    return f'{round(x):,}' if abs(x-round(x)) < 1e-6 else f'{x:,.6g}'


def write_statistics(path, input_file, total_entries, available, missing):
    lines = [
        'PHAST / CORAL RECONSTRUCTION EVENT STATISTICS',
        '='*58,
        f'Input ROOT: {input_file}',
        f'Plotter version: {PLOTTER_VERSION}',
        f'Tree: {TREE_PATH}',
        f'Tree entries (actual processed events): {total_entries:,}',
        '',
        'Selection / reconstruction cutflow (independent of MC truth):',
    ]
    flow = available.get(CUTFLOW_PATH)
    if flow is not None:
        obj = flow['obj']
        for i, val in enumerate(flow['counts'], start=1):
            label = str(obj.GetXaxis().GetBinLabel(i)) or f'Stage {i}'
            ratio = f' ({100*val/total_entries:.3f}% of tree entries)' if total_entries else ''
            lines.append(f'  {label}: {whole_or_float(val)}{ratio}')
        if flow['counts'].size:
            diff = float(flow['counts'][0]) - float(total_entries)
            if abs(diff) > 0:
                lines.append(f'  WARNING: generated count minus tree entries = {diff:,.0f}')
    else:
        lines.append('  Cutflow histogram not available.')
    lines += [
        '',
        'Reconstructed photon candidates = neutral ECAL2 clusters (NOT confirmed photons).',
        'Only reconstructed quantities are summarized below.',
    ]
    for key in [f'{DIR_GAMMA}/n_gamma_reco_candidates',
                f'{DIR_GAMMA}/n_gamma_reco_candidates_Egt200MeV']:
        h = available.get(key)
        if h is None:
            lines.append(f'  MISSING: {key}')
            continue
        n_zero = float(h['counts'][0]) if h['counts'].size else 0.0
        n_any = float(np.sum(h['counts'][1:]) + h['overflow'])
        n_clusters = float(np.dot(np.arange(len(h['counts'])), h['counts']))
        lines.append(f'  {key.split("/")[-1]}:')
        lines.append(f'    Events with zero candidates: {whole_or_float(n_zero)}')
        lines.append(f'    Events with >=1 candidate (includes overflow bin): {whole_or_float(n_any)}')
        lines.append(f'    Candidates counted in displayed bins: {whole_or_float(n_clusters)}')
    ph = available.get(f'{DIR_GAMMA}/E_gamma_reco_candidates')
    if ph is not None:
        lines.append(f'  Photon-candidate energy histogram: {whole_or_float(ph["entries"])} cluster fills')
        lines.append(f'  Photon-candidate energy bin width: {(ph["edges"][1] - ph["edges"][0]):.6g} GeV')
    lines += [
        '',
        'Histogram statistics (all bins equidistant; E_mu, Q2 and angles rebuilt from reco tree):',
    ]
    for k in sorted(available):
        if k == CUTFLOW_PATH:
            continue
        h = available[k]
        width = float(h['edges'][1] - h['edges'][0])
        lines.append(f'  {k}: bins={len(h["counts"])} width={width:.9g} [x-axis units] '
                     f'entries={whole_or_float(h["entries"])} '
                     f'in_range={whole_or_float(float(sum(h["counts"])))} '
                     f'underflow={whole_or_float(h["underflow"])} '
                     f'overflow={whole_or_float(h["overflow"])}')
    if missing:
        lines += ['', 'MISSING HISTOGRAMS (no plot generated):']
        lines += [f'  {name}' for name in sorted(missing)]
    lines += [
        '',
        'Important: E_gamma_reco in the tree is populated only after MC matching.',
        'Therefore photon-candidate energy plots use the saved ECAL2 histogram.',
        'Linear-x and log-x photon plots use EXACTLY the same equidistant binning.',
        'Distributions use logarithmic y by default; zero/nonpositive bins are not plotted on log-y.',
        'Log-x uses positive bin centers; the original first bin spans 0-0.5 GeV.',
        'Additional zooms (both based on theta_mu_reco_z < 5 mrad and finite Q2 >= 0):',
        '  muon_angle_lab_zoom.png: theta_mu_reco_z, 0-5 mrad, 250 uniform bins.',
        '  momentum_transfer_zoom.png: Q2_reco, 0-0.18 GeV^2, 250 uniform bins.',
        '  Full-range muon_angle_lab.png and momentum_transfer.png are unchanged.',
        'No sub-bin information near 0 GeV can be recovered from the saved histogram.',
        '',
    ]
    path.write_text('\n'.join(lines), encoding='utf-8')


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--input', required=True, help='Merged PHAST ROOT file')
    ap.add_argument('--outdir', required=True, help='Folder for PNGs and event_statistics.txt')
    ap.add_argument('--distribution-yscale', choices=('log', 'linear'), default='log',
                    help='y-axis for energy, angle, Q2, vertex, track and photon-candidate plots')
    ap.add_argument('--cutflow-yscale', choices=('log', 'linear'), default='linear',
                    help='y-axis for cutflow plot; default linear')
    args = ap.parse_args()
    src = Path(args.input).expanduser()
    outdir = Path(args.outdir).expanduser()
    if not src.is_file():
        ap.error(f'Input ROOT file does not exist: {src}')
    outdir.mkdir(parents=True, exist_ok=True)
    try:
        import ROOT
    except ImportError as exc:
        raise SystemExit('ERROR: PyROOT unavailable. Source the configured LCG ROOT setup first.\n' + str(exc))
    ROOT.gROOT.SetBatch(True)
    root_file = ROOT.TFile.Open(str(src), 'READ')
    if not root_file or root_file.IsZombie():
        raise SystemExit('ERROR: Cannot read merged PHAST ROOT file.')
    tree = root_file.Get(TREE_PATH)
    if not tree or not tree.InheritsFrom('TTree'):
        raise SystemExit(f'ERROR: Missing tree: {TREE_PATH}')
    entries = int(tree.GetEntries())  # metadata; no 7M-event scan
    print(f'PLOTTER_VERSION={PLOTTER_VERSION} (distribution y scale: {args.distribution_yscale})')
    print(f'Merged PHAST tree: {entries:,} entries')
    keys = {v[1] for v in PLOTS} | {CUTFLOW_PATH}
    available = {}
    missing = []
    for key in sorted(keys - REBUILT_TREE_HISTS):
        h = read_hist(root_file, key)
        if h is None:
            print(f'WARNING: Missing histogram: {key}')
            missing.append(key)
        else:
            available[key] = h
    available.update(make_full_muon_histograms(tree))
    made = 0
    for name, hist_key, title, x_label, y_label, logx in PLOTS:
        h = available.get(hist_key)
        if h is None:
            continue
        target = outdir / f'{name}.png'
        make_plot(h, target, title, x_label, y_label, logx,
                  logy=args.distribution_yscale == 'log')
        print(f'  {target}')
        made += 1
    if CUTFLOW_PATH in available:
        target = outdir / 'reconstruction_cutflow.png'
        make_cutflow(available[CUTFLOW_PATH], target,
                     logy=args.cutflow_yscale == 'log')
        print(f'  {target}')
        made += 1
    # E_mu, Q2 and muon angles were re-histogrammed from the saved TTree
    # to recover tails missing from the old fixed-bin histograms.
    report = outdir / 'event_statistics.txt'
    write_statistics(report, src, entries, available, missing)
    print(f'  {report}')

    root_file.Close()
    if made == 0:
        raise SystemExit('ERROR: ROOT file contains no expected reconstruction histograms.')
    print(f'Distribution y scale: {args.distribution_yscale}; cutflow y scale: {args.cutflow_yscale}')
    print(f'Finished: {made} PNGs + event_statistics.txt (reco TTree read for E_mu/Q2/angles).')


if __name__ == '__main__':
    main()
