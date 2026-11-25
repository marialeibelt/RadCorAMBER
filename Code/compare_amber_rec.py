#!/usr/bin/env python3
"""TGEANT truth vs CORAL/PHAST RECO: distribution comparisons with common bins.

Uses per-event double-precision kinematic branches so the comparison does not
inherit different pre-existing TH1 binnings. This is NOT event matching.
"""
from __future__ import annotations
import argparse
from array import array
import math
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

BEFORE_COLOR = '#717b83'   # TGEANT primary kinematics
AFTER_COLOR = '#2ca02c'    # CORAL / PHAST reconstruction
plt.rcParams.update({
    'font.family': 'serif', 'font.size': 12, 'axes.titlesize': 18,
    'axes.labelsize': 14, 'xtick.labelsize': 12, 'ytick.labelsize': 12,
    'axes.grid': True, 'grid.linestyle': ':', 'grid.alpha': 0.28,
    'savefig.dpi': 300,
})


def validate_tree(f, name, required):
    tree = f.Get(name)
    if not tree or not tree.InheritsFrom('TTree'):
        raise RuntimeError(f'Missing TTree {name}')
    missing = [k for k in required if not tree.GetBranch(k)]
    if missing:
        raise RuntimeError(f'Missing branches in {name}: {missing}')
    return tree


def root_hist_to_arrays(hist):
    n = hist.GetNbinsX()
    return {
        'edges': np.array([hist.GetXaxis().GetBinLowEdge(i) for i in range(1, n+2)], dtype=float),
        'counts': np.array([hist.GetBinContent(i) for i in range(1, n+1)], dtype=float),
        'errors': np.array([hist.GetBinError(i) for i in range(1, n+1)], dtype=float),
        'entries': float(hist.GetEntries()),
        'underflow': float(hist.GetBinContent(0)),
        'overflow': float(hist.GetBinContent(n+1)),
    }


def book_histograms(root, tree, mapping, edges_by_var, flag):
    """Lazy-book all TH1 actions so ROOT can fill several histograms in one scan."""
    rdf = root.RDataFrame(tree)
    filtered = rdf.Filter(f'{flag} == 1')
    jobs = {}
    for name, branch in mapping.items():
        edges = edges_by_var[name]
        filt = filtered.Filter(f'std::isfinite({branch}) && {branch} >= 0.0')
        if name == 'Q2_log':
            filt = filtered.Filter(f'std::isfinite({branch}) && {branch} > 0.0')
        model = root.RDF.TH1DModel(f'compare_{flag}_{name}', '', len(edges)-1,
                                  array('d', edges.tolist()))
        jobs[name] = filt.Histo1D(model, branch)
    try:
        root.RDF.RunGraphs(list(jobs.values()))
    except (AttributeError, TypeError):
        print('INFO: Filling histograms one by one (RunGraphs unavailable).', flush=True)
    return {name: root_hist_to_arrays(job.GetValue()) for name, job in jobs.items()}


def draw_comparison(before, after, output, title, xlabel, mode='raw', logx=False):
    fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
    ymax = 0.
    for h, color, label in [(before, BEFORE_COLOR, 'TGEANT (before)'),
                            (after, AFTER_COLOR, 'CORAL/PHAST (after)')]:
        counts = h['counts'].astype(float).copy()
        if mode == 'density':
            area = float(np.sum(counts))
            widths = np.diff(h['edges'])
            counts = np.divide(counts, widths * area, out=np.zeros_like(counts),
                               where=(widths > 0)) if area else np.zeros_like(counts)
        if np.any(counts > 0):
            ax.stairs(counts, h['edges'], color=color, lw=2, label=label, zorder=3)
            ymax = max(ymax, float(np.max(counts)))
        else:
            ax.plot([], [], color=color, label=label)
    if logx:
        ax.set_xscale('log')
    ax.set_xlim(before['edges'][0] if not logx else max(before['edges'][0], 1e-300),
                before['edges'][-1])
    ax.set_ylim(bottom=0, top=max(ymax*1.12, 1e-12))
    ax.set_ylabel('Normalized density' if mode == 'density' else 'Events / bin')
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=16)
    ax.legend(loc='upper right', frameon=False, fontsize=12)
    ax.set_axisbelow(True)
    for s in ax.spines.values(): s.set_linewidth(.9)
    fig.savefig(output, dpi=300)
    plt.close(fig)


def read_existing_hist(f, name):
    h = f.Get(name)
    if not h or not h.InheritsFrom('TH1') or h.InheritsFrom('TH2'):
        return None
    return root_hist_to_arrays(h)


def photon_diagnostic(before, after, outdir):
    hb = read_existing_hist(before, 'h_photon_energy')
    ha = read_existing_hist(after, 'UserEvent300/05_Photon_CORAL_Reco/E_gamma_reco_candidates')
    if hb is None or ha is None:
        print('WARNING: Missing stored photon/cluster histogram; skipping photon diagnostic.')
        return False
    for log in (False, True):
        fig, ax = plt.subplots(figsize=(10.5, 5.8), constrained_layout=True)
        for h, col, label in [(hb, BEFORE_COLOR, 'TGEANT: real photons'),
                               (ha, AFTER_COLOR, 'PHAST: neutral ECAL2 clusters')]:
            if np.any(h['counts']):
                if log:
                    # Uniform bins starting at zero: use strictly-positive bin centers;
                    # a log axis does NOT recover information below the first bin width.
                    x = 0.5*(h['edges'][1:] + h['edges'][:-1])
                    ax.plot(x[x > 0], h['counts'][x > 0], color=col, lw=1.9, label=label)
                else:
                    ax.stairs(h['counts'], h['edges'], color=col, lw=1.9, label=label)
            else:
                ax.plot([], [], color=col, label=label)
        if log:
            ax.set_xscale('log')
            ax.set_xlim(min(max(hb['edges'][1]/2, 1e-4), max(ha['edges'][1]/2, 1e-4)),
                        min(hb['edges'][-1], ha['edges'][-1]))
        else:
            ax.set_xlim(0, min(hb['edges'][-1], ha['edges'][-1]))
        ax.set_ylim(bottom=0)
        ax.set_xlabel(r'$E\ (\mathrm{GeV})$')
        ax.set_ylabel('Particles / clusters per bin')
        ax.set_title('Photon Emission vs. ECAL2 Cluster Response', pad=16)
        ax.legend(loc='upper right', frameon=False)
        ax.text(.02, .02, 'Different objects; not an efficiency or matched-photon comparison',
                transform=ax.transAxes, fontsize=10, color='0.4', va='bottom')
        fig.savefig(outdir / f'photon_vs_cluster_energy_{"logx" if log else "linear"}.png', dpi=300)
        plt.close(fig)
    return True


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument('--before', required=True)
    p.add_argument('--after', required=True)
    p.add_argument('--outdir', required=True)
    p.add_argument('--before-tree', default='TGEANTTruth')
    p.add_argument('--after-tree', default='UserEvent300/PRM')
    p.add_argument('--q2-log-min', type=float, default=1e-10)
    p.add_argument('--q2-log-max', type=float, default=.05)
    p.add_argument('--q2-log-bins', type=int, default=450)
    p.add_argument('--photon-diagnostic', type=int, choices=(0,1), default=1)
    a = p.parse_args()
    if not (0 < a.q2_log_min < a.q2_log_max and a.q2_log_bins >= 2):
        p.error('Require 0 < Q2_LOG_MIN < Q2_LOG_MAX and Q2_LOG_BINS >= 2')
    for source in (a.before, a.after):
        if not Path(source).is_file():
            p.error(f'Missing ROOT file: {source}')
    try:
        import ROOT
    except ImportError as e:
        raise SystemExit(f'ERROR: PyROOT not available after ROOT_SETUP: {e}')
    ROOT.gROOT.SetBatch(True)
    before = ROOT.TFile.Open(a.before, 'READ')
    after = ROOT.TFile.Open(a.after, 'READ')
    if not before or before.IsZombie() or not after or after.IsZombie():
        raise SystemExit('ERROR: Cannot open one of the ROOT files')
    t_before = validate_tree(before, a.before_tree,
        ['validMuonKinematics','muonEnergy_GeV','thetaMuon_mrad','Q2_GeV2'])
    t_after = validate_tree(after, a.after_tree,
        ['has_muon_reco','E_mu_reco','theta_mu_reco_z','Q2_reco'])
    n_before, n_after = int(t_before.GetEntries()), int(t_after.GetEntries())
    print(f'TGEANT entries: {n_before:,}; PHAST entries: {n_after:,}', flush=True)
    if n_before and abs(n_before-n_after)/n_before > 0.01:
        print('WARNING: Tree entry counts differ by >1%. Are these the same simulation subset?', flush=True)
    out = Path(a.outdir).expanduser()
    out.mkdir(parents=True, exist_ok=True)
    edges = {
        'E_mu': np.linspace(60, 105, 181),
        'theta_mu': np.linspace(0, 5, 251),
        'Q2': np.linspace(0, .05, 251),
        'Q2_log': np.geomspace(a.q2_log_min, a.q2_log_max, a.q2_log_bins+1),
    }
    before_mapping = {'E_mu': 'muonEnergy_GeV', 'theta_mu': 'thetaMuon_mrad',
                      'Q2': 'Q2_GeV2', 'Q2_log':'Q2_GeV2'}
    after_mapping = {'E_mu': 'E_mu_reco', 'theta_mu': 'theta_mu_reco_z',
                     'Q2': 'Q2_reco', 'Q2_log':'Q2_reco'}
    print('Histogramming TGEANT truth branches ...', flush=True)
    hb = book_histograms(ROOT, t_before, before_mapping, edges, 'validMuonKinematics')
    print('Histogramming PHAST reco branches ...', flush=True)
    ha = book_histograms(ROOT, t_after, after_mapping, edges, 'has_muon_reco')
    kinds = {
        'E_mu': ('Scattered Muon Energy', r'$E_\mu\ (\mathrm{GeV})$'),
        'theta_mu': ('Muon Scattering Angle', r'$\theta_\mu^{\mathrm{lab}}\ (\mathrm{mrad})$'),
        'Q2': ('Momentum Transfer', r'$Q^2\ (\mathrm{GeV}^2)$'),
        'Q2_log': ('Momentum Transfer — Fine Log Binning', r'$Q^2\ (\mathrm{GeV}^2)$'),
    }
    for name, (title, label) in kinds.items():
        for mode in ['raw', 'density']:
            path = out / f'{name}_before_after_{mode}.png'
            draw_comparison(hb[name], ha[name], path, title, label, mode=mode, logx=name=='Q2_log')
            print('  ', path, flush=True)
    diag = photon_diagnostic(before, after, out) if a.photon_diagnostic else False
    lines = [
        'TGEANT BEFORE vs PHAST AFTER -- RECONSTRUCTION COMPARISON',
        '='*60,
        f'TGEANT file: {a.before}', f'TGEANT tree: {a.before_tree}',
        f'TGEANT tree entries: {n_before:,}',
        f'PHAST file: {a.after}', f'PHAST tree: {a.after_tree}',
        f'PHAST tree entries: {n_after:,}',
        f'Difference (TGEANT - PHAST): {n_before-n_after:+,}',
        '',
        'Only valid muon kinematics / reconstructed-muon events are histogrammed.',
        'Counts are unnormalized; density plots are normalized separately to unit integral.',
        'Both sample files MUST correspond to the same input event subset.',
        'These plots are distribution-level comparisons, NOT event-by-event matching.',
        'No MC-truth fields from PHAST are used.',
        '',
        'Per-variable statistics (in range / underflow / overflow):',
    ]
    for name in kinds:
        b,h = hb[name], ha[name]
        lines.append(f'  {name}:')
        for label, z in [('TGEANT',b),('PHAST',h)]:
            lines.append(f'    {label}: filled={int(z["entries"]):,} '
                         f'in_range={int(np.sum(z["counts"])):,} '
                         f'underflow={int(z["underflow"]):,} '
                         f'overflow={int(z["overflow"]):,}')
    lines += [
        '',
        'Photon diagnostics made: ' + str(diag),
        'WARNING: Generated photons and neutral ECAL2 clusters are different objects.',
        'Photon diagnostic is NOT an efficiency/response measurement.',
        'WARNING: The TGEANT primary-muon selection and unit conversion should',
        'be independently validated before interpreting physical shifts.',
        'A binning lower edge of Q2=1e-10 GeV^2 is NOT detector resolution.',
        '',
    ]
    report = out/'comparison_statistics.txt'
    report.write_text('\n'.join(lines), encoding='utf-8')
    print('  ',report)
    before.Close(); after.Close()


if __name__=='__main__':
    try:
        main()
    except (OSError,RuntimeError) as e:
        sys.exit(f'ERROR: {e}')
