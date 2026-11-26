#!/usr/bin/env python3
"""Plot McMule NLO ROOT batches without merging them into one huge ROOT file.

Reads TTree 'Output' using ROOT.RDataFrame in compiled C++ (through PyROOT).
Main plots show actual event counts per equally wide bin, with no
post-processing/rebinning or normalization. For unit-weight ROOT events
these are actual event numbers. Empty bins are drawn explicitly at zero.

Required branches from lhe_to_root.cpp:
  weight, beamEnergy, beamMomentumX/Y/Z,
  scatteredPID (vector<int>), scatteredEnergy/MomentumX/Y/Z (vector<double>).

Usage: python3 analyze_mcmule_evtgen.py --help
"""

from __future__ import annotations

import argparse
from array import array
from datetime import datetime
import math
import os
from pathlib import Path
import sys

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import numpy as np

# Same clean aesthetic as McMule's plotting.py, for ordinary equal-width bins.
plt.rcParams.update({
    'font.family': 'serif',
    'font.size': 12,
    'axes.grid': True,
    'grid.alpha': 0.30,
    'grid.linestyle': ':',
    'savefig.dpi': 300,
    'text.usetex': False,
})

COLORS = {'NLO': '#2ca02c', 'LO': '#1f77b4'}
MUON_MASS_GEV = 0.105658375
PROTON_MASS_GEV = 0.938272088
REQUIRED = (
    'weight', 'beamEnergy', 'beamMomentumX', 'beamMomentumY', 'beamMomentumZ',
    'scatteredPID', 'scatteredEnergy', 'scatteredMomentumX', 'scatteredMomentumY', 'scatteredMomentumZ',
)


def observables():
    """Raw-count histograms; uniform bins on linear X axes, log-spaced bins on log X axes.

    The plotting ranges are NOT analysis cuts. The text report includes aggregate event totals and outside-range counts.
    """
    return [
        ('theta_mu', 'Muon Scattering Angle',
         r'$\theta_\mu^{\mathrm{lab}}\ (\mathrm{mrad})$',
         '', np.linspace(0, 4.5, 121), 'linear', 'auto'),
        ('muon_energy', 'Scattered Muon Energy',
         r'$E_\mu^{\mathrm{lab}}\ (\mathrm{GeV})$',
         '', np.linspace(0, 101, 405), 'linear', 'auto'),   # linear X axis (energies cluster near 100 GeV)
        ('muon_energy_endpoint', 'Scattered Muon Energy (endpoint)',
         r'$E_\mu^{\mathrm{lab}}\ (\mathrm{GeV})$',
         '', np.linspace(99.0, 100.10, 111), 'linear', 'auto'),
        ('Q2', 'Momentum Transfer',
         r'$Q^2\ (\mathrm{GeV}^2)$',
         '', np.linspace(0, 0.18, 151), 'linear', 'auto'),
        ('muon_energy_loss', 'Muon Energy Loss',
         r'$E_{\mathrm{beam}}-E_\mu^{\mathrm{lab}}\ (\mathrm{GeV})$',
         '', np.linspace(0, 100, 401), 'linear', 'auto'),
        ('phi_mu', 'Muon Azimuth',
         r'$\phi_\mu\ (\mathrm{rad})$',
         '', np.linspace(-math.pi, math.pi, 73), 'linear', 'linear'),
        ('photon_energy', 'Real Photon Energy',
         r'$E_\gamma^{\mathrm{lab}}\ (\mathrm{GeV})$',
         '', np.linspace(0, 100, 401), 'linear', 'auto'),
        # Additional plot: logarithmic X axis requires strictly positive X edges.
        # Log-spaced bins, so every bin has the same width on the log axis.
        ('photon_energy_logx', 'Real Photon Energy (log X)',
         r'$E_\gamma^{\mathrm{lab}}\ (\mathrm{GeV})$',
         '', np.geomspace(0.001, 100, 401), 'log', 'auto'),  # 1 MeV - 100 GeV
        ('photon_theta', 'Real Photon Polar Angle',
         r'$\theta_\gamma^{\mathrm{lab}}\ (\mathrm{mrad})$',
         '', np.linspace(0, 1000*math.pi, 361), 'linear', 'auto'),
        ('photon_theta_forward', 'Real Photon Polar Angle (forward region)',
         r'$\theta_\gamma^{\mathrm{lab}}\ (\mathrm{mrad})$',
         '', np.linspace(0, 12.5, 126), 'linear', 'auto'),
    ]


def source_column(plot_key):
    return {
        'theta_mu': '__theta_lab_mrad',
        'muon_energy': '__E_lab_GeV',
        'muon_energy_endpoint': '__E_lab_GeV',
        'Q2': '__Q2_GeV2',
        'muon_energy_loss': '__Eloss_lab_GeV',
        'phi_mu': '__phi_mu_rad',
        'photon_energy': '__ph_E_lab_GeV',
        'photon_energy_logx': '__ph_E_lab_GeV',
        'photon_theta': '__ph_theta_lab_mrad',
        'photon_theta_forward': '__ph_theta_lab_mrad',
    }[plot_key]


def check_and_chain(ROOT, filenames, tag):
    chain = ROOT.TChain('Output')
    per_file = []
    for filename in filenames:
        filepath = Path(filename)
        if not filepath.is_file():
            raise FileNotFoundError(f'{tag} ROOT file does not exist: {filepath}')
        f = ROOT.TFile.Open(str(filepath), 'READ')
        if not f or f.IsZombie():
            raise RuntimeError(f'Cannot read {tag} ROOT file: {filepath}')
        tree = f.Get('Output')
        if not tree:
            f.Close()
            raise RuntimeError(f"No TTree 'Output' in {filepath}")
        branches = {str(br.GetName()) for br in tree.GetListOfBranches()}
        missing = [s for s in REQUIRED if s not in branches]
        if missing:
            f.Close()
            raise RuntimeError(f'Wrong branch names in {filepath}: missing {missing}. '
                               f'Available: {sorted(branches)}')
        n = int(tree.GetEntries())
        per_file.append((str(filepath), n))
        f.Close()
        if n <= 0:
            raise RuntimeError(f'Empty ROOT tree: {filepath}')
        added = int(chain.Add(str(filepath)))
        if added <= 0:
            raise RuntimeError(f'Unable to add ROOT file to TChain: {filepath}')
    if not per_file:
        raise RuntimeError(f'No ROOT files selected for {tag}')
    if chain.GetEntries() != sum(n for _, n in per_file):
        raise RuntimeError(f'TChain entry-count mismatch for {tag}')
    chain.GetEntry(0)
    incoming_E = float(chain.beamEnergy)
    incoming_pz = float(chain.beamMomentumZ)
    return chain, per_file, incoming_E, incoming_pz


def choose_frame(frame, incoming_energy, beam_energy_lab):
    if frame != 'auto':
        return frame
    # Old McMule LHE ROOT converter usually stores COM (beam E ~6.83 GeV).
    # If the incoming beam energy is close to 100 GeV, it is already LAB.
    if incoming_energy < 0.5 * beam_energy_lab:
        return 'cms'
    if incoming_energy > 0.75 * beam_energy_lab:
        return 'lab'
    raise RuntimeError(
        f'Ambiguous input frame: incoming beamEnergy={incoming_energy:.6g} GeV. '
        'Pass --input-frame cms or --input-frame lab explicitly.'
    )


def build_rdf(ROOT, chain, frame, beam_energy_lab):
    # Direction convention: use incoming muon beam as the positive longitudinal
    # direction. This works regardless of whether McMule stores beam pz as +/-.
    if beam_energy_lab <= MUON_MASS_GEV:
        raise ValueError('beam energy must exceed muon rest mass')
    p_lab = math.sqrt(beam_energy_lab**2 - MUON_MASS_GEV**2)
    beta = p_lab / (beam_energy_lab + PROTON_MASS_GEV)
    gamma = 1.0 / math.sqrt(1.0 - beta*beta)
    sgn = '(beamMomentumZ >= 0.0 ? 1.0 : -1.0)'
    # The beam branches are scalars. The outgoing branches, however, are
    # std::vector / ROOT::RVec and can contain a muon, proton and photon(s).
    # Select the actual final-state muon using its PDG code (abs(PID)==13).
    # Never assume scatteredMomentumX[0] is the muon!
    ppar = f'(__mu_Pz_GeV * {sgn})'
    if frame == 'cms':
        elab = f'({gamma:.17g} * (__mu_E_GeV + {beta:.17g} * {ppar}))'
        plablong = f'({gamma:.17g} * ({ppar} + {beta:.17g} * __mu_E_GeV))'
    elif frame == 'lab':
        elab = '__mu_E_GeV'
        plablong = ppar
    else:
        raise ValueError(frame)

    rdf = ROOT.RDataFrame(chain)
    for component, branch in (
        ('E', 'scatteredEnergy'),
        ('Px', 'scatteredMomentumX'),
        ('Py', 'scatteredMomentumY'),
        ('Pz', 'scatteredMomentumZ'),
    ):
        rdf = rdf.Define(f'__mu_{component}_GeV',
                         f'mcmule_analysis::muon_value(scatteredPID, {branch})')

    rdf = rdf.Define('__E_lab_GeV', elab)
    rdf = rdf.Define('__p_lab_parallel_GeV', plablong)
    rdf = rdf.Define('__theta_lab_mrad',
        '1000.0 * std::atan2(std::hypot(__mu_Px_GeV, __mu_Py_GeV), '
        '__p_lab_parallel_GeV)')
    rdf = rdf.Define('__phi_mu_rad',
        'std::atan2(__mu_Py_GeV, __mu_Px_GeV)')
    rdf = rdf.Define('__Q2_GeV2',
        'std::max(0.0, '
        '(beamMomentumX-__mu_Px_GeV)*(beamMomentumX-__mu_Px_GeV) + '
        '(beamMomentumY-__mu_Py_GeV)*(beamMomentumY-__mu_Py_GeV) + '
        '(beamMomentumZ-__mu_Pz_GeV)*(beamMomentumZ-__mu_Pz_GeV) - '
        '(beamEnergy-__mu_E_GeV)*(beamEnergy-__mu_E_GeV))')
    rdf = rdf.Define('__Eloss_lab_GeV', f'({beam_energy_lab:.17g} - __E_lab_GeV)')
    rdf = rdf.Define('__weight_abs', 'std::abs(weight)')
    rdf = rdf.Define('__weight_sq', 'weight * weight')
    rdf = rdf.Define('__n_photons', 'mcmule_analysis::count_photons(scatteredPID)')
    return rdf


def build_photon_rdf(rdf, frame, beam_energy_lab):
    """Select one-photon real-emission events; skip NLO0 events without photons."""
    p_lab = math.sqrt(beam_energy_lab**2 - MUON_MASS_GEV**2)
    beta = p_lab / (beam_energy_lab + PROTON_MASS_GEV)
    gamma = 1.0 / math.sqrt(1.0-beta*beta)
    sgn = '(beamMomentumZ >= 0.0 ? 1.0 : -1.0)'
    rdf = rdf.Filter('__n_photons == 1')
    for comp, branch in (
        ('E', 'scatteredEnergy'),
        ('Px', 'scatteredMomentumX'),
        ('Py', 'scatteredMomentumY'),
        ('Pz', 'scatteredMomentumZ'),
    ):
        rdf = rdf.Define(f'__ph_{comp}_GeV',
                         f'mcmule_analysis::photon_value(scatteredPID, {branch})')
    ppar = f'(__ph_Pz_GeV * {sgn})'
    if frame == 'cms':
        elab = f'({gamma:.17g} * (__ph_E_GeV + {beta:.17g} * {ppar}))'
        plablong = f'({gamma:.17g} * ({ppar} + {beta:.17g} * __ph_E_GeV))'
    elif frame == 'lab':
        elab = '__ph_E_GeV'
        plablong = ppar
    else:
        raise ValueError(frame)
    rdf = rdf.Define('__ph_E_lab_GeV', elab)
    rdf = rdf.Define('__ph_p_lab_parallel_GeV', plablong)
    rdf = rdf.Define('__ph_theta_lab_mrad',
        '1000.0 * std::atan2(std::hypot(__ph_Px_GeV, __ph_Py_GeV), '
        '__ph_p_lab_parallel_GeV)')
    return rdf


def extract_hist(ROOT, rdf, hist_defs, frame, beam_energy_lab, photon_output=None):
    """One pass through the sample for total stats and histograms.

    Photon histograms contain raw counts of events with one real photon;
    NLO0 events without photons do not contribute.
    """
    one_gamma = rdf.Filter('__n_photons == 1')
    photon_rdf = build_photon_rdf(rdf, frame, beam_energy_lab)
    actions = {
        'entries': rdf.Count(),
        'negative_entries': rdf.Filter('weight < 0').Count(),
        'weight_sum': rdf.Sum('weight'),
        'weight_abs_sum': rdf.Sum('__weight_abs'),
        'weight_sq_sum': rdf.Sum('__weight_sq'),
        'one_photon_events': one_gamma.Count(),
        'one_photon_weight_sum': one_gamma.Sum('weight'),
        'multiple_photon_events': rdf.Filter('__n_photons > 1').Count(),
        'zero_photon_events': rdf.Filter('__n_photons == 0').Count(),
    }
    hist_actions = {}
    hist_entry_actions = {}
    for key, _, _, _, edges, _, _ in hist_defs:
        bins = array('d', np.asarray(edges, dtype=float).tolist())
        model = ROOT.RDF.TH1DModel(f'root_{key}', '', len(edges)-1, bins)
        # Weighted histogram supplies the PNG content.
        source = photon_rdf if key.startswith('photon_') else rdf
        hist_actions[key] = source.Histo1D(model, source_column(key), 'weight')
        # An additional UNWEIGHTED histogram counts the actual event entries,
        # irrespective of their possibly negative or non-unit weights.
        count_model = ROOT.RDF.TH1DModel(f'root_entries_{key}', '', len(edges)-1, bins)
        hist_entry_actions[key] = source.Histo1D(count_model, source_column(key))
    snapshot = None
    if photon_output is not None:
        # LAB convention: longitudinal momentum along the incoming beam;
        # retain original transverse axes and raw input vectors as well.
        ph = photon_rdf
        for comp, branch in (('E','scatteredEnergy'), ('Px','scatteredMomentumX'),
                             ('Py','scatteredMomentumY'), ('Pz','scatteredMomentumZ')):
            ph = ph.Define(f'__pr_{comp}_GeV', f'mcmule_analysis::proton_value(scatteredPID, {branch})')
        beta = math.sqrt(beam_energy_lab**2-MUON_MASS_GEV**2)/(beam_energy_lab+PROTON_MASS_GEV)
        gamma = 1/math.sqrt(1-beta*beta)
        pp = '(__pr_Pz_GeV * (beamMomentumZ >= 0.0 ? 1.0 : -1.0))'
        ph = ph.Define('__pr_E_lab_GeV', f'{gamma:.17g} * (__pr_E_GeV + {beta:.17g} * {pp})' if frame == 'cms' else '__pr_E_GeV')
        ph = ph.Define('__pr_p_lab_parallel_GeV', f'{gamma:.17g} * ({pp} + {beta:.17g} * __pr_E_GeV)' if frame == 'cms' else pp)
        aliases = {
            'eventID': 'rdfentry_', 'photonE': '__ph_E_lab_GeV',
            'photonPx': '__ph_Px_GeV', 'photonPy': '__ph_Py_GeV',
            'photonPz': '__ph_p_lab_parallel_GeV',
            'photonTheta': '__ph_theta_lab_mrad',
            'photonPhi': 'std::atan2(photonPy, photonPx)',
            'photonP': 'std::hypot(std::hypot(photonPx, photonPy), photonPz)',
            'muonE': '__E_lab_GeV', 'muonPx': '__mu_Px_GeV',
            'muonPy': '__mu_Py_GeV', 'muonPz': '__p_lab_parallel_GeV',
            'muonTheta': '__theta_lab_mrad', 'muonPhi': '__phi_mu_rad',
            'Q2_muon': '__Q2_GeV2',
            'protonE': '__pr_E_lab_GeV', 'protonPx': '__pr_Px_GeV',
            'protonPy': '__pr_Py_GeV', 'protonPz': '__pr_p_lab_parallel_GeV',
            'Q2_proton': f'protonPx*protonPx + protonPy*protonPy + protonPz*protonPz - (protonE-{PROTON_MASS_GEV:.17g})*(protonE-{PROTON_MASS_GEV:.17g})',
        }
        # Newer lhe_to_root.cpp trees already contain some of these names
        # (eventID, photonTheta, photonPhi, muonTheta). Redefine them so the
        # snapshot keeps exactly this script's definitions (eventID = chain entry).
        existing = set(str(c) for c in ph.GetColumnNames())
        for name, expr in aliases.items():
            ph = ph.Redefine(name, expr) if name in existing else ph.Define(name, expr)
        # Preserve raw vectors in addition to derived LAB quantities.
        columns = ROOT.std.vector('string')()
        for name in [*REQUIRED, *aliases]:
            columns.push_back(name)
        opts = ROOT.RDF.RSnapshotOptions()
        opts.fLazy = True
        opts.fMode = 'RECREATE'
        snapshot = ph.Snapshot('RealPhotonEvents', str(photon_output), columns, opts)
    counts = {k: float(v.GetValue()) for k, v in actions.items()}
    if snapshot is not None:
        snapshot.GetValue()
        f = ROOT.TFile.Open(str(photon_output), 'READ')
        if not f or f.IsZombie():
            raise RuntimeError('Photon snapshot could not be opened')
        tree = f.Get('RealPhotonEvents')
        if not tree or tree.GetEntries() != int(counts['one_photon_events']):
            f.Close()
            raise RuntimeError('Photon snapshot entry count mismatch')
        f.Close()

    if counts['multiple_photon_events']:
        print(f"WARNING: {int(counts['multiple_photon_events'])} events contain multiple photons; "
              "photon plots currently select events with exactly one photon.", flush=True)
    if not counts['one_photon_events']:
        print('WARNING: No photon (PID=22) found in ROOT output; photon plots will be empty. '
              'Check lhe_to_root.cpp to ensure photons are preserved.', flush=True)
    hists = {}
    for key, _, _, _, edges, _, _ in hist_defs:
        h = hist_actions[key].GetValue()
        val = np.array([float(h.GetBinContent(i)) for i in range(1, len(edges))])
        err = np.array([float(h.GetBinError(i)) for i in range(1, len(edges))])
        under = float(h.GetBinContent(0))
        over = float(h.GetBinContent(len(edges)))
        count_hist = hist_entry_actions[key].GetValue()
        count_bins = [int(round(count_hist.GetBinContent(i)))
                      for i in range(1, len(edges))]
        count_in_range = sum(count_bins)
        count_under = int(round(count_hist.GetBinContent(0)))
        count_over = int(round(count_hist.GetBinContent(len(edges))))
        count_filled = int(round(count_hist.GetEntries()))
        count_unaccounted = count_filled - (count_in_range + count_under + count_over)
        # RAW counts: do not divide by event totals or by bin widths.
        # For all weights == +1, each GetBinContent is an event count.
        # For a weighted sample, it is the sum of weights in that bin.
        hists[key] = {
            'edges': np.asarray(edges, dtype=float),
            'bin_count': val,
            'count_error': err,
            'underflow_weight': under,
            'overflow_weight': over,
            'in_range_weight': float(np.sum(val)),
            'zero_bin_count': int(np.count_nonzero(val == 0)),
            'events_filled': count_filled,
            'events_in_range': count_in_range,
            'events_underflow': count_under,
            'events_overflow': count_over,
            'events_unaccounted': count_unaccounted,
        }
    return counts, hists


def _plot_step(ax, hist, label, color):
    """Plot every bin, including true zero bins, with the actual bin edges."""
    edges = np.asarray(hist['edges'], dtype=float)
    y = np.asarray(hist['bin_count'], dtype=float)
    err = np.asarray(hist['count_error'], dtype=float)
    if not np.all(np.isfinite(y)) or not np.all(np.isfinite(err)):
        raise ValueError(f"Non-finite histogram bin content/error in {label}")

    # Step values at exactly zero show up on linear and symlog axes; never
    # replace them with NaNs or artificial tiny positive numbers.
    ax.stairs(y, edges, color=color, linewidth=1.65,
              baseline=None, label=label)

    low = np.maximum(0.0, y - err) if np.all(y >= 0) else y - err
    high = y + err
    ax.fill_between(edges, np.r_[low, low[-1]], np.r_[high, high[-1]],
                    step='post', color=color, alpha=0.15, linewidth=0)


def _choose_y_scale(requested, histograms):
    """Use symlog (rather than log) for wide ranges so zero bins stay visible."""
    if requested in ('linear', 'symlog'):
        return requested
    values = np.concatenate([np.asarray(h['bin_count'], dtype=float)
                             for h in histograms])
    values = values[np.isfinite(values)]
    positive = values[values > 0]
    if len(positive) < 2:
        return 'linear'
    if np.any(values < 0):
        return 'symlog'
    mid = float(np.median(positive))
    peak = float(np.max(positive))
    p10, p90 = np.quantile(positive, [0.1, 0.9])
    # Matplotlib's sym-log scale is logarithmic above one count, but unlike
    # plain log permits bins with y=0 and an exact baseline.
    if (mid > 0 and peak/mid >= 50) or (p10 > 0 and p90/p10 >= 25):
        return 'symlog'
    return 'linear'


def make_plot(key, title, xlabel, unit, requested_xscale, requested_yscale,
              data_by_sample, outdir, ylabel='Sum of event weights / bin'):
    fig, ax = plt.subplots(figsize=(7.6, 4.7), constrained_layout=True)
    available = [(tag, data_by_sample[tag][key]) for tag in ('LO', 'NLO')
                 if tag in data_by_sample and (not key.startswith('photon_') or tag == 'NLO')]
    if not available:
        plt.close(fig)
        return
    mode = _choose_y_scale(requested_yscale, [h for _, h in available])
    if mode == 'symlog':
        # All bins including 0 are representable; avoid huge 10^-300 axes.
        ax.set_yscale('symlog', linthresh=1, linscale=0.8)
        if any(np.any(h['bin_count'] < 0) for _, h in available):
            ax.axhline(0, color='0.45', lw=0.9, ls='--')
    for tag, current in available:
        _plot_step(ax, current, tag, COLORS[tag])
    if mode == 'symlog' and all(np.all(h['bin_count'] >= 0) for _, h in available):
        ax.set_ylim(bottom=0)
    # Only linear X binning is used. No variable-width or adaptive bins.
    ax.set_title(title, pad=11)
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)
    ax.set_xscale(requested_xscale)
    if requested_xscale == 'log' and available[0][1]['edges'][0] > 0:
        # Keep a strictly positive x range for logarithmic plots.
        ax.set_xlim(left=available[0][1]['edges'][0],
                    right=available[0][1]['edges'][-1])
    if not any(np.any(h['bin_count'] != 0) for _, h in available):
        ax.text(0.5, 0.5, 'No events in plotted range',
                transform=ax.transAxes, ha='center', va='center', fontsize=11)
    ax.legend(loc='best', frameon=False)
    fig.savefig(outdir / f'{key}.png', dpi=300, bbox_inches='tight')
    plt.close(fig)


def write_summary(filename, batch_files, stats, hist_defs, frame_info, missing, hist_data):
    out = []
    out.append('McMule ROOT event-generation analysis')
    out.append('Generated ' + datetime.now().strftime('%Y-%m-%d %H:%M:%S'))
    out.append('NOTE: Main figures show raw events per bin (sum of event weights).')
    out.append('NOTE: Raw event counts do not encode the absolute NLO/LO cross-section normalization.')
    out.append('NOTE: Bin errors use ROOT sumw2 and neglect correlations between events.')
    out.append('NOTE: For NLO CU, generator resampling can correlate events.')
    out.append('NOTE: Photon plots count ONLY events with exactly one final-state photon (PDG 22).')
    out.append('NOTE: Counts are physical event numbers ONLY if all input weights equal +1.')
    out.append('NOTE: If every ROOT event has weight +1, per-job generator normalization may still differ.')
    out.append('')
    for tag in ('NLO', 'LO'):
        if tag not in stats:
            continue
        s = stats[tag]
        info = frame_info[tag]
        out.append(f'{tag} sample')
        out.append(f'  input frame (auto-detected or chosen): {info[0]}')
        out.append(f'  first beamEnergy: {info[1]:.10g} GeV')
        out.append(f'  first beamMomentumZ: {info[2]:.10g} GeV')
        out.append(f"  entries: {int(s['entries']):,}")
        out.append(f"  negative-weight entries: {int(s['negative_entries']):,}")
        out.append(f"  sum(weights): {s['weight_sum']:.12g}")
        out.append(f"  sum(abs(weights)): {s['weight_abs_sum']:.12g}")
        out.append(f"  sum(weights^2): {s['weight_sq_sum']:.12g}")
        out.append(f"  zero-photon events: {int(s['zero_photon_events']):,}")
        out.append(f"  one-photon events (PDG 22): {int(s['one_photon_events']):,}")
        out.append(f"  >1 photon events: {int(s['multiple_photon_events']):,}")
        out.append(f"  one-photon sum(weights): {s['one_photon_weight_sum']:.12g}")
        out.append(f"  events with >=1 real photon (PDG 22): "
                   f"{int(s['one_photon_events'] + s['multiple_photon_events']):,}")
        out.append('  INTERPRETATION: 0-photon events contain Born-like/virtual '
                   'contributions and cannot be split into LO vs virtual from this ROOT file.')
        out.append('  INTERPRETATION: photon presence identifies final-state real-emission '
                   'events, not a separately normalized physical R cross section.')
        out.append('')
        out.append('  ROOT inputs:')
        for f, n in batch_files[tag]:
            out.append(f'    {n:>12,}  {f}')
        out.append('')
    if missing:
        out.append('MISSING FILES (allowed by --allow-missing):')
        out.extend('  ' + x for x in missing)
        out.append('')
    out.append('Total events per histogram (no per-bin output):')
    out.append('  In range = events actually shown in the plot; outside = underflow + overflow.')
    out.append('  Photon plots select events with exactly one final-state photon (PID 22).')
    for key, title, xlabel, unit, edges, _, _ in hist_defs:
        for tag in ('NLO', 'LO'):
            if tag not in hist_data:
                continue
            h = hist_data[tag][key]
            outside = h['events_underflow'] + h['events_overflow']
            out.append(f"  {tag} | {title}: {h['events_in_range']:,} events in plot"
                       f" (total selected: {h['events_filled']:,}; outside: {outside:,})")
    out.append('')
    filename.write_text('\n'.join(out)+'\n', encoding='utf-8')


def main(argv=None):
    p = argparse.ArgumentParser(description=__doc__)
    inputs = p.add_mutually_exclusive_group(required=True)
    inputs.add_argument('--batch-root', type=Path)
    inputs.add_argument('--nlo-root', nargs='+', type=Path, help='One or more existing NLO ROOT files')
    p.add_argument('--n-runs', type=int, default=10)
    p.add_argument('--base-seed', type=int, default=10000)
    p.add_argument('--xi', default='10')
    p.add_argument('--prefix', default='08_11_evtgen_xi1_NLO_CU')
    p.add_argument('--outdir', required=True, type=Path)
    p.add_argument('--lo-root', type=Path, default=None,
                   help='Optional existing LO ROOT for shape-only comparison')
    p.add_argument('--input-frame', choices=['auto', 'cms', 'lab'], default='auto')
    p.add_argument('--lo-frame', choices=['auto', 'cms', 'lab'], default='auto')
    p.add_argument('--beam-energy-gev', type=float, default=100.0)
    p.add_argument('--allow-missing', action='store_true')
    p.add_argument('--threads', type=int, default=1)
    args = p.parse_args(argv)

    try:
        import ROOT
    except ImportError as e:
        raise SystemExit('PyROOT is unavailable. Run the command script after sourcing LCG_109. '
                         'Check: python3 -c "import ROOT, numpy, matplotlib"') from e
    ROOT.gROOT.SetBatch(True)
    # PyROOT RDataFrame sees the outgoing particles as RVec<double>.
    # Compile a scalar muon extractor once before building the RDF graph.
    ROOT.gInterpreter.Declare(r"""
        #include <ROOT/RVec.hxx>
        #include <cmath>
        #include <algorithm>
        #include <stdexcept>
        #include <cstddef>
        namespace mcmule_analysis {
          inline int count_photons(const ROOT::VecOps::RVec<int>& pids) {
            return static_cast<int>(std::count(pids.begin(), pids.end(), 22));
          }
          template <typename T>
          double photon_value(const ROOT::VecOps::RVec<int>& pids,
                              const ROOT::VecOps::RVec<T>& values) {
            if (pids.size() != values.size()) {
              throw std::runtime_error("McMule ROOT: photon PID/values size mismatch");
            }
            bool found = false;
            double result = 0.0;
            for (std::size_t i = 0; i < pids.size(); ++i) {
              if (pids[i] == 22) {
                if (found) throw std::runtime_error("Multiple photons in 1-photon selection");
                found = true;
                result = static_cast<double>(values[i]);
              }
            }
            if (!found) throw std::runtime_error("Photon requested in photon-free event");
            return result;
          }
          template <typename T>
          double proton_value(const ROOT::VecOps::RVec<int>& pids,
                              const ROOT::VecOps::RVec<T>& values) {
            if (pids.size() != values.size()) throw std::runtime_error("Proton vector size mismatch");
            int n = 0; double v = 0;
            for (std::size_t i=0; i<pids.size(); ++i) if (pids[i]==2212) { ++n; v=values[i]; }
            if (n != 1) throw std::runtime_error("Expected exactly one outgoing proton");
            return v;
          }
          template <typename T>
          double muon_value(const ROOT::VecOps::RVec<int>& pids,
                            const ROOT::VecOps::RVec<T>& values) {
            if (pids.size() != values.size()) {
              throw std::runtime_error("McMule ROOT: PID and four-vector sizes differ");
            }
            bool found = false;
            double result = 0.0;
            for (std::size_t i = 0; i < pids.size(); ++i) {
              if (std::abs(pids[i]) == 13) {
                if (found) {
                  throw std::runtime_error("McMule ROOT: multiple final-state muons");
                }
                found = true;
                result = static_cast<double>(values[i]);
              }
            }
            if (!found) {
              throw std::runtime_error("McMule ROOT: no final-state muon in event");
            }
            return result;
          }
        }
    """)
    if args.threads > 1:
        ROOT.EnableImplicitMT(args.threads)

    args.outdir.mkdir(parents=True, exist_ok=True)
    nlo_files, missing = [], []
    for i in range(1, args.n_runs+1) if args.batch_root else []:
        seed = args.base_seed + i-1
        filepath = (args.batch_root / f'part_{i}' / 'out' /
                    f'{args.prefix}_part_{i}_seed{seed}.root')
        if filepath.is_file():
            nlo_files.append(str(filepath))
        else:
            missing.append(str(filepath))
    if args.nlo_root:
        nlo_files = [str(x) for x in args.nlo_root]
    if missing and not args.allow_missing:
        raise SystemExit(f'Missing {len(missing)} of {args.n_runs} NLO ROOT files. '
                         'Run conversion first, or pass --allow-missing to process '
                         'the completed parts only.\n' + '\n'.join(missing))
    if not nlo_files:
        raise SystemExit('No NLO ROOT files found at the configured paths.')

    sample_files = {'NLO': nlo_files}
    if args.lo_root is not None:
        sample_files['LO'] = [str(args.lo_root)]
    hist_defs = observables()
    data, filelist, stats, frame_info = {}, {}, {}, {}
    for tag in ('NLO', 'LO'):
        if tag not in sample_files:
            continue
        print(f'\nLoading {tag}: {len(sample_files[tag])} ROOT files', flush=True)
        chain, files, energy, pz = check_and_chain(ROOT, sample_files[tag], tag)
        frame_arg = args.input_frame if tag == 'NLO' else args.lo_frame
        frame = choose_frame(frame_arg, energy, args.beam_energy_gev)
        print(f'  first incoming muon: E={energy:.7g} GeV, pz={pz:.7g} GeV; '
              f'input frame={frame}', flush=True)
        rdf = build_rdf(ROOT, chain, frame, args.beam_energy_gev)
        counts, plots = extract_hist(ROOT, rdf, hist_defs, frame, args.beam_energy_gev,
                                     args.outdir/'mcmule_evtgen.partial.root' if tag == 'NLO' else None)
        print(f"  entries={int(counts['entries']):,}, "
              f"negative={int(counts['negative_entries']):,}, "
              f"sum(weights)={counts['weight_sum']:.10g}, "
              f"one-photon events={int(counts['one_photon_events']):,}", flush=True)
        stats[tag], data[tag], filelist[tag] = counts, plots, files
        frame_info[tag] = (frame, energy, pz)

    # If every event has weight +1 (unweighted CU sample), bins are plain event counts.
    unit_weights = all(st['negative_entries'] == 0 and st['weight_sum'] == st['entries']
                       and st['weight_sq_sum'] == st['entries'] for st in stats.values())
    ylabel = 'Events / bin' if unit_weights else 'Sum of event weights / bin'
    for key, title, xlabel, unit, edges, xscale, yscale in hist_defs:
        make_plot(key, title, xlabel, unit, xscale, yscale, data, args.outdir, ylabel)
    # Store histograms and photon event tree together in ONE analysis ROOT file.
    partial = args.outdir/'mcmule_evtgen.partial.root'
    f = ROOT.TFile.Open(str(partial), 'UPDATE')
    if not f or f.IsZombie():
        raise RuntimeError('Cannot update analysis ROOT file')
    try:
        for tag, plots in data.items():
            directory = f.mkdir(tag)
            directory.cd()
            for key, h in plots.items():
                if tag == 'LO' and key.startswith('photon_'):
                    continue
                edges = h['edges']
                hist = ROOT.TH1D(key, key, len(edges)-1, array('d', edges.tolist()))
                hist.Sumw2()
                for i, (v, e) in enumerate(zip(h['bin_count'], h['count_error']), 1):
                    hist.SetBinContent(i, float(v))
                    hist.SetBinError(i, float(e))
                hist.SetBinContent(0, h['underflow_weight'])
                hist.SetBinContent(len(edges), h['overflow_weight'])
                hist.SetEntries(h['events_filled'])
                hist.Write()
            f.cd()
        ROOT.TNamed('kinematics_convention', 'GeV; theta mrad; phi rad; Q2 GeV^2; longitudinal axis along beam; raw vectors in input frame; eventID is chain entry').Write()
        ROOT.TNamed('NLO_input_files', '\n'.join(nlo_files)).Write()
        ROOT.TNamed('NLO_input_frame', frame_info['NLO'][0]).Write()
    finally:
        f.Close()
    partial.replace(args.outdir/'mcmule_evtgen.root')
    write_summary(args.outdir/'event_statistics.txt', filelist,
                  stats, hist_defs, frame_info, missing, data)
    # Do not emit CSV files. Remove only the two legacy CSVs produced by the
    # earlier versions of THIS analysis, so old counts cannot be mistaken for
    # current results when reusing the same output directory.
    for obsolete in ('histogram_bin_counts.csv', 'histogram_event_counts.csv'):
        (args.outdir / obsolete).unlink(missing_ok=True)
    print('\nWrote figures to', args.outdir, flush=True)
    print('Events per histogram and sample statistics:',
          args.outdir/'event_statistics.txt', flush=True)
    print('Created', len(hist_defs), 'main PNG figures', flush=True)
    if 'LO' in data:
        print('LO/NLO main plots show RAW event counts, not rate-normalized.', flush=True)
    else:
        print('NLO raw event counts only; empty bins shown at y=0.', flush=True)


if __name__ == '__main__':
    main()