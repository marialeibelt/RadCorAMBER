#!/usr/bin/env python3
"""Determine the first occupied ultra-fine photon-energy bin without allocating a histogram.

Input: McMule LHE->ROOT event-level ROOT files, TTree Output, vector branches
scatteredPID and scatteredEnergy (GeV). Photon is PDG 22.

No ROOT histogram is constructed: scan photon energies and report the bin
containing the smallest nonnegative photon energy, plus the smallest positive
photon energy for physical diagnostics. All reported energies are in MeV.
"""
from __future__ import annotations

import argparse
from decimal import Decimal, ROUND_FLOOR, localcontext
import glob
from pathlib import Path
import sys


def bin_edges(energy_mev: float, width_mev: str) -> tuple[int, Decimal, Decimal]:
    """Return zero-based bin index and [lower,upper) edges in MeV.

    Decimal(str(...)) prevents accidental display rounding in the final report.
    Underlying ROOT floating-point precision still determines actual accuracy.
    """
    with localcontext() as ctx:
        ctx.prec = 50
        energy = Decimal(str(energy_mev))
        width = Decimal(width_mev)
        if energy < 0 or width <= 0:
            raise ValueError('Energy must be nonnegative; bin width must be positive')
        index = int((energy / width).to_integral_value(rounding=ROUND_FLOOR))
        lower = width * index
        return index, lower, lower + width


def main() -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input-glob', required=True,
                        help='McMule LHE->ROOT event-file glob; may include **')
    parser.add_argument('--outdir', required=True)
    parser.add_argument('--max-mev', type=float, default=50.0,
                        help='Exclusive upper photon energy bound in MeV')
    parser.add_argument('--bin-width-mev', default='1e-10',
                        help='Hypothetical linear-bin width in MeV; no actual bins allocated')
    args = parser.parse_args()
    try:
        width = Decimal(args.bin_width_mev)
    except Exception:
        parser.error('BIN_WIDTH_MEV must be a positive decimal number')
    if not (width.is_finite() and width > 0 and args.max_mev > 0):
        parser.error('Require positive finite bin width and max energy')

    paths = sorted(p for p in set(glob.glob(args.input_glob, recursive=True)) if Path(p).is_file())
    if not paths:
        parser.error(f'No ROOT files match {args.input_glob!r}')
    out = Path(args.outdir)
    out.mkdir(parents=True, exist_ok=True)

    try:
        import ROOT
    except ImportError as exc:
        raise RuntimeError('PyROOT unavailable: source ROOT/LCG environment first') from exc
    ROOT.gROOT.SetBatch(True)
    chain = ROOT.TChain('Output')
    for f in paths:
        if chain.Add(f) <= 0:
            raise RuntimeError(f'Unable to add ROOT file: {f}')
    entries = int(chain.GetEntries())
    if entries <= 0:
        raise RuntimeError('No events in ROOT tree Output')
    chain.LoadTree(0)
    for branch in ['scatteredPID', 'scatteredEnergy']:
        if not chain.GetBranch(branch):
            raise RuntimeError(f'ROOT branch {branch} missing; use McMule event-level ROOT, not PHAST')

    # Uses ROOT's vector-aware Draw selection: one row per selected photon.
    # This never creates any histogram; only stores energies of selected photons.
    # SetEstimate ensures the drawing buffer will not truncate the extracted values.
    chain.SetEstimate(max(1_000_000, entries))
    max_gev = args.max_mev / 1000.0
    condition = f'scatteredPID==22 && scatteredEnergy>=0 && scatteredEnergy<{max_gev:.17g}'
    print(f'Processing {entries:,} ROOT events across {len(paths)} input files...')
    count = int(chain.Draw('scatteredEnergy*1000', condition, 'goff'))
    if count < 0:
        raise RuntimeError('ROOT Draw failed')
    if count > chain.GetEstimate():
        raise RuntimeError('ROOT value-buffer limit exceeded')
    if count == 0:
        lines = [
            'McMule generated real-photon energy: first occupied hypothetical bin',
            f'Input ROOT files: {len(paths)}',
            f'ROOT events: {entries}',
            f'Energy range: 0 <= E_gamma < {args.max_mev:g} MeV',
            f'Hypothetical linear bin width: {width} MeV',
            'Selected photons in this range: 0',
            'First occupied bin: NONE (no photons in selected range)',
        ]
    else:
        vals = chain.GetV1()
        # One pass; no dense bins, no huge arrays of empty bins.
        min_all = float('inf')
        min_positive = float('inf')
        zero_count = 0
        valid_count = 0
        for i in range(count):
            v = float(vals[i])
            if not (0 <= v < args.max_mev):
                continue
            valid_count += 1
            if v < min_all:
                min_all = v
            if v == 0.0:
                zero_count += 1
            elif v < min_positive:
                min_positive = v
        if valid_count == 0:
            raise RuntimeError('Selected photon energies are outside range after unit conversion')
        index, lower, upper = bin_edges(min_all, str(width))
        lines = [
            'McMule generated real-photon energy: first occupied hypothetical bin',
            f'Input ROOT files: {len(paths)}',
            f'ROOT events: {entries}',
            f'Energy range: 0 <= E_gamma < {args.max_mev:g} MeV',
            f'Hypothetical linear bin width: {width} MeV',
            f'Photons in range: {valid_count}',
            f'Photons exactly at E_gamma=0: {zero_count}',
            f'Minimum E_gamma including zero: {min_all:.17g} MeV',
            f'First occupied bin (index starting at 0): {index}',
            f'First occupied bin: [{lower}, {upper}) MeV',
        ]
        if min_positive != float('inf'):
            pi, pl, pu = bin_edges(min_positive, str(width))
            lines.extend([
                f'Minimum strictly positive E_gamma: {min_positive:.17g} MeV',
                f'First positive-energy bin (index starting at 0): {pi}',
                f'First positive-energy bin: [{pl}, {pu}) MeV',
            ])
        else:
            lines.append('No strictly positive photon energy in selected range')
        lines.extend([
            'The first occupied bin is determined without building any histogram.',
            'Underlying energy values are ROOT floating-point numbers; this does not create new numerical precision.',
            'This is the sample minimum, not a physical photon-energy cutoff.',
        ])
    report = out / 'first_photon_energy_bin.txt'
    report.write_text('\n'.join(lines) + '\n', encoding='utf-8')
    print('\n'.join(lines))
    print(f'\nSaved: {report}')
    return 0


if __name__ == '__main__':
    sys.exit(main())
