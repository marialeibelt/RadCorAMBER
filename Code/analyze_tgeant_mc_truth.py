#!/usr/bin/env python3
"""Stream AMBER TGEANT ASCII .tgeant.gz to MC-truth/MC-hit ROOT histograms.

Test parsing without PyROOT with --summary-only. The normal mode writes a
ROOT file through PyROOT (available from the user's LCG ROOT environment).
The TPC pad ADC waveforms are intentionally skipped, NOT summed as energies.
"""
import argparse
import gzip
import math
import sys
from collections import Counter

MU_MASS = 0.105658375  # GeV
PROTON_MASS = 0.938272013  # GeV; mass quoted in uploaded TGEANT beamdata
GEV_PER_MEV = 0.001
MM_PER_CM = 10.0


def fresh_event():
    return {
        'vertex': None, 'particles': [], 'n_hit': Counter(), 'edep': Counter(),
        'n_tpc_proton': 0, 'edep_tpc_proton': 0.0, 'n_adc': 0,
    }


def read_line(fh):
    line = fh.readline()
    if not line:
        raise ValueError('Unexpected end of gzip file inside event block')
    return line.strip()


def parse_beamdata(fh, header, event):
    # BEAMDATA: x y z t [other run-specific fields]; coordinates in mm.
    fields = header.split()
    if len(fields) < 4:
        raise ValueError('Malformed #BEAMDATA header')
    event['vertex'] = tuple(float(x) / MM_PER_CM for x in fields[1:4])
    # A long auxiliary beam line (e.g. 100 zero entries in the uploaded file).
    read_line(fh)
    n = int(read_line(fh))
    if n < 0 or n > 200:
        raise ValueError(f'Unexpected beam particle count: {n}')
    particles = []
    for _ in range(n):
        p = read_line(fh).split()
        # Uploaded format: prefix, PDG, 3 placeholders, px, py, pz, E, mass,
        # where (px,py,pz,E,mass) are in GeV.
        if len(p) < 10:
            raise ValueError(f'Bad BEAMDATA particle line: {p!r}')
        particles.append({
            'pdg': int(p[1]),
            'px': float(p[5]), 'py': float(p[6]), 'pz': float(p[7]),
            'energy': float(p[8]), 'mass': float(p[9]),
        })
    event['particles'] = particles
    # TGEANT's beamdata also contains a list of Geant trajectories; these
    # are NOT extra generator-level external particles.
    n_tracks = int(read_line(fh))
    if n_tracks < 0 or n_tracks > 1_000_000:
        raise ValueError(f'Unexpected number of beam trajectories: {n_tracks}')
    for _ in range(n_tracks):
        read_line(fh)


def parse_hits(fh, header, event):
    block, n = header.split()[:2]
    count = int(n)
    if count < 0 or count > 10_000_000:
        raise ValueError(f'Unexpected number of {block} hits: {count}')
    event['n_hit'][block] += count
    if block == '#TPCPadHit':
        for _ in range(count):
            metadata = read_line(fh).split()
            # Last column = number of ADC values in the NEXT line.
            if metadata:
                try: event['n_adc'] += int(metadata[-1])
                except ValueError: pass
            read_line(fh)  # Full waveform skipped, can be >10000 characters.
        return
    for _ in range(count):
        values = read_line(fh).split()
        # HIT sections list eDep at zero-based column 8 (MeV).
        if len(values) < 9:
            raise ValueError(f'Bad {block} hit (expected >=9 fields): {values!r}')
        edep_mev = float(values[8])
        event['edep'][block] += edep_mev
        if block == '#TPCIonization' and int(values[4]) == 2212:
            event['n_tpc_proton'] += 1
            event['edep_tpc_proton'] += edep_mev


def event_stream(filenames, max_events=None):
    n_processed = 0
    for name in filenames:
        print(f'Reading {name}', file=sys.stderr)
        with gzip.open(name, 'rt', errors='replace') as f:
            event = None
            for line in f:
                if not line or line[0] == '*':
                    continue
                if line.startswith('#EVENT\t') or line.startswith('#EVENT '):
                    if event is not None:
                        raise ValueError(f'Encountered #EVENT before #EVENT_END in {name}')
                    event = fresh_event()
                elif line.startswith('#EVENT_END'):
                    if event is None:
                        raise ValueError(f'Unexpected #EVENT_END in {name}')
                    yield event
                    n_processed += 1
                    event = None
                    if max_events is not None and n_processed >= max_events:
                        return
                elif event is not None and line.startswith('#BEAMDATA'):
                    parse_beamdata(f, line, event)
                elif event is not None and line.startswith(('#TRACKING', '#TRIGGER', '#CALORIMETER', '#TPCIonization', '#TPCPadHit')):
                    parse_hits(f, line, event)
                elif line.startswith('#') and event is not None:
                    raise ValueError(f'Unknown block in {name}: {line[:160]!r}')
            if event is not None:
                raise ValueError(f'Unfinished event at end of {name}')


def p_norm(p):
    return math.sqrt(p['px'] ** 2 + p['py'] ** 2 + p['pz'] ** 2)


def truth_muon(event):
    muons = [p for p in event['particles'] if abs(p['pdg']) == 13]
    # VERIFIED FOR THE UPLOADED LO EXAMPLE: the FIRST is the propagated
    # incoming beam muon, the SECOND is the generated outgoing muon.
    # If the simulated process has a different record layout, adapt here.
    return muons[1] if len(muons) >= 2 else None


def q2_nominal(p, beam_energy):
    pbeam = math.sqrt(beam_energy ** 2 - MU_MASS ** 2)
    # Q² = 2 (Ein Eout - pin·pout - m_mu²); nominal +z beam, NOT event-by-event.
    return 2.0 * (beam_energy * p['energy'] - pbeam * p['pz'] - MU_MASS ** 2)


def histograms(ROOT):
    b = {}
    def h1(key, title, bins, low, high):
        b[key] = ROOT.TH1D(key, title, bins, low, high)
        b[key].Sumw2()
    def h2(key, title, nx, lx, hx, ny, ly, hy):
        b[key] = ROOT.TH2D(key, title, nx, lx, hx, ny, ly, hy)
        b[key].Sumw2()
    # Identical axis definitions to the PHAST module, but DIFFERENT level!
    h1('h_vertex_x', 'MC truth interaction x;x [cm];Events', 200, -20, 20)
    h1('h_vertex_y', 'MC truth interaction y;y [cm];Events', 200, -20, 20)
    h1('h_vertex_z', 'MC truth interaction z;z [cm];Events', 300, -500, 500)
    h2('h_vertex_xy', 'MC truth vertex x/y;x [cm];y [cm]',120,-12,12,120,-12,12)
    h1('h_muon_p', 'Generated outgoing muon momentum;p [GeV/c];Events',240,0,120)
    h1('h_muon_E', 'Generated outgoing muon energy;E [GeV];Events',240,0,120)
    h1('h_muon_theta', 'Generated outgoing muon angle to +z;theta [mrad];Events',240,0,12)
    h1('h_muon_theta_fine', 'Generated outgoing muon angle (zoom);theta [mrad];Events',250,0,2)
    h1('h_muon_phi', 'Generated outgoing muon azimuth;phi [rad];Events',128,-math.pi,math.pi)
    h1('h_q2_nominal_beam_approx','MC nominal +z beam Q2 approximation;Q2 [GeV2];Events',250,0,0.05)
    h1('h_q2_nominal_beam_approx_fine','MC nominal +z beam Q2 approximation (zoom);Q2 [GeV2];Events',250,0,0.001)
    h1('h_nu_nominal_beam_approx','MC nominal beam energy difference, NOT eDep;Ebeam-Eout [GeV];Events',220,-10,100)
    h2('h_muon_E_vs_theta','Generated outgoing muon;theta [mrad];E [GeV]',120,0,12,160,0,120)
    h2('h_q2_vs_theta_nominal_beam_approx','Nominal beam Q2 vs theta;theta [mrad];Q2 [GeV2]',120,0,12,120,0,0.05)
    # MC-truth-only and GEANT hit-level histograms (not directly PHAST reco!)
    h1('h_n_truth_particles','Particles in beamdata;N;Events',16,-0.5,15.5)
    h1('h_n_truth_photons','Generated photon count;N;Events',10,-0.5,9.5)
    h1('h_truth_photon_E','Generated photons;E_gamma [GeV];Photons',240,0,100)
    h1('h_q2_proton_truth','Proton-side Q2 (initial proton at rest);Q2 [GeV2];Events',250,0,0.05)
    h1('h_q2_proton_truth_fine','Proton-side Q2 (zoom, proton initially at rest);Q2 [GeV2];Events',250,0,0.001)
    for key, lbl in [('tracking','TRACKING'),('trigger','TRIGGER'),('calo','CALORIMETER'),('tpc_ion','TPCIonization'),('tpc_pad','TPCPadHit')]:
        bins = 500 if key in ('tracking', 'tpc_pad') else 100
        xmax = 500 if key in ('tracking', 'tpc_pad') else 100
        h1('h_n_'+key+'_hits',lbl+' hit count;N hits;Events',bins,-0.5,xmax-0.5)
    h1('h_n_tpc_adc_samples','TPCPadHit waveform samples;N ADC values;Events',400,0,400000)
    for key,label,xmax in [('tracking','Tracking',50),('trigger','Trigger',100),('calo','Calorimeter',50000),('tpc_ion','TPC ionization',100)]:
        h1('h_'+key+'_edep_sum_MeV',label+' GEANT eDep sum (not reconstructed);Edep [MeV];Events',300,0,xmax)
    h1('h_tpc_proton_edep_sum_MeV','TPC proton GEANT eDep sum;Edep [MeV];Events',300,0,10)
    h1('h_n_tpc_proton_hits','TPC ionization hits of PDG 2212;N;Events',100,-0.5,99.5)
    return b


def fill(h, name, *values):
    h[name].Fill(*values)


def fill_event(event,h,beam_energy,muon_min_p):
    if event['vertex'] is not None:
        vx,vy,vz = event['vertex']
        for name,v in [('h_vertex_x',vx),('h_vertex_y',vy),('h_vertex_z',vz)]:fill(h,name,v)
        fill(h,'h_vertex_xy',vx,vy)
    particles = event['particles']
    fill(h,'h_n_truth_particles',len(particles))
    photons = [p for p in particles if p['pdg']==22]
    fill(h,'h_n_truth_photons',len(photons))
    for p in photons:fill(h,'h_truth_photon_E',p['energy'])
    p = truth_muon(event)
    if p is not None and p_norm(p) > muon_min_p:
        theta=math.atan2(math.hypot(p['px'],p['py']),p['pz'])*1000.0
        phi=math.atan2(p['py'],p['px'])
        pmag=p_norm(p)
        for key,val in [('h_muon_p',pmag),('h_muon_E',p['energy']),('h_muon_theta',theta),('h_muon_phi',phi)]:fill(h,key,val)
        fill(h,'h_muon_theta_fine',theta)
        q2=q2_nominal(p,beam_energy)
        fill(h,'h_q2_nominal_beam_approx_fine',q2)
        fill(h,'h_q2_nominal_beam_approx',q2)
        fill(h,'h_nu_nominal_beam_approx',beam_energy-p['energy'])
        fill(h,'h_muon_E_vs_theta',theta,p['energy'])
        fill(h,'h_q2_vs_theta_nominal_beam_approx',theta,q2)
    # Applies only when a recoil proton is in beamdata and was initially at rest.
    protons=[p for p in particles if p['pdg']==2212]
    if len(protons)==1:
        q2p=2*PROTON_MASS*(protons[0]['energy']-PROTON_MASS)
        fill(h,'h_q2_proton_truth',q2p)
        fill(h,'h_q2_proton_truth_fine',q2p)
    for key, block in [('tracking','#TRACKING'),('trigger','#TRIGGER'),('calo','#CALORIMETER'),('tpc_ion','#TPCIonization'),('tpc_pad','#TPCPadHit')]:
        fill(h,'h_n_'+key+'_hits',event['n_hit'][block])
    fill(h,'h_n_tpc_adc_samples',event['n_adc'])
    for key, block in [('tracking','#TRACKING'),('trigger','#TRIGGER'),('calo','#CALORIMETER'),('tpc_ion','#TPCIonization')]:
        fill(h,'h_'+key+'_edep_sum_MeV',event['edep'][block])
    fill(h,'h_tpc_proton_edep_sum_MeV',event['edep_tpc_proton'])
    fill(h,'h_n_tpc_proton_hits',event['n_tpc_proton'])


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument('--output', help='Output ROOT file (required unless --summary-only)')
    p.add_argument('--max-events', type=int, default=0, help='0 means all events')
    p.add_argument('--nominal-beam-energy', type=float, default=100., help='GeV, for DIAGNOSTIC Q2 only')
    p.add_argument('--muon-min-p', type=float, default=70., help='GeV/c; common muon candidate cut with PHAST')
    p.add_argument('--summary-only', action='store_true', help='Parse/validate without ROOT; does not produce histograms')
    p.add_argument('inputs',nargs='+',help='.tgeant.gz files')
    args=p.parse_args()
    if args.max_events<0: p.error('--max-events must be >=0')
    if not args.summary_only and not args.output: p.error('--output required unless --summary-only')
    h=None; ROOT=None; out=None
    if not args.summary_only:
        import ROOT as ROOT_module
        ROOT=ROOT_module
        ROOT.gROOT.SetBatch(True)
        out=ROOT.TFile.Open(args.output,'RECREATE')
        if not out or out.IsZombie(): raise OSError('Cannot create '+args.output)
        out.mkdir('MC_TRUTH').cd()
        h=histograms(ROOT)
    n=0; section_counts=Counter(); muE=[]; sample_q2=[]; first_edep=None
    for ev in event_stream(args.inputs, args.max_events if args.max_events else None):
        n+=1
        if h is not None: fill_event(ev,h,args.nominal_beam_energy,args.muon_min_p)
        for k,v in ev['n_hit'].items():section_counts[k]+=v
        mu=truth_muon(ev)
        if mu and len(muE)<5:
            muE.append(mu['energy'])
            sample_q2.append(q2_nominal(mu,args.nominal_beam_energy))
        if first_edep is None:first_edep=dict(ev['edep'])
        if n%1000==0:print(f'  processed {n} events',file=sys.stderr)
    if not n:raise ValueError('No completed events found')
    if out:
        ROOT.TNamed('MC_TRUTH_readme','TGEANT beamdata + GEANT hits; NOT CORAL reconstructed; beamdata position converted mm->cm; hit eDep in MeV; nominal +z beam for approximate Q2; unweighted events').Write()
        out.Write()
        out.Close()
        print('Saved ROOT histogram file: '+args.output)
    print('Events:',n)
    print('First outgoing muon energies [GeV]:',muE)
    print('First nominal-beam Q2 [GeV2]:',sample_q2)
    print('First event eDep section sums [MeV]:',first_edep)
    print('Total number of detector hits:',dict(section_counts))


if __name__=='__main__':
    main()
