# AMBER PRM: TGEANT MC truth/hits ↔ CORAL/PHAST reconstructed histogram comparison

This package keeps the **event reconstruction** in CORAL, processes the reconstructed **mDST** with PHAST, and independently produces MC-truth/hit-level histograms from the TGEANT `*.tgeant.gz` ASCII output. A third ROOT macro plots the **common variables** from both histogram ROOT files.

**Do not process TGEANT hits as if they were reconstructed PHAST tracks.** Histograms with matching names can be overlaid, but in general are at different analysis levels. MC-vs-reco overlaid shapes require similar event selections and identical coordinate systems. For actual efficiency/resolution, process the simulated hits with CORAL, match reconstructed objects to their Monte Carlo truth within the **same events**, and include lost events.

## Contents

- `Code/analyze_tgeant_mc_truth.py`: streaming parser for TGEANT ASCII `.tgeant.gz`; writes `MC_TRUTH` directory with ROOT histograms using PyROOT.
- `scripts/command_analyze_tgeant_mc_truth.sh` and `inputs/input_analyze_tgeant_mc_truth.conf`: TGEANT analysis runner/config.
- `Code/u12345_AMBER_PRM_histograms.cc`: PHAST UserEvent module for CORAL-reconstructed mDST histograms, stored in `AMBER_PRM` ROOT directory. It assumes a negative outgoing track with p>70 GeV as the muon *candidate*, not positive PID.
- `scripts/command_analyze_phast_prm.sh`, `Code/plot_analyze_phast_prm.C`, `inputs/input_analyze_phast_prm.conf`: PHAST histogram plotting, not PHAST running/compiling.
- `Code/plot_compare_tgeant_phast.C`: overlay common 1D hists normalized individually to unit area, and show common 2D hists side-by-side; saves PDF/PNG.
- `scripts/command_compare_tgeant_phast.sh` and `inputs/input_compare_tgeant_phast.conf`: comparison runner/config.
- `sample_validation.txt`: result of parsing the supplied 1,000-event example in summary mode.

## Run TGEANT analysis on Themis

1. Edit `inputs/input_analyze_tgeant_mc_truth.conf` to set your actual `input_glob` and `mc_histogramfile`. It defaults to a **single sample file pattern** (the precise server location may differ). `max_files=0` means all matching files. You may use `--max-events` in the Python CLI or config `max_events` for a quick test.
2. From this package directory run:

   ```bash
   bash scripts/command_analyze_tgeant_mc_truth.sh
   ```

   This sources your existing LCG_109/ROOT environment, checks PyROOT, streams every matching gzip file, and writes one ROOT histogram file. The large TPC Pad ADC waveform rows are **skipped**, not materialized as event-wise Python objects, to keep memory use low. ADC samples are counted, not interpreted as energy.

3. Without CERN ROOT/PyROOT you can validate any file with:

   ```bash
   python3 Code/analyze_tgeant_mc_truth.py --summary-only --max-events 10 /path/file.tgeant.gz
   ```

## PHAST side

Use the `u12345_AMBER_PRM_histograms.cc` as a PHAST UserEvent in your local PHAST build, according to your actual build/user-slot procedure, and run PHAST on CORAL-generated mDSTs. This is **not** a standalone C++ ROOT program. It creates an `AMBER_PRM` histogram directory in the PHAST histogram ROOT output. Adapt the candidate selection if you have reliable scattered-muon PID in your PRM reconstruction. Methods were not compiled against your local PHAST installation and could need adaptation.

Optional PHAST plotting (after the PHAST histogram ROOT file exists):

```bash
bash scripts/command_analyze_phast_prm.sh
```

## Compare histogram ROOT files

1. Enter the existing MC and PHAST histogram ROOT paths into `inputs/input_compare_tgeant_phast.conf`.
2. Execute:

   ```bash
   bash scripts/command_compare_tgeant_phast.sh
   ```

The macro compares `h_vertex_x/y/z`, `h_muon_p/E/theta/phi` (including fine-angle zoom), `h_q2_nominal_beam_approx` (including zoom), `h_nu_nominal_beam_approx`, and 2D vertex and muon/Q² distributions. It does **not** compare PHAST track counts to raw hit multiplicities or GEANT energy deposition to a reconstructed calorimeter cluster without calibration/reconstruction.

## Important parsing details verified against uploaded example

The uploaded `LO_scaled_07_09_batch2956_1000_run000.tgeant.gz` has **1,000 events** and the following sections: `#EVENT`, `#BEAMDATA`, `#TRACKING`, `#TRIGGER`, `#CALORIMETER`, `#TPCIonization`, `#TPCPadHit`, `#EVENT_END`.

- In `#BEAMDATA`, first three numbers following the marker are vertex coordinates **in mm**, converted to **cm** for matching the current PHAST histogram axes; check that your PHAST version uses cm. A large auxiliary line and numbered trajectory list follow the generator particle list.
- In this **uploaded LO example** the **first** muon entry in the beam particle list is the *transported incident* muon and the **second** muon is the *generated outgoing* muon. The parser explicitly chooses the second one. If your NLO generation or file layout differs, update `truth_muon()` accordingly before processing.
- In detector-hit sections, the 9th whitespace-separated field (index 8) is **deposited energy (`eDep`) in MeV**. The stored hit particle energies/momenta are likewise in MeV/MeV/c; generated beamdata momenta are in GeV/c.
- For the shared muon kinematics, both analyzers use a 70 GeV/c scattered-muon candidate threshold.
- The nominal-beam `Q²` is calculated from generated outgoing muon and a **100 GeV incident muon exactly along +z**. It is diagnostic ONLY, because the recorded incident transported muon may already have lost energy upstream. It is not the physical event-wise generator `Q²` when the nominal-beam assumptions fail. `h_q2_proton_truth` is separately computed from recoil proton energy assuming an initially stationary proton; for radiative events muon-defined and hadron-defined Q² can differ.
- `nuApprox = E_nominal - E_mu_out` is **not** the detector's energy loss. Detector eDep is summed from individual tracking, trigger, calorimeter and TPC ionization hits in **MeV**. Do NOT compare these uncalibrated GEANT hit sums directly with reconstructed calorimeter energy or PHAST `XX0` (material radiation lengths).
- Event weights: The numeric fields on `#EVENT` are **not assumed to be statistical MC weights** without verifying their meaning in the TGEANT exporter; one event contributes one count.
- The waveforms from `#TPCPadHit` are **ADC amplitudes**, not MeV. Do not convert their sum directly to deposited energy.
- Matching names and axes is not sufficient to ensure equivalent observables. In particular a reconstructed primary vertex only exists for successfully reconstructed events; MC-truth vertices exist for all generated events. Detector/trigger acceptance, MC weights, beam directions, energy changes, resolution and tracking/vertex selection can change the distributions.

### Test status

The Python gzip parser **was executed on the full uploaded sample**, correctly completed **1,000 events**, and counted all sections. The summary is saved in `sample_validation.txt`. ROOT/PyROOT and your PHAST libraries are not installed in this execution environment, so the actual ROOT writing, ROOT macros and the PHAST UserEvent **could not be compiled or executed here**. ROOT/PyROOT must be tested on Themis before full-scale runs.
