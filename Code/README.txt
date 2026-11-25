Correction for AMBER McMule event-generator c4ee45b9 (Nov 2026)

Files:
  analyze_mcmule4vecs.py - revised analysis + physical cross section normalization
  plotting.py           - revised plotting utilities and NLO labels

Validated with Python AST/py_compile and a lightweight mocked pymule execution.
Requires the user's McMule+pymule+LaTeX environment for an actual integration replay.

Verified in provided Fortran code:
  src/mcmule_header.f95 sets physical alpha before integration.
  src/integrands.f95 mp2mpNLO0 includes tree + one-loop.
  pre/auto_mat_el.f95 applies include_tree_level.

Cross sections in the script (pymule conv gives microbarn):
  lo    = mergefks(sigma('mp2mp0')) * conv
  full  = mergefks(sigma('mp2mpNLO0'),sigma('mp2mpR')) * conv
  nlo   = full - lo             # signed correction, not a separate sample
  onlyR = mergefks(sigma('mp2mpR')) * conv

IMPORTANT: Do not use the signed delta NLO cross section as a positive sampling fraction
for independent LO and NLO event files. Before applying it to event scaling,
inspect the generation and ROOT scaling scripts and weights.

User-selected plot settings remain: theta_mu K shown up to 3.9 mrad,
E_mu K shown up to 95 GeV, photon comparisons show no K; NLO correction on top.
Photon R-only diagnostics are generated separately.
phi_gamma labels are in radians rather than milliradians.

The mock test does not prove the physical correctness of all reconstruction
cuts, real emission cut handling, error correlations, or event mixing.
