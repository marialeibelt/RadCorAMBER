# TGEANT event generators — overview, PRM usage, and the McMule route

Reference document. Nothing here is an instruction to implement anything.

---

## 1. How generators are wired

Selection is **by string only** — there is no enum and no switch/factory. Each plugin
self-registers via a file-scope `static T4X* x = new T4X();` at the top of its `.cc`, and
`T4PrimaryGenerator` looks up `<beamPlugin>` in the `T4BeamPluginList::beamPlugins` map
([T4PrimaryGenerator.cc:8-20](src/simulation/beam/src/T4PrimaryGenerator.cc#L8-L20)). Adding a
generator therefore touches no existing dispatch code.

Two layers:

- **`T4BeamBackEnd`** ([T4BeamBackend.hh](src/simulation/beam/include/T4BeamBackend.hh)) — makes
  the primary vertex. The default `generateEvent()` / `initialize()` handle beam-file reading,
  pile-up and target extrapolation, so most plugins inherit them unchanged.
- **`T4ProcessBackEnd`** ([T4ProcessBackend.hh](src/simulation/beam/include/T4ProcessBackend.hh))
  — a `G4VDiscreteProcess` attached to the beam particle and fired at the vertex inside the
  target. The only pure virtual is `PostStepDoIt()`; it fills
  `T4EventManager::getInstance()->getBeamData()` and pushes secondaries.

Everything lives in [src/simulation/beam/](src/simulation/beam/), with settings structs in
[T4SGlobals.hh](lib/settings/include/T4SGlobals.hh) and XML I/O in
[T4SSettingsFileXML.cc](lib/settings/src/T4SSettingsFileXML.cc).

## 2. Available generators

### Physics generators (have a `T4ProcessBackEnd`)

| `<beamPlugin>` | id | Physics | Advantages | Disadvantages |
|---|---|---|---|---|
| `Elastic` | 9 | internal elastic ℓp, 1γ exchange | **what PRM uses**; fast, self-contained, form factor + radius are settings | **no internal radiative corrections**; `expansion` form factor is flagged broken in the ctor ([T4ElasticProcess.cc:41](src/simulation/beam/src/T4ElasticProcess.cc#L41)); G_E and G_M share one radius parameter |
| `ROOTFile` | 10 | full events read from a `TTree` | any external generator, arbitrary final-state multiplicity, no recompile to change physics | no event weights; no pile-up; EOF wraps silently |
| `ascii` | 8 | final states from a text file | trivial format, zero dependencies | momenta only, no vertex, no weights, no beam info |
| `LEPTO` | 2 | inclusive DIS from a binary LEPTO file | full DIS kinematics carried into `beamData` | Fortran unformatted, slurped entirely into RAM at startup |
| `HEPGEN` / `EventGen` | 3 / 5 | exclusive µp: DVCS, π⁰→γγ, ρ⁰→ππ, φ→KK | in-process, well integrated | needs `USE_HEPGEN`; exclusive channels only. `EventGen` throws the vertex uniformly in the target instead of tracking the beam in |
| `PYTHIA` / `PYTHIA8` | 1 / 42 | hadronic, Drell-Yan; p/n target mixing | mature, tunable via a card file | `USE_PYTHIA6` / `USE_PYTHIA8`; not applicable to elastic µp |
| `DJANGOH` | 19 | DIS incl. QED radiative corrections | the only in-tree generator with rad. corr. | `USE_DJANGOH`; DIS, not elastic |
| `Primakoff` | 6 | pion polarisability / Primakoff | has a `radcorr` switch | `USE_PRIMAKOFF`; unrelated to PRM |
| `PhaseSpace` | 1619 | `TGenPhaseSpace` decay chains, JSON-driven | flexible topologies, optional TH2-driven m_X/t′ | pure phase space, no dynamics |
| `DummyGen` | 4 | kills the beam, no secondaries | vertex/target-stepping debugging | no physics |

### Guns (no process attached)

`BeamOnly`, `User`, `Cosmics`, `Cgeantino`, `EcalCalib`, `ElectronBeam`, `PhotonBeam`,
`PhotonCone`, `VisualizationMode`.

### Not present

No HepMC reader, no LHE/LHEF reader, no `G4GeneralParticleSource`. Confirmed by grep over `src/`
and `lib/`.

## 3. What PRM uses

All PRM production configs (`share/prm_2018`, `prm_2021`, `prm_2022`, `uts_test_2023`) set:

```xml
<beamParticle>13</beamParticle>
<beamPlugin>Elastic</beamPlugin>
<beamEnergy>100000</beamEnergy>      <!-- 190000 in prm_2018 -->
<useBeamfile>true</useBeamfile>
<useTargetExtrap>true</useTargetExtrap>
<triggerPlugin>PRM2022</triggerPlugin>
```

The three `settings_prm_*PhysicsStudies_*.xml` files are the exception — they use `BeamOnly`
(a plain gun) for TPC dE/dx studies.

`T4ElasticProcess` ([T4ElasticProcess.cc](src/simulation/beam/src/T4ElasticProcess.cc)) produces
µp → µp, one-photon exchange, target at rest, recoil by momentum conservation:

- Q² sampling: `flat` (actually log-flat, `getQ2Flat()` line 397), `exponential`, `inverse`, or
  `crosssection` (accept/reject on the Rosenbluth-style σ, with or without lepton mass).
- Form factors: dipole with `charge_radius = 0.8775 fm`, `G_E0 = 1.0`, `G_M0 = 2.792847356`.
- Kinematics: `photonE = Q²/2M`, `cosθ_lab` from the invariants, uniform φ, `rotateUz` onto the
  incoming beam direction.
- `SetNumberOfSecondaries(2)` — µ′ and p, always.

**Radiative corrections are entirely absent from this path.** Grepping `radiativ|radcorr` finds
only the Primakoff generator's `radcorr` flag and Geant4's `G4MuBremsstrahlung` in the physics
list. So external (material) radiation is simulated during transport, but vertex,
vacuum-polarisation and soft-photon corrections are not modelled at all. That is the gap McMule
would fill.

## 4. Route for McMule events — the `ROOTFile` plugin

TGEANT already has the bridge. `T4ROOTFileBeam` / `T4ROOTFileProcess` / `T4ROOTFileEvent` were
added in Oct 2019 by C. Dreisbach (`d3977d9a`, "Added event reader based on root file") in the
same series as the `T4GasH2Target` work — built for exactly this purpose for PRM. A sample input
is still shipped at `share/prm_2022/generated.root` (produced elsewhere as
`generatorOutput.root`) and is already wired into `share/prm_2022/settings_prm.xml`.

Nothing in the tree mentions ESEPP or McMule by name, so this plugin is the generic, unnamed
"external generator" entry point rather than an ESEPP-specific reader.

### Expected TTree schema

`T4ROOTFileEvent::setTree()` binds by branch name and **fatals on any missing branch**
([T4ROOTFileEvent.cc:27-94](src/simulation/beam/src/T4ROOTFileEvent.cc#L27-L94)):

| branch | type | unit |
|---|---|---|
| `vertexX`, `vertexY`, `vertexZ` | `double` | mm |
| `beamPID` | `int` | PDG |
| `beamEnergy` | `double` | GeV |
| `beamMomentumX/Y/Z` | `double` | GeV |
| `scatteredPID` | `vector<int>*` | PDG |
| `scatteredEnergy` | `vector<double>*` | GeV |
| `scatteredMomentumX/Y/Z` | `vector<double>*` | GeV |

Extra branches are ignored — `generated.root` also carries `scatteredPhi/Theta/Mass`.

```xml
<BeamSettings>
  <beamPlugin>ROOTFile</beamPlugin>
</BeamSettings>
<ROOTFileGeneratorSettings>
  <path>$TGEANT/share/prm_2022/mcmule.root</path>
  <tree>Output</tree>
</ROOTFileGeneratorSettings>
```

Parsed at [T4SSettingsFileXML.cc:1454-1461](lib/settings/src/T4SSettingsFileXML.cc#L1454-L1461),
only when `beamPlugin == "ROOTFile"`.

### Things to know before feeding McMule through it

- **Variable multiplicity works.** `T4ROOTFileProcess::PostStepDoIt()` loops over
  `scatteredPID->size()`, so a radiated photon is simply a third entry. No code change needed for
  µ′ + p + γ events.
- **`vertex*` is the beam *start* position, not the interaction point.**
  `T4ROOTFileBeam::generateEvent()`
  ([T4ROOTFileBeam.cc:36-63](src/simulation/beam/src/T4ROOTFileBeam.cc#L36-L63)) hands it to
  `T4Extrapolate::extrapAndCalcTargetDist()` and shoots the gun from there; the real vertex is
  sampled inside the target afterwards. It retries the next tree entry until extrapolation hits
  the target. PRM conventions: `beamZStart = -8000`, `beamFileZConvention = -3200`.
- **Beam phase space comes from the file, not the beam file.** McMule fires along a fixed axis,
  so unless the converter samples position/slope from the PRM beam file
  (`$BEAMFILES/output_compass_mu100.dat`, format in
  [T4BeamFileBackend.hh](src/simulation/beam/include/T4BeamFileBackend.hh)) and rotates the final
  state onto it, every event has an identical, unphysical beam.
- **No pile-up.** `T4ROOTFileBeam` overrides `generateEvent()` and never calls
  `generatePileUp()`. PRM production runs with `usePileUp=true`, `beamFlux=0.002` — a `ROOTFile`
  run has no beam pile-up in the TPC.
- **EOF wraps silently.** `readNextEvent()`
  ([T4ROOTFileEvent.cc:113-128](src/simulation/beam/src/T4ROOTFileEvent.cc#L113-L128)) resets to
  entry 0 at the end, so asking for more events than the file holds reuses them with no warning.
  Contrast `LEPTO`, which warns and aborts the run
  ([T4LeptoFile.cc:53-62](src/simulation/beam/src/T4LeptoFile.cc#L53-L62)).
- **Minor:** `new T4Extrapolate` is allocated per event at
  [T4ROOTFileBeam.cc:41](src/simulation/beam/src/T4ROOTFileBeam.cc#L41) and never freed.

### Open question — event weights

McMule is fixed-order and produces **weighted events, including negative weights**. The schema
has no weight branch and `T4ROOTFileProcess` implicitly treats every event as weight 1. Options,
listed for a later decision:

- **(a) Unweight in the converter.** Hit-or-miss on |w|; nothing in TGEANT changes. Clean only if
  the negative-weight fraction is negligible after binning, otherwise it biases the result.
- **(b) Carry the weight through.** Add a `weight` branch, read it in `T4ROOTFileEvent`, park it
  in an unused `T4BeamData` field — `aux` (double) or a free `uservar[i]` slot
  ([T4Event.hh:107-120](lib/events/include/T4Event.hh#L107-L120)); neither is touched by the
  elastic or ROOTFile paths. Would need a check that `TGEANT2ROOT` propagates the chosen field.
- **(c) Weight-binned samples.** Split McMule output into several near-unweighted files by weight
  bracket, run each with its own normalisation. No code change, more bookkeeping.

### Alternative that was considered and not chosen

A dedicated `T4McMule` / `T4McMuleProcess` plugin pair modelled on `T4Elastic`, which would
inherit the default `generateEvent()` and thereby keep the real beam file, pile-up, target
extrapolation and H2-only vertex sampling — taking only the µ′/p/γ kinematics from McMule and
rotating them onto each real beam track. More faithful for production, but more new C++.
