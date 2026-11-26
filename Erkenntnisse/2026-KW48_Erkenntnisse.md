# Erkenntnisse KW48 (23.–29.11.2026) – McMule-Strahlungskorrekturen AMBER

Stand: 2026-11-26. Pro Woche eine neue Datei in `Erkenntnisse/`; dazugehörige ToDo-Liste: `meetings/2026-KW48-49_ToDo.md`.

## Überblick

Mit korrigiertem NLO0 (A enthalten) liegen 18.3 M McMule-NLO-Events für µp-Streuung bei 100 GeV vor, fertig geparst, analysiert und in 18,309 TGEANT-Batches gesplittet. Ziel der Arbeit: Strahlungskorrekturen für AMBER verlässlich simulieren, insbesondere reelle Photonen ab etwa 200 MeV.

| Schritt | Ergebnis | Ablage |
| --- | --- | --- |
| Integration NLO0 + R (ξ_c = 1) | NLO0 χ² = 0.68, R χ² = 0.91 | `/nfs/momos/user/mleibelt/evtgen_24_11_NLO0corrected` |
| Eventgenerierung CU | 13 Jobs, Seeds 75448–75460 | `…/evtgen_CU_13x5000_seed75448/part_*` |
| Parsing LHE → ROOT | 18,304,498 Events, 1,083 mit Photon, 25 negative verworfen | `part_*/out/24_11_evtgen_xi1_NLO_CU_part_*.root` |
| Analyze | 10 Plots, `event_statistics.txt`, `mcmule_evtgen.root` | `McMule_Results/Figures/26_11_evtgen_NLO0corrected_18M` |
| Split | 18,309 Batches à 1000 Events (Summe geprüft) | `…/evtgen_CU_13x5000_seed75448/NLO_batches_1K` |
| First photon bin | kleinste E_γ im Sample 2.37 MeV | `McMule_Results/Figures/26_11_McMule_NLO0corrected_first_photon_bin` |

Noch laufend (Stand 26.11.): ξ_c = 0.1-Integration (Job 4879860) und AMBER-Kinematik NLO0 + R (Job 4879859).

## McMule-Grundlagen: Pieces, FKS und Divergenzen

McMule braucht keinen physikalischen Cutoff: Die IR-Divergenzen werden mit FKS-Subtraktion lokal entfernt und heben sich analytisch gegen die Schleife weg.

| Piece | Inhalt |
| --- | --- |
| LO (`mp2mp0`) | Born-Wirkungsquerschnitt |
| F, A | virtuelle Korrekturen (Schleifen) inkl. integriertem weichen Gegenterm (Bedeutung von A noch bestätigen) |
| NLO0 | LO + F + A (seit dem Fix im event-generator-Build cda91c8 enthält NLO0 auch A) |
| R | reelle Abstrahlung eines Photons, FKS-subtrahiert |
| Voll-NLO | NLO0 + R |

**Divergenzen.** UV-Divergenzen gibt es nur in Schleifen; sie verschwinden durch Renormierung (QED, On-Shell-Schema). IR-Divergenzen gibt es in der Schleife (weiches virtuelles Photon) und in R (E_γ → 0); beide sind gleich groß mit umgekehrtem Vorzeichen (Bloch-Nordsieck). Kollineare Divergenzen fehlen, weil das Myon eine Masse hat; es bleiben große Logarithmen ln(m_μ).

**ξ und ξ_c.** ξ ist die normierte Photonenergie im Schwerpunktsystem:

    ξ = 2 E_γ / √s ,     ξ_c = xinormcut · (1 − (m_μ + M_p)² / s)

**So funktioniert die Subtraktion** (Code: `SIGMA_1` in `mcmule-event-generator/src/integrands.f95`):

1. Für jeden Phasenraumpunkt mit Photon (ξ) entsteht ein Event mit Gewicht ~ ξ · |M_real|².
2. Ist ξ < ξ_c, entsteht zusätzlich ein Gegen-Event: gleiche Kinematik, aber ξ = 0 (kein Photon), Gewicht = minus der weiche Grenzfall.
3. Zusammen verhalten sich beide wie [f(ξ) − f(0)]/ξ und bleiben für ξ → 0 endlich.
4. Der abgezogene Teil f(0)/ξ wird in d = 4 − 2ε Dimensionen von 0 bis ξ_c analytisch integriert. Seine 1/ε-Pole heben die IR-Pole der Schleife auf; das Ergebnis steckt in NLO0.

ξ_c ist damit ein Trennparameter, kein Cutoff: Photonen unter ξ_c bleiben in R (mit Gegen-Event), und NLO0 + R ist unabhängig von ξ_c. Der `softcut` (oft 10⁻¹⁰) ist nur ein numerischer Schutz gegen Rundungsfehler bei ξ → 0; für `mp2mpR` ist er in diesem Build aus (−1). Im CU-Event-Sample wirkt ξ_c trotzdem fast wie ein Cutoff, weil sich Photon- und Gegen-Events unter ξ_c beim Resampling weitgehend wegheben.

In der Integration (4vecs) landet das Gegen-Event mit negativem Gewicht in den Bins ohne Photon; deshalb hat R negative Bins und Event/Gegen-Event-Paare in Nachbarbins.

## Eventgenerierung: N, CU und die Photonen

Das CU-Sample ist für Myon- und Proton-Observablen gut, als Photonspektrum bei ξ_c = 1 aber vermutlich nicht verlässlich: Nach dem Resampling bleiben etwa 4× weniger Photonen übrig.

**N und CU.** N ist ein gewichtetes Sample. CU = Cell Resampling (C) + Unweighting (U) und liefert Events mit Gewicht ±1. N ist keine Voraussetzung für CU; beide lesen nur die `.mcmule`-Integrationsdateien.

1. Gegen-Events werden zu eigenen Events ohne Photon mit negativem Gewicht (N-Datei: 1,375 negative Events, keins davon mit Photon).
2. Cell Resampling fasst Events mit ähnlicher Kinematik zu Zellen zusammen und verrechnet die Gewichte. Die Distanz steht in `EVENT_DISTANCE` im Userfile und nutzt nur E_μ und θ_μ (Schwellen 5 MeV und 15 µrad).
3. Unweighting macht daraus Events mit ±1. Im 18-M-Sample blieben 25 negative Events, die `lhe_to_root` verwirft.

**Gewichte und Plots.** Weil alle Gewichte +1 sind, zeigen die Analyse-Plots Events pro Bin (y-Achse "Events / bin"). Für dσ/dx: Anzahl · (σ_gesamt / N_Events) / Bin-Breite.

| Sample (ξ_c = 1) | Events | Photon-Events | Photon-Anteil am Gewicht |
| --- | --- | --- | --- |
| N (gewichtet, kein Resampling; evtgen_08_11 part_1) | 68,120,755 | 1,762 | 0.027 % |
| CU (nach Resampling + Unweighting; evtgen_08_11 part_1) | 1,136,391 | 75 | 0.0066 % |
| CU, neues 18-M-Sample | 18,304,498 | 1,083 | 0.0059 % |

**Vermutung.** Weil `EVENT_DISTANCE` das Photon ignoriert, verrechnet das Resampling Photon-Events mit Gegen-Events gleicher Myon-Kinematik. Myon-Observablen bleiben dabei richtig, Photon-Observablen nicht. Dazu passt das Spektrum: Bei weichen Photonen sollte dN/dlog E_γ unter etwa 10 GeV flach sein; im Sample fällt es unter etwa 1 GeV ab. Der Abfall oberhalb von etwa 20 GeV ist dagegen das kinematische Limit.

Zum Vergleich: Ältere ξ_c = 0.1-Dateien hatten etwa 0.4 % Photon-Anteil am Gewicht, rund 60× mehr als bei ξ_c = 1.

## F-Vergleich: event-generator- vs. devel-Branch

Der Branch erklärt die früher gesehene Abweichung von ~10⁻⁶ mb nicht: F_devel und F_evtgen stimmen auf 3·10⁻⁸ relativ überein (1.1 σ). Die Abweichung ist sehr wahrscheinlich eine statistische Fluktuation der LO/NLO0-Integration.

**Achtung α = 1.** Der devel-Build rechnet mit α = 1 (`mcmule-devel/pre/global_def.f95`). F ist O(α³), also muss F_devel mit α³ = 1/137.036³ multipliziert werden; das rohe Verhältnis ist genau 1/α³ = 2.5734·10⁶.

**Welcher Vergleich ist richtig?** Im alten event-generator-Build (c4ee45b) enthielt NLO0 kein A, dort muss NLO0 − LO = F gelten. Im neuen Build (cda91c8, Lauf 24_11) enthält NLO0 auch A, dort muss NLO0 − LO = F + A gelten. Der manuelle Zusatz von A beim alten Lauf (NLO0 + A − LO vs. F + A) ist gleichwertig zum ersten Vergleich.

| Check | Wert [µb] | Fehler [µb] | Pull |
| --- | --- | --- | --- |
| F_devel − F_evtgen | −6.1·10⁻⁹ | 5.8·10⁻⁹ | −1.1 σ |
| (NLO0 − LO) − F, alter Lauf (= 1.7·10⁻⁶ mb) | +1.73·10⁻³ | 6.0·10⁻⁴ | +2.9 σ |
| (NLO0 − LO_rep) − F | −3.7·10⁻⁴ | 6.7·10⁻⁴ | −0.6 σ |
| (NLO0_rep − LO) − F | +8.8·10⁻⁴ | 6.4·10⁻⁴ | +1.4 σ |
| (NLO0_rep − LO_rep) − F | −1.23·10⁻³ | 7.1·10⁻⁴ | −1.7 σ |
| (NLO0 mit A − LO) − (F + A) | +1.05·10⁻³ | 5.9·10⁻⁴ | +1.8 σ |
| (NLO0 mit A − LO_rep) − (F + A) | −1.06·10⁻³ | 6.6·10⁻⁴ | −1.6 σ |
| LO_rep − LO (zwei unabhängige LO-Läufe) | +2.10·10⁻³ | 6.5·10⁻⁴ | +3.3 σ |

Relativ zu LO (2.6·10⁴ µb) sind diese Abweichungen von der Größenordnung 10⁻⁸ bis 10⁻⁷, also im Bereich der Integrationsfehler. Dass zwei LO-Läufe um 3.3 σ differieren, deutet auf leicht unterschätzte VEGAS-Fehler hin.

**Bin für Bin.** F_devel vs. F_evtgen: χ²/ndf = 1.01 (θ_μ) und 0.83 (Q²). NLO0 mit A − LO vs. F + A: alle Pulls im Rahmen.

**Warum ist NLO0 − LO so unruhig, F aber glatt?** F wird direkt integriert: jeder Bin ist eine eigene, kleine Zahl mit kleinem Fehler (relativ ~3·10⁻⁴). NLO0 − LO ist dagegen die Differenz zweier großer Zahlen, die fast gleich sind: Bei 1.0–1.25 mrad ist LO ≈ 3.6·10⁴ µb/mrad, F + A dagegen nur ≈ 560 µb/mrad (1.5 %). Die absoluten Fehler von LO und NLO0 (je ~0.06 % ≈ 20 µb/mrad) bleiben in der Differenz voll erhalten und sind relativ zu F + A schon ~4 %. Oberhalb von ~1.3 mrad wachsen die LO/NLO0-Fehler auf ~2 % (≈ 300 µb/mrad), während F + A auf ≈ 240 µb/mrad fällt – die Differenz ist dort von Rauschen dominiert (±100 % und mehr). Es ist also keine Physik, sondern Auslöschung (Cancellation): eine kleine Korrektur aus zwei großen, unabhängig verrauschten Zahlen.

Das ist reine Fehlerfortpflanzung: σ(NLO0 − LO) = √(σ_NLO0² + σ_LO²), und dieser absolute Fehler wird durch die kleine Zahl F + A geteilt. Misst man die Abweichung in Einheiten ihres eigenen Fehlers (unteres Panel, Pull), ist sie über den ganzen Bereich flach und mit 0 verträglich (fast alle Punkte zwischen −2 und +2). Werte von mehreren 100 % (bei Q² und bei NLO0 − LO vs. F allein sogar bis ~2000 %) sind deshalb kein Widerspruch, sondern heißen nur: Die Differenz ist in diesen Bins nicht genau genug bestimmt. Bei Q² beginnt der verrauschte Bereich bei Q² ≈ 0.017 GeV², das entspricht θ_μ ≈ 1.3 mrad (Q² ≈ E²θ²).

**Lesehilfe für die Vergleichsplots.** Mittleres Panel: relative Abweichung in Prozent, 100 · (X − Y)/Y; das Band ist ±1σ des kombinierten Fehlers, ebenfalls in Prozent. Unteres Panel: Abweichung geteilt durch ihren Fehler, (X − Y)/σ ("Pull"); sind X und Y verträglich, liegen ~68 % der Punkte zwischen −1 und +1 und ~95 % zwischen −2 und +2.

Plots: `/nfs/momos/user/mleibelt/McMule_Results/Figures/26_11_meeting/` (`*_F_devel_vs_F_evtgen_*`, `*_NLO0rep_minus_LOrep_vs_F_*`, `*_NLO0withA_minus_LO_vs_FplusA_*`). Zahlen: `/nfs/momos/user/mleibelt/McMule_Results/Vals/26_11_meeting/26_11_meeting_F_branches_comparison.txt`. Skript: `Code/compare_F_branches.py`.

## Fluktuationen um 1.4 mrad: VEGAS-Sampling

Die Fluktuationen kommen aus dem Integrator, nicht aus der Physik: Bei LO und NLO0 springt der relative Fehler pro θ_μ-Bin irgendwo zwischen 1.3 und 1.6 mrad abrupt um einen Faktor 20–70. Ab dort sampelt das adaptive VEGAS-Gitter deutlich weniger dicht (= weniger Sample-Punkte pro Bin).

| Lauf | Piece | Rel. Fehler 0.8–1.2 mrad | Rel. Fehler 1.6–2.0 mrad | Faktor | Sprung bei |
| --- | --- | --- | --- | --- | --- |
| 18_11_complete (alt) | LO | 2.9·10⁻⁴ | 1.1·10⁻² | 39 | ~1.55 mrad |
| 18_11_complete (alt) | NLO0 | 2.4·10⁻⁴ | 1.8·10⁻² | 74 | ~1.36 mrad |
| 18_11_NLO0R0 | LO | 5.4·10⁻⁴ | 1.6·10⁻² | 29 | ~1.31 mrad |
| 18_11_NLO0R0 | NLO0 | 5.9·10⁻⁴ | 1.0·10⁻³ | 1.7 | erst ~2.0 mrad |
| 24_11_NLO0corrected | NLO0 | 7.2·10⁻⁴ | 1.4·10⁻² | 20 | ~1.3–1.5 mrad |

F und A zeigen keinen solchen Sprung (Faktor 1.3–2). R wächst glatt, hat aber oberhalb von ~1.5 mrad einzelne Bins mit riesigem Fehler: Event/Gegen-Event-Paare mit großen Gewichten in Nachbarbins (z. B. 1.537/1.545 mrad im F+R-Plot).

**Warum die Fluktuation im alten Delta_NLO0-Plot bei 1.4 mrad begann:** Im alten Lauf springt NLO0 bei 1.36 mrad und LO bei 1.55 mrad. NLO0 − LO erbt ab 1.36 mrad die großen NLO0-Fehler; genau dort setzt die Fluktuation ein.

**Mehr Statistik oder Verschiebung?** Alle Läufe haben dieselbe Statistik (ITMX 50 × 80 M). Die Lage des Sprungs hängt davon ab, wie sich das VEGAS-Gitter im jeweiligen Lauf (Seed) anpasst, und verschiebt sich deshalb zwischen 1.3 und 2.0 mrad. Mehr Statistik senkt die Fehler oberhalb des Sprungs, beseitigt die Stufe aber nicht.

**Warum macht VEGAS das?** VEGAS teilt jede Integrationsvariable in eine feste Zahl von Intervallen (in McMule `ndmx = 50`, fest in `mcmule-event-generator/src/vegas_m.f95`) und wirft in jedes Intervall gleich viele Punkte. Nach jeder Iteration verschiebt es die Intervallgrenzen so, dass jedes Intervall ungefähr gleich viel zum **Gesamt**integral beiträgt. Der Wirkungsquerschnitt fällt mit dem Streuwinkel extrem steil ab (~1/Q⁴), deshalb landen fast alle 50 Intervalle bei sehr kleinen Winkeln. Für den ganzen Bereich oberhalb von ~1.3 mrad bleibt nur noch ein breites Intervall übrig, also wenige Punkte pro Histogramm-Bin. Das ist Absicht: VEGAS optimiert den Fehler des integrierten Wirkungsquerschnitts, nicht den Fehler einzelner Bins im Ausläufer. Wo genau die letzte Grenze landet, hängt von der zufälligen Anpassung im jeweiligen Lauf ab.

**Was helfen könnte (zu prüfen):** mehr VEGAS-Gitterintervalle; die Integration in θ_μ-Bereiche aufteilen (z. B. getrennte Läufe unter/über 1.3 mrad per Userfile-Cut); für die AMBER-Kinematik (0.316–2 mrad) entschärft der engere Cut das Problem bereits.

Plot: `/nfs/momos/user/mleibelt/McMule_Results/Figures/26_11_meeting/26_11_meeting_vegas_relerr_th3.png`. Tabelle: `/nfs/momos/user/mleibelt/McMule_Results/Vals/26_11_meeting/26_11_meeting_vegas_relerr_jumps.tsv`. Skript: `Code/plot_vegas_error_jump.py`. F+R-Check (1.3–1.5 mrad): neuer Lauf χ²/ndf = 1.39 (p = 0.11), alter R-Lauf 2.82 (p = 0.004).

## Detektor, Trees und Userfile

Photon-Hits werden geradlinig auf ECAL2 bei 36.452 m projiziert, sowohl im Event-Tree als auch in den 4vec-Histogrammen (statt bisher 30 m).

| Objekt (TGEANT `share/prm_2022/alignment_prm.xml`) | x [mm] | y [mm] | z [mm] |
| --- | --- | --- | --- |
| Target / Vertex | 0 | 0 | −3200 |
| ECAL1 | −0.126 | −3.207 | 11074 |
| ECAL2 | 450 | −3.968 | 33252 |

    x_ECAL2 = d · tan θ_γ · cos φ_γ ,   y_ECAL2 = d · tan θ_γ · sin φ_γ ,   d = 33252 − (−3200) mm = 36.452 m

Relativ zur ECAL2-Mitte: x − 450 mm, y + 3.968 mm. Photonen treffen bei x ≈ 0, die ECAL2-Mitte liegt bei +450 mm; die Lage des Strahllochs ist noch zu prüfen.

**Event-Tree (`Code/lhe_to_root.cpp`).** Tree `Output` (alle Events) hat neu: `eventID, nPhotons, hasPhoton, photonEnergy, photonTheta, photonPhi, photonAngleToBeam, photonAngleToMuon, photonECAL2X/Y, photonECAL2XRel/YRel, muonTheta, Q2mu, Q2p, protonKineticEnergy, protonTheta, energyBalance` (Laborsystem; GeV, GeV², mrad, mm; −9999 = kein Photon). Tree `Photons`: ein Eintrag pro Photon. Validierung: ohne Photon Q²_μ = Q²_p bis 3·10⁻¹⁴ und Energiebilanz = E_γ bis 4·10⁻¹⁴ GeV; mit Photon Q²_μ ≫ Q²_p (Strahlungssignatur).

**Analyse (`Code/analyze_mcmule_evtgen.py`).** Schreibt `mcmule_evtgen.root` mit Histogrammen und Tree `RealPhotonEvents` (Events mit genau einem Photon). Log-x-Plots nutzen log-Bins; die Myon-Energie ist linear. y-Achse "Events / bin", wenn alle Gewichte +1 sind.

**Userfiles (`userfile.f95`, `userfile_AMBER_kinematic.f95`).** d_detec = 36.452 m, x5/y5-Grenzen ±0.5 m, neue Histogramme `x5rel`, `y5rel`, `Qsqp` (Q² aus Proton), `th35` (Winkel Photon–Myon). Gilt für künftige gen4vecs-Läufe; alte Versionen in `userfiles/Archiv/`. In der Integration kann man x5/y5 nicht nachträglich berechnen, weil `.mcmule` nur Histogramme enthält.

## Speicherbedarf für die Produktion

Begriff: **TGEANT** ist die Geant4-Detektorsimulation von AMBER, im Recipe "Evtgen with AMBER". Sie macht aus den McMule-Events simulierte Detektor-Hits (`.tgeant.gz`). Coral rekonstruiert daraus Spuren (mDST), Phast analysiert die mDST.

Gemessen an den bisherigen Läufen (Batches à 1000 Events):

| Stufe | pro Event | 18.3 M | 700 M | 7000 M |
| --- | --- | --- | --- | --- |
| McMule LHE (CU) | ~0.7 kB | 13 GB | 0.5 TB | 5 TB |
| McMule ROOT (Parts) | ~0.1 kB | 1.9 GB | 70 GB | 0.7 TB |
| McMule ROOT (1K-Batches für TGEANT) | ~0.1 kB | 1.9 GB | 70 GB | 0.7 TB |
| TGEANT `.tgeant.gz` (127 MB/Batch) | ~127 kB | 2.3 TB | 89 TB | 890 TB |
| Coral mDST (`.root`) | ~154 kB | 2.8 TB | 108 TB | 1.1 PB |
| Coral Hists + Opts | ~0.2 kB | 3.7 GB | 140 GB | 1.4 TB |
| PHAST-Output pro Batch | ~0.11 kB | 2 GB | 80 GB | 0.8 TB |

**Was dauerhaft bleibt, wenn `.tgeant.gz` (und LHE) gelöscht werden:**

- Ohne mDST (nur McMule-ROOT, Coral-Hists, PHAST): ~0.5 kB pro Event → **~0.35 TB für 700 M**, **~3.5 TB für 7000 M**.
- Mit mDST (auch eine `.root`-Datei!): ~155 kB pro Event → ~108 TB für 700 M, ~1.1 PB für 7000 M. Die mDST sind der eigentliche Speicherfresser nach TGEANT; dauerhaft behalten ist nur realistisch, wenn daraus vorher reduzierte PHAST-Trees geschrieben werden.

**Arbeitsspeicher für Chunks in Bearbeitung:** TGEANT + mDST gleichzeitig ≈ 0.28 MB pro Event, also ~280 GB pro 1 M Events. Mit 1.1 TB frei auf momos: Chunks ≤ ~3 M Events, oder Streaming (jede `.tgeant.gz` direkt nach ihrem Coral-Job löschen, jede mDST direkt nach Phast).

**Nur TGEANT, ohne Rekonstruktion:**

| Sample | McMule (LHE + ROOT + Batches) | TGEANT | Summe |
| --- | --- | --- | --- |
| 18.3 M (vorhanden) | ~17 GB (schon belegt) | ~2.3 TB | ~2.3 TB neu |
| 36.6 M (doppelt) | ~34 GB | ~4.7 TB | ~4.7 TB neu |

**McMule-Output aufräumen?** Die LHE-Dateien (13 GB für 18.3 M) kann man nach dem Parsing löschen, weil die ROOT-Trees alle Vierervektoren und Gewichte enthalten. Vorher sollte man den LHE-Header (Wirkungsquerschnitt, Generator-Einstellungen) sichern. Die 1K-Batches (1.9 GB) können nach TGEANT weg, weil sie nur eine Kopie der Part-Dateien sind. Verglichen mit TGEANT ist das aber wenig: Der McMule-Output ist nicht das Speicherproblem.

Rechenzeit: TGEANT braucht ~30 min pro 1000er-Batch, also ~350,000 CPU-Stunden für 700 M Events – vermutlich der eigentliche Engpass neben dem Speicher.

Auf momos liegen noch ~1.8 TB alte TGEANT-Ausgaben (`TGEANT_runs/output/LO_07_09_xi1` 602 GB, `LO_07_09_xi1_1` 350 GB, …). Werden vorerst noch gebraucht (Maria, 26.11.).

## Stolperfallen und Einheiten

| Falle | Was gilt |
| --- | --- |
| µb vs. mb | 1.7·10⁻³ µb = 1.7·10⁻⁶ mb; relativ zu LO 6.6·10⁻⁸. pymule `mergefks(...) * conv` liefert µb. |
| α im devel-Build | α = 1 intern → F_devel × α³ |
| NLO0 mit/ohne A | alter Build (c4ee45b): ohne A; neuer Build (cda91c8, ab 24.11.): mit A |
| d_detec | 4vec-Läufe vor dem 25.11.2026: 30 m; danach 36.452 m. Alte x5, y5 × 1.21507 (dσ/dx ÷ 1.21507) |
| `eventID` doppelt | in `Output`/`Photons`: Nummer innerhalb einer Datei; in `RealPhotonEvents`: Nummer über alle Dateien (Chain-Eintrag) |
| φ in rad | Grad = rad · 180/π; φ_μ muss flach sein (Plausibilitäts-Check) |
| Events vs. dσ | CU-Plots zeigen Events/Bin, 4vec-Plots dσ/dx in µb |
| NLO0-ξ im `.menu` | pymule setzt für NLO0 automatisch ξ = 1.0; bei ξ_c = 0.1 von Hand korrigieren |
| `module load gcc/13.2.0` | geht nur nach `module load spack`; Claude-Shells lesen die `.bashrc` nicht (über `~/.bash_profile` gelöst) |
| Seeds | Eventgenerierungs-Seeds stehen in `inputs/used_seeds_evtgen.txt` und werden automatisch vergeben (`auto`) |

## Offene Fragen und nächste Schritte

Das wichtigste offene Ziel: reelle Photonen ab etwa 200 MeV physikalisch verlässlich simulieren, mit einem in log E_γ flachen Spektrum unter ~10 GeV. Aktuelle Liste: `meetings/2026-KW48-49_ToDo.md`.

Offene Frage: Steht A für einen bestimmten Schleifenbeitrag (bisher hier als „virtuelle Korrektur“ beschrieben)?
