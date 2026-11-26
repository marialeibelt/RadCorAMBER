# TODO / Zukunft

Langfristige Liste. Stand: 2026-11-26. Claude behaelt diese Liste im Blick und meldet sich, wenn etwas relevant wird.

## Offen aus Meeting 2026-11-26 (Details: `meetings/2026-11-26_ToDo.md`)
- [ ] **xi_c-Test (Punkt 6):** Integration xi_c = 0.1, Job 4879860, Ordner `evtgen_25_11_xi01_NLO0corrected`.
      Ziel: NLO0 + R bei xi_c = 0.1 vs. 1.0 (`evtgen_24_11_NLO0corrected`) innerhalb der Fehler gleich?
      Danach Eventgenerierung mit xi_c = 0.1 fuer Photon-/ECAL2-Studien.
- [ ] **AMBER-Kinematik (Punkt 7):** Job 4879859, Ordner `evtgen_25_11_AMBER_kinematic_NLO0corrected`.
      Ziel: CS (`commands/command_getcs.sh`) mit A-Fix, Vergleich mit `evtgen_14_11_xi1_AMBER_kinematic`.
- [ ] **Photon-Unterdrueckung durch CU (Punkt 9):** N vs. CU im xi_c = 0.1 Lauf vergleichen; ggf. Photon in
      `EVENT_DISTANCE` (Userfile) aufnehmen; Frage an McMule-Entwickler.
- [ ] **Flaches log-E_gamma-Spektrum (Meeting-Notiz 10):** unter ~10 GeV sollte dN/dlogE flach sein, ist es bei
      xi_c = 1 nicht. Ziel: 200-MeV-Photonen physikalisch verlaesslich. Mit xi_c = 0.1 pruefen, Marco fragen.
- [ ] **Naechster Recipe-Schritt:** TGEANT mit `INPUTDIR=/nfs/momos/user/mleibelt/evtgen_24_11_NLO0corrected/evtgen_CU_13x5000_seed75448/NLO_batches_1K`
      (18,309 Batches, 18.3 M Events) - braucht Marias OK.
- [x] Ordner `26_11_evtgen_NLO0corrected_14M_logbins` entfernt (26.11.).
- [ ] Optional: eigenen Branch `eventID` in `lhe_to_root.cpp` umbenennen (z.B. `entryInFile`), um die doppelte
      Bedeutung mit `eventID` (Chain-Eintrag) in `analyze_mcmule_evtgen.py` zu vermeiden.

## Erledigt (25./26.11.)
- [x] NLO0-Integration `evtgen_24_11_NLO0corrected` (NLO0 chi2 0.68, R chi2 0.91), Eventgenerierung 13 x CU
      -> 18.3 M Events, Parsing, Analyze, Split, First photon bin. Siehe `meetings/2026-11-26_ToDo.md`.

## Grosses Ziel: automatisierte Statistik-Produktion
Idee: viele kleine "Chunks" (z.B. ein paar Mio. Events) laufen jeweils komplett durch alle
Schritte; nur die reduzierten Trees/Histogramme bleiben, grosse Zwischendateien werden geloescht.

- [ ] **Eine Config pro Produktion** statt globale `inputs/*.conf` per vim aendern:
      z.B. `productions/<name>/production.conf`, alle Skripte bekommen den Pfad als Argument.
      Voraussetzung fuer alles Weitere.
- [ ] **Chunk-Pipeline als SGE-Abhaengigkeitskette**: `qsub -hold_jid` verbindet
      evtgen -> lhe_to_root -> analyze -> split -> TGEANT -> Coral/Phast -> merge -> analyze.
      Jeder Schritt schreibt eine `.ok`-Markerdatei (wie schon `COMPLETE_*.ok`), damit ein Chunk
      nach Absturz am richtigen Schritt weitermacht.
- [ ] **Validierung pro Schritt automatisieren** (das, was jetzt von Hand geprueft wird):
      Eventzahlen, Gewichtssumme/Normierung, Dateigroesse > 0, Fehler im Log (FATAL, double free).
      Erst wenn Check ok -> Loeschen erlaubt.
- [ ] **Loeschen nur nach Check + Freigabe**: Schalter `CONFIRM_DELETE=1` (fragt) / `0` (voll
      automatisch, erst wenn der Ablauf vertraut ist).
- [ ] **Speicher-Budget**: vor Start eines Chunks `df` pruefen, nur so viele Chunks parallel wie
      Platz da ist (z.B. max. 300 GB in Arbeit). momos: 3 TB, aktuell 1.1 TB frei.
- [ ] **Merge ueber Chunks**: `hadd` der Ergebnis-Trees/Histogramme; Normierung (Wirkungsquerschnitt,
      Gewichte) pro Chunk mitschreiben, damit das Kombinieren korrekt ist.
- [ ] Buchhaltung: eine Tabelle pro Produktion (Chunk, Seed, Events, Status, Groesse, Datum).

- [ ] lhe_to_root.cpp verwirft Events mit negativem Gewicht (Test: 3 von 1.14 M -> vernachlaessigbar,
      aber bei kleinem xi_c erneut pruefen).
- [ ] ECAL2-Strahlloch: Photonen landen bei x ~ 0, ECAL2-Zentrum liegt bei x = +450 mm
      -> Lage des Strahllochs in TGEANT pruefen, dann Akzeptanz-Flag in den Tree.

## Aufraeumen / Kleinkram
- [ ] Alte Produktionen pruefen und ggf. loeschen: `evtgen_08_11` (103 GB), `evtgen_31_08` (24 GB)
      - erst nachdem klar ist, dass sie (alte NLO0 ohne A) nicht mehr gebraucht werden.
- [ ] `McMule_Gen4vecs.o*` Logs landen im Repo-Root (`-wd` in submit_mcmule_gen4vecs.sh) ->
      eigenen `logs/`-Ordner nehmen.
- [ ] Recipe: "Compare McMule & AMBER after Reconstruction" nutzt dieselbe Config/dasselbe Skript
      wie "Compare AMBER Before & After" -> klaeren/korrigieren.
- [ ] Recipe-Schritt "Compare McMule & AMBER TGEANT" ist noch To-Do.
- [ ] Neue Integration: Schritt "correct xi cut in .menu file for NLO0" automatisieren
      (command_create.sh koennte das direkt setzen).

## Meeting-ToDos
Aktuelle Liste: `meetings/2026-KW48-49_ToDo.md` (vorher `meetings/2026-11-26_ToDo.md`). Erkenntnisse: `Erkenntnisse/2026-KW48_Erkenntnisse.md`
