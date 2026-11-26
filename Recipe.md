# Recipe

## short workflow
* Generate 4-Vectors with McMule: (19 h)
  
* Check if seed is good (chi^2)

* Analyze 4vecs: (10 s)
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_mcmule4vecs.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_mcmule4vecs.sh 
---
* Generate Events with McMule: (1h 15min ~14 M events)
  * copy all .mcmule from /out folder to /input folder
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_evtgen_mcmule.conf

  * Seeds: automatisch ueber inputs/used_seeds_evtgen.txt ("auto" = naechster freier Seed; benutzte Seeds werden abgelehnt)
  * Output landet in: $basefolder/$anafolder/evtgen_<MODE>_<N>x<nenter>_seed<seed>/part_*

  * cd /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/mcmule_eventgen
  * bash submit_nlo_jobs.sh CU 5000 20 auto 10   (nur CU; N wird fuer CU NICHT gebraucht)
    * alter Ablauf N + CU: bash submit_nlo_jobs.sh NCU 5000 20 auto 10 (Dieser Fall dauert: 45 min + 30 min; 50 GB)

  * output anschauen:
    * head -n 200 yyy.lhe 
    * tail -n 200 yyy.lhe 
  
* Parse .lhe to .root: (< 5 min)
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_mcmule_lhe_to_root.conf
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/mcmule_eventgen/command_mcmule_lhe_to_root.sh
  
* Analyze: (2 min)
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_root_NLO_batch.conf
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/mcmule_eventgen/command_analyze_mcmule_evtgen.sh

* Split in batches: (3.5 min)
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_split_root.conf
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_split_root.sh

* Find first Photon Energy bin:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_find_first_mcmule_photon_bin.conf
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_find_first_mcmule_photon_bin.sh
---
* Evtgen with AMBER: ( 6h pro 700,000 events; 30min pro batch)
  * Change: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_tgeant_amber.conf
  * bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/submit_tgeant.sh
  * check run: qstat

* Warm-up Cache: time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_warmup_cache.sh 

* Analyze TGEANT output: (7.5 h)
    * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_tgeant_amber_batches.conf
    * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_tgeant_amber_batches.sh

..................... To Do .....................

* Compare McMule & AMBER TGEANT:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_mcmule_amber.conf 
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_compare_mcmule_amber.sh 

.................................................

* Coral+Phast: (7000 batches in 45 min)
  * Histogram definitions: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/phast/user/u300_PRM_mu_p.cc
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_coral_phast.conf 
  * bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/submit_coral_phast.sh
  * monitor: https://mon.e18.ph.tum.de/cmk/check_mk/dashboard.py?name=cluster

* Merge output: (12 min)
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_merge_coralphast.conf
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_merge_coralphast_autocount.sh

* Merge specific files:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_merge_coralphast_files.conf 
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_merge_coralphast_files.sh 

* Analyze Coral+Phast:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_coralphast.conf 
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_coralphast.sh 

---

* Compare AMBER Before & After Reconstruction:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_amber_rec.conf 
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_compare_amber_rec.sh 

* Compare McMule & AMBER after Reconstruction:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_amber_rec.conf 
  * time bash /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_compare_amber_rec.sh 


## Generate 4-Vectors with McMule
* If you want to Generate Events afterwards, only do 1 seed, but 2 xi cuts
* Change inputstrings: nano /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputstrings
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_create.sh
* !!! Go into the newly created folder: nano mp2mp.toml
  * Fuer NLO0 und R, use:
    * mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"
  * Fuer other pieces use:
    * mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-devel/build/src/mcmule"
* correct xi cut in .menu file for NLO0 (automatically set to 1.0)
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_mcmule_gen4vecs.conf
* qsub /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/submit_mcmule_gen4vecs.sh

## Q2 Exact vs. Ultrarel. Limit
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_Q2approx_test.conf
* time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_Q2approx_test.sh
    
## Coral step
* Delete old output: grep -a -nE 'Number of Processed Events|events had been written|TTrackAppend|FATAL|double free' \
/tmp/ge93juy/coral_test/coral.log
* run: printf '0\n' | ./coral.exe /tmp/ge93juy/coral_test/test.opt \
  > /tmp/ge93juy/coral_test/coral.log 2>&1
* Check: grep -a -nE 'Number of Processed Events|events had been written|TTrackAppend|FATAL|double free' \
/tmp/ge93juy/coral_test/coral.log
* set number of events directly with:  sed -i 's/^events to read 100/events to read 0/' \
/tmp/ge93juy/coral_test/test.opt

## Good To Know
* AMBER PRM Geometrie (TGEANT: $TGEANT/share/prm_2022/alignment_prm.xml, alle Werte in mm):
  * Target / Vertex:  z = -3200
  * ECAL1:  x = -0.126, y = -3.207, z = 11074
  * ECAL2:  x = 450,    y = -3.968, z = 33252   <- relevant fuer Photonen
  * Abstand Vertex -> ECAL2: 33252 - (-3200) = 36452 mm = 36.452 m
  * Photon-Hit auf ECAL2 (gerade Linie, kein Feld): x = d*tan(th)*cos(phi), y = d*tan(th)*sin(phi), d = 36.452 m
    relativ zur ECAL2-Mitte: x - 450 mm, y + 3.968 mm
  * Userfile d_detec: bis 25.11.2026 = 30 m, danach 36.452 m.
    Alte 4vec-Laeufe umrechnen: x5, y5 * 36.452/30 = * 1.21507 (Dichte dsigma/dx entsprechend / 1.21507)
* Connect to Windows: xfreerdp /v:TUPHE18-Pool36.e18.ph.tum.de /u:ge93juy /d:ads /dynamic-resolution
* shift + prnt for scrrenshot of an area saving to Pictures
* vim:
  * Nur speichern: :w
  * Speichern und schließen: :wq
  * Ohne Speichern schließen: :q!
  * Aus Insert-Mode raus: ESC
  * Alles Markieren und loeschen: ggVGd
  * Alles Markieren und kopieren: ggVGy
* check storage: df -h
  * specific  loc: v /nfs/freenas/tuph/e18/project/prm/mleibelt/
  * df -h /nfs/momos/user/mleibelt/
* Git:
  * Rebase: git add <datei>
git rebase --continue
git pull --rebase
* LATEX:
  * Bereinigen und neu kompilieren: latexmk -C
* TGEANT: 
  * Did you change .cc,.cpp,.hh code from TGEANT?
    * build newly: cd /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build
    * make -j4
    * $TGEANT/build/bin/TGEANT /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml
  * Beam settings file: /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml
  * Header File:  nano src/simulation/beam/include/T4ROOTFileEvent.hh
* Terminal:
  * nautilus . &
  * kill process in terminal: kill -9 PID
  * reverse search: strg + r
  * start of line: strg + a
  * check on process status: top
  * leave status: q
  * tail worker* (check status worker files; if more lines needed: -n numberoflines)
  * head worker* (analogous to above but beginning)
* ctrl + / (mehrere Zeilen kommentiern)
* SERVER:
  * ssh -p 22222 ge93juy@login.e18.ph.tum.de
  * ssh themis.e18.ph.tum.de -o "VerifyHostKeyDNS=yes"
  * screen -r
  * ctr A+D
* Rename: F2
* Suche: F3
* Tastatur wechseln: Windows+spacebar




## Links:
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/608536905/E18+Alma+Server+Hardware
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864971/E18+Storage+Server
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864924/GridEngine
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864795/home

