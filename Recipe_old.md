# Recipe

## short workflow
* Generate 4-Vectors with McMule
  
* Check if seed is good (chi^2)
  
* Analyze 4vecs:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_mcmule4vecs.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_mcmule4vecs.sh 
  
* Generate Events with McMule
  
* Get Cross Sections: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_getcs.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_getcs.sh 
  
* Parse .lhe to .root: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_lhe_to_root.conf
    * LO: time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_lhe_to_root_MCMULE_unweighted_LO.sh
    * NLO: time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_lhe_to_root_MCMULE_unweighted_NLO.sh
  
* Analyze raw Output:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_root.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_root.sh 
  
* Scale with CS: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_scale.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_scale.sh
  
* Analyze:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_root.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_root.sh 

* Split in batches
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_split_root.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_split_root.sh 

* Evtgen with AMBER
  * Change: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_tgeant_amber.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/submit_tgeant.sh
  * check run: qstat
  
* Warm-up Cache: time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_warmup_cache.sh 

* Analyze TGEANT output:
  * Batches:
    * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_tgeant_amber_batches.conf
    * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_tgeant_amber_batches.sh 
  * Einzelne Dateien:
    * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_tgeant_amber.conf
    * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_tgeant_amber.sh 

* Compare McMule & AMBER TGEANT:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_mcmule_amber.conf 
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_compare_mcmule_amber.sh 

* Coral+Phast
  * Histogram definitions: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/phast/user/u300_PRM_mu_p.cc
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_coral_phast.conf 
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/submit_coral_phast.sh
  * monitor: https://mon.e18.ph.tum.de/cmk/check_mk/dashboard.py?name=cluster
* Get # expected events 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_count_phast_expected_events.conf 
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/count_phast_expected_events.sh 
* Merge output
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_merge_coralphast.conf
  * time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_merge_coralphast.sh

## Generate 4-Vectors with McMule
* If you want to Generate Events afterwards, only do 1 seed, but 2 xi cuts
* Change inputstrings: nano /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputstrings
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_create.sh
* !!! Go into the newly created folder
* nano mp2mp.toml

!!! that is old mcmule with alpha =1 !!!
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"

!!! that is old mcmule alpha =1!!!
* correct xi cut in .menu file for NLO0
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_mcmule_gen4vecs.conf
* time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_mcmule_gen4vecs.sh
  
## Generate Events wth McMule
* copy all .mcmule from /out folder to /input folder
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_evtgen_mcmule.conf
* time /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/rungen.sh

!!!! immer neue seeds !!!!

* bash submit_nlo_jobs.sh CU 5000 20 10000 10

!!!! immer neue seeds !!!!


* output anschauen:
  * head -n 200 yyy.lhe 
  * tail -n 200 yyy.lhe 

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

