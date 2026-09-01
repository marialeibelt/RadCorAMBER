# Recipe

## short workflow
* Generate 4-Vectors with McMule
  
* Check if seed is good (chi^2)
  
* Analyze 4vecs:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_mcmule4vecs.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_mcmule4vecs.sh 
  
* Generate Events with McMule
  
* Get Cross Sections: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_getcs.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_getcs.sh 
  
* Parse .lhe to .root: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_lhe_to_root.conf
    * LO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_lhe_to_root_MCMULE_unweighted_LO.sh
    * NLO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_lhe_to_root_MCMULE_unweighted_NLO.sh
  
* Analyze raw Output:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_root.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_root.sh 
  
* Scale with CS: 
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_scale.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_scale.sh
  
* Analyze:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_root.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_root.sh 
  
* Evtgen with AMBER
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_tgeant_amber.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_tgeant_amber.sh 
  
* Warm-up Cache: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_warmup_cache.sh 

* Analyze TGEANT output:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_tgeant_amber.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_analyze_tgeant_amber.sh 

* Compare McMule & AMBER TGEANT:
  * Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_mcmule_amber.conf 
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_compare_mcmule_amber.sh 

## Generate 4-Vectors with McMule
* If you want to Generate Events afterwards, only do 1 seed, but 2 xi cuts
* Change inputstrings: nano /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputstrings
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_create.sh
* !!! Go into the newly created folder
* nano mp2mp.toml
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"
* correct xi cut in .menu file for NLO0
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_mcmule_gen4vecs.conf
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_mcmule_gen4vecs.sh
  
## Generate Events wth McMule
* copy all .mcmule from /out folder to /input folder
* Change input: vim /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_evtgen_mcmule.conf
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/rungen.sh

* output anschauen:
  * head -n 200 yyy.lhe 
  * tail -n 200 yyy.lhe 

  
## Good To Know
* vim:
  * Nur speichern: :w
  * Speichern und schließen: :wq
  * Ohne Speichern schließen: :q!
  * Aus Insert-Mode raus: ESC
  * Alles Markieren und loeschen: ggVGd
  * Alles Markieren und kopieren: ggVGy
* check storage: df -h
  * specific  loc: df -h /nfs/freenas/tuph/e18/project/prm/mleibelt/
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
  * ssh -XY ge93juy@login.e18.ph.tum.de -p 22222
  * ssh -XY hercules
  * screen -r
  * ctr A+D


## Links:
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/608536905/E18+Alma+Server+Hardware
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864971/E18+Storage+Server
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864924/GridEngine
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864795/home

