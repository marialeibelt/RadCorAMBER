# Recipe

## short workflow
* Evtgen with McMule
* Parse .lhe to .root: 
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_lhe_to_root.conf
    * LO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_lhe_to_root_MCMULE_unweighted_LO.sh
    * NLO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_lhe_to_root_MCMULE_unweighted_NLO.sh
* Analyze raw Output:
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_analyze_root.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_analyze_root.sh 
* Scale with CS: 
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_scale.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_scale.sh
* Analyze:
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_analyze_root.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_analyze_root.sh 
* Evtgen with AMBER
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_tgeant_amber.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_tgeant_amber.sh 
* Warm-up Cache: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_warmup_cache.sh 
* Analyze TGEANT output:
  * Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_analyze_tgeant_amber.conf
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_analyze_tgeant_amber.sh 


## To Evtgen
* create NLO0,R,0 normally:
  * Change inputstrings: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputstrings
  * /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_create.sh
* go in newly created folder
* give gen-mcmule-path in .toml file
* correct xi cut in .menu file for NLO0
* compile with: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_mcmule_gen.sh
* copy all .mcmule from /out folder to /input folder
* gen.sh:
  * change path to newly created folder
  * change names NLO,0,R .mcmule files
  * check correct xi cut value
  * check that npieces is correct
  
* Go in Analysis Folder!
* rungen.sh:
  * put correct xi values
  * change random seeds
* run: chmod +x /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/rungen.sh
       nohup /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/rungen.sh &

* after source code change you need to rebuild:
  * go into mcmule-event-generator folder 
  * ninja -C build

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
* LATEX:
  * Bereinigen und neu kompilieren: latexmk -C
* nautilus . &
* kill process in terminal: kill -9 PID
* reverse search: strg + r
* start of line: strg + a
* check on process status: top
* leave status: q
* tail worker* (check status worker files; if more lines needed: -n numberoflines)
* head worker* (analogous to above but beginning)
* ctrl + / (mehrere Zeilen kommentiern)
* check storage: df -h
* check storage of specific  loc: df -h /nfs/freenas/tuph/e18/project/prm/mleibelt/
* SERVER:
  * ssh -XY ge93juy@login.e18.ph.tum.de -p 22222
  * ssh -XY hercules

  * screen -r
  * ctr A+D
* Links:
  * https://collab.dvb.bayern/spaces/TUMtuphe18/pages/608536905/E18+Alma+Server+Hardware
  * https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864971/E18+Storage+Server
  * https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864924/GridEngine
  * https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864795/home


## General
* source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-environment/bin/activate
* /usr/bin/python3.9 Code/getCS.py
* more eventgen space: cd /nfs/momos/user/mleibelt/
* Beam settings file: /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml
  * nano /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml
* pROBLEM fILE: nano src/simulation/beam/src/T4ROOTFileEvent.cc
* Header File:  nano src/simulation/beam/include/T4ROOTFileEvent.hh


## Git
* Rebase: git add <datei>
git rebase --continue
git pull --rebase


## To Run
* bash create_command 
* Go into the newly created folder
* nano mp2mp.toml
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-release/mcmule"
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"
* . ../compileandrun_command


## Run with Python
* python3 radcor.py


## .lhe parsing
* Themis: x86_64-el9-gcc13-opt
  * module load gcc/13.2.0
  * source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
  * unweighted: 
    * g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/lhe_to_root.cpp $(root-config --cflags --libs) -o lhe_to_root
    * LO: ./lhe_to_root /nfs/momos/user/mleibelt/05_08_evtgen/out/gen-N-0-10-29367.lhe /nfs/momos/user/mleibelt/05_08_evtgen/out/05_08_evtgen_LO_lab.root
    * NLO: ./lhe_to_root /nfs/momos/user/mleibelt/05_08_evtgen/out/gen-C-1-01-67187.lhe /nfs/momos/user/mleibelt/05_08_evtgen/out/05_08_evtgen_NLO_lab.root
  * weighted: 
    * g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/lhe_to_root_weighted.cpp $(root-config --cflags --libs) -o lhe_to_root_weighted
    * /nfs/momos/user/mleibelt/05_08_evtgen_25_08/lhe_to_root_weighted \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/gen-N-0-10-29334.lhe \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_weighted_tgeffe-2_LO.root

<3

                                        !!!NEW!!!
* Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_lhe_to_root.conf
  * LO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_lhe_to_root_MCMULE_unweighted_LO.sh
  * NLO: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_lhe_to_root_MCMULE_unweighted_NLO.sh


## Scale ROOT
* g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/scale_root.cpp $(root-config --cflags --libs) -o /nfs/momos/user/mleibelt/05_08_evtgen_25_08/scale_root
* /nfs/momos/user/mleibelt/05_08_evtgen_25_08/scale_root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_LO_lab.root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_NLO_lab.root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/14_08_05_08_evtgen_cross_sections.txt /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out
* /nfs/momos/user/mleibelt/05_08_evtgen_25_08/scale_root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_LO_lab.root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_NLO_lab.root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/14_08_05_08_evtgen_cross_sections.txt /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out
  * string loFileName  = argv[1];
    string nloFileName = argv[2];
    string csFileName  = argv[3];
    string outputDir   = argv[4];
    string suffix = argv[5];
* /nfs/momos/user/mleibelt/05_08_evtgen_25_08/scale_root \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_LO_lab.root \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_NLO_lab.root \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/14_08_05_08_evtgen_cross_sections.txt \
/nfs/momos/user/mleibelt/05_08_evtgen_25_08/out \
25_08
<3

                                        !!!NEW!!!
* Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_scale.conf
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_scale.sh


## Analyze ROOT
* Load ROOT:
  * module load gcc/13.2.0
  * source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
* g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/analyze_root.cpp $(root-config --cflags --libs) -o analyze_root
* ./analyze_root LO.root NLO.root output_folder
  * /nfs/momos/user/mleibelt/05_08_evtgen_25_08/analyze_root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_weighted_tgeffe-2.root /nfs/momos/user/mleibelt/05_08_evtgen_25_08/out/05_08_evtgen_25_08_NLO_lab.root /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/output_05_08_analysis_25_08

<3

                                        !!!NEW!!!
* Change input: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/input_analyze_root.conf
* /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/command_analyze_root.sh


## Finally in TGEANT
* Did you change .cc,.cpp,.hh code from TGEANT?
  * build newly: cd /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build
  * make -j4
  * $TGEANT/build/bin/TGEANT /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml


* module unload gcc/13.2.0
* source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/setup.sh
* export GEANT4_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64/cmake/Geant4
* export CLHEP_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib/CLHEP-2.4.6.4
* export LD_LIBRARY_PATH=/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build/lib:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib:$LD_LIBRARY_PATH
* export TGEANT=/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build
* source build/thisgeant.sh
* build/bin/TGEANT /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml 

* !!!NEW!!!-Skript:
  * source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/setup_tgeant_amber.sh
  * (/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build/bin/TGEANT /nfs/momos/user/mleibelt/TGEANT_runs/settings_Mary_root.xml)
  * /nfs/momos/user/mleibelt/TGEANT_runs/run_Mary_root.sh \
10 \
/nfs/momos/user/mleibelt/05_08_evtgen/out/LO_scaled_21_08.root


## Analyze TGEANT Output
* Load ROOT:
  * module load gcc/13.2.0
  * source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
* Compile: g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/analyze_AMBER_eventgen.cpp \
$(root-config --cflags --libs) \
-lz \
-o analyze_AMBER_eventgen
* Run: /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/analyze_AMBER_eventgen \
/nfs/momos/user/mleibelt/TGEANT_runs/output/PRM_run006.tgeant.gz \
/nfs/momos/user/mleibelt/TGEANT_runs/output/PRM_run006_histograms.root

* Plot: 
  * (Cashe warmmachen: time root -l -b -q -e 'TCanvas c;')
  * Default input und output folder: root -l -b -q /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/plot_analyze_AMBER_eventgen.cpp
  * root -l -b -q \
'/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/plot_analyze_AMBER_eventgen.cpp("/nfs/momos/user/mleibelt/TGEANT_runs/output/PRM_run007_histograms.root","/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Figures/AMBER_eventgen")'
