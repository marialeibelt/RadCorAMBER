# Recipe

## To Run

* bash create_command 
* Go into the newly created folder
* nano mp2mp.toml
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-release/mcmule"
* mcmule = "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"
* . ../compileandrun_command

## Run with Python

* python3 radcor.py

## To Evtgen
* create NLO0,R,0 normally
* go in newly created folder
* give gen-mcmule-path in .toml file
* correct xi cut in .menu file for NLO0
* compile with: . ../compileandrun_gen_command
* create folder called "input"
* copy all .mcmule from /out folder to /input folder
* gen.sh:
  * change path to newly created folder
  * change names NLO,0,R .mcmule files
  * check correct xi cut value
  * check that npieces is correct
* rungen.sh:
  * put correct xi values
  * change random seeds
* run: chmod +x ../eventgen/rungen.sh
       nohup ../eventgen/rungen.sh &

* after source code change you need to rebuild:
  * go into mcmule-event-generator folder 
  * ninja -C build

* output anschauen:
  * head -n 200 yyy.lhe 
  * tail -n 200 yyy.lhe 

## .lhe parsing
* python3 /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/lhe_to_root.py xxx/out/yyy.lhe -o xxx/out/xxx.root
* python3 /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/lhe_to_root.py 05_08_evtgen/out/gen-C-1-01-67187.lhe -o 05_08_evtgen/out/05_08_evtgen.root

## Analyze ROOT
* load ROOT: module load root/6.28.06_py3.9.18_cxx17
* g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/analyze_root.cpp $(root-config --cflags --libs) -o analyze_root
* ./analyze_root

### Bsp:
* Compile: g++ /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/lhe_to_root.cpp $(root-config --cflags --libs) -o lhe_to_root
* Run: ./lhe_to_root 05_08_evtgen/out/gen-C-1-01-67187.lhe 05_08_evtgen/out/05_08_evtgen.root
* ./analyze_root

## LATEX
* Bereinigen und neu kompilieren: latexmk -C
  
## To Know
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


## SERVER:
* ssh -XY ge93juy@login.e18.ph.tum.de -p 22222
* ssh -XY hercules

* screen -r
* ctr A+D

## Links:
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/608536905/E18+Alma+Server+Hardware
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864971/E18+Storage+Server
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864924/GridEngine
* https://collab.dvb.bayern/spaces/TUMtuphe18/pages/71864795/home
