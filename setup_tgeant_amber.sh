module unload gcc/13.2.0

source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/setup.sh

export GEANT4_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64/cmake/Geant4
export CLHEP_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib/CLHEP-2.4.6.4
export TGEANT=/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build

export LD_LIBRARY_PATH=/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build/lib:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib:$LD_LIBRARY_PATH

source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build/thisgeant.sh
