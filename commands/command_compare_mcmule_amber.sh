#!/bin/bash

# Read input parameters
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_compare_mcmule_amber.conf

#Environment
module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

# Compile
g++ "$codefolder/compare_mcmule_amber.cpp" \
    $(root-config --cflags --libs) \
    -o "$codefolder/compare_mcmule_amber"

# Run
"$codefolder/compare_mcmule_amber" \
    "$inputfile_mcmule" \
    "$inputfile_amber" \
    "$outfolder"
