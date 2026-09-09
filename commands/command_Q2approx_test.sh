#!/bin/bash

# Read input parameters
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_Q2approx_test.conf

#Environment
module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

mkdir -p "$outfolder/$figuredir"

# Compile
g++ "$codefolder/Q2approx_test.cpp" \
    $(root-config --cflags --libs) \
    -o "$codefolder/Q2approx_test"

# Run
"$codefolder/Q2approx_test" \
    "$inputfile_LO" \
    "$inputfile_NLO" \
    "$outfolder/$figuredir"