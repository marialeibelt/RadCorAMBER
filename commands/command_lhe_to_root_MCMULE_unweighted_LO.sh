#!/bin/bash

# Read input parameters
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_lhe_to_root.conf

# Environment
module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

# Compile
g++ "$codefolder/lhe_to_root.cpp" \
    $(root-config --cflags --libs) \
    -o "$codefolder/lhe_to_root"

# Run LO
"$codefolder/lhe_to_root" \
    "$folder/$inputfile_LO" \
    "$folder/$outputfile_LO"
