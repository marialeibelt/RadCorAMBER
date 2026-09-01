#!/bin/bash

# Read input parameters
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_scale.conf

#Environment
module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

# Compile
g++ "$codefolder/scale_root.cpp" \
    $(root-config --cflags --libs) \
    -o "$codefolder/scale_root"

# Run
"$codefolder/scale_root" \
    "$folder/$inputfile_LO" \
    "$folder/$inputfile_NLO" \
    "$CSfile" \
    "$outfolder" \
    "$addindex"
