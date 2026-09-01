#!/bin/bash

# Read input
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_mcmule_gen.conf

# Compile
/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/compilefile_gen

# Go to individual output folder
cd "$OUTPUTFOLDER" || exit 1

# Start McMule
nohup pymule batch shell menu-mp2mp.menu &

# Create input folder
mkdir -p input
