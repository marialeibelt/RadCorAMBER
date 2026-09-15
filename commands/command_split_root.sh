#!/bin/bash

set -e

echo "============================================"
echo "Split ROOT file into batches"
echo "============================================"

# --------------------------------------------------
# Environment
# --------------------------------------------------

module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh


CONFIG="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_split_root.conf"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file does not exist:"
    echo "  $CONFIG"
    exit 1
fi

source "$CONFIG"


# --------------------------------------------------
# Code / executable
# --------------------------------------------------

sourcefile="${codefolder}/split_root.cpp"
executable="${codefolder}/split_root"


# --------------------------------------------------
# Print settings
# --------------------------------------------------

echo "Input file:    $inputfile"
echo "Output folder: $outputfolder"
echo "Batch size:    $batchsize"
echo "Source file:   $sourcefile"
echo ""


# --------------------------------------------------
# Create output folder
# --------------------------------------------------

mkdir -p "$outputfolder"


# --------------------------------------------------
# Compile
# --------------------------------------------------

echo "Compiling split_root.cpp..."

g++ "$sourcefile" $(root-config --cflags --libs) -o "$executable"

echo "Compilation finished."
echo ""


# --------------------------------------------------
# Run
# --------------------------------------------------

echo "Splitting ROOT file..."
echo ""

"$executable" \
    "$inputfile" \
    "$outputfolder" \
    "$batchsize"

echo ""
echo "============================================"
echo "Finished"
echo "============================================"