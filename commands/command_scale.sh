#!/bin/bash

set -eo pipefail

# ============================================
# Input
# ============================================

CONF="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_scale.conf"

source "$CONF"

# ============================================
# Logging
# ============================================

LOGDIR="$(dirname "$outfolder")/log"
mkdir -p "$LOGDIR"

TIMESTAMP=$(date +"%Y-%m-%d_%H-%M-%S")
LOGFILE="$LOGDIR/scale_${addindex}_${TIMESTAMP}.log"

# Everything from here goes to terminal AND logfile
exec > >(tee -a "$LOGFILE") 2>&1

echo "============================================"
echo "Scale ROOT"
echo "============================================"
echo "Started:       $(date)"
echo "Hostname:      $(hostname)"
echo "Working dir:   $(pwd)"
echo
echo "Config file:   $CONF"
echo
echo "codefolder:    $codefolder"
echo "folder:        $folder"
echo "outfolder:     $outfolder"
echo "inputfile_LO:  $inputfile_LO"
echo "inputfile_NLO: $inputfile_NLO"
echo "CSfile:        $CSfile"
echo "addindex:      $addindex"
echo
echo "Logfile:       $LOGFILE"
echo "============================================"

# ============================================
# Check inputs
# ============================================

for file in \
    "$folder/$inputfile_LO" \
    "$folder/$inputfile_NLO" \
    "$CSfile" \
    "$codefolder/scale_root.cpp"
do
    if [[ ! -f "$file" ]]; then
        echo "ERROR: File does not exist:"
        echo "       $file"
        exit 1
    fi
done

mkdir -p "$outfolder"

echo
echo "Input files found successfully."
echo

# ============================================
# Environment
# ============================================

echo "Loading environment..."

module load gcc/13.2.0
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh

echo
echo "g++:"
which g++
g++ --version | head -1

echo
echo "root-config:"
which root-config
root-config --version

# ============================================
# Compile
# ============================================

echo
echo "============================================"
echo "Compiling scale_root.cpp"
echo "============================================"

set -x

g++ "$codefolder/scale_root.cpp" \
    $(root-config --cflags --libs) \
    -o "$codefolder/scale_root"

set +x

echo
echo "Compilation successful."

# ============================================
# Run
# ============================================

echo
echo "============================================"
echo "Running scaling"
echo "============================================"

echo "Command:"
echo "$codefolder/scale_root \\"
echo "  $folder/$inputfile_LO \\"
echo "  $folder/$inputfile_NLO \\"
echo "  $CSfile \\"
echo "  $outfolder \\"
echo "  $addindex"
echo

set -x

"$codefolder/scale_root" \
    "$folder/$inputfile_LO" \
    "$folder/$inputfile_NLO" \
    "$CSfile" \
    "$outfolder" \
    "$addindex"

set +x

echo
echo "============================================"
echo "Finished successfully: $(date)"
echo "Logfile: $LOGFILE"
echo "============================================"