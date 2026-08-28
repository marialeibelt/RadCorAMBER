#!/bin/bash

set -e


# ============================================================
# Read input configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/input_analyze_tgeant_amber.conf"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file does not exist:"
    echo "  $CONFIG"
    exit 1
fi

source "$CONFIG"


# ============================================================
# Paths
# ============================================================

REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"

ANALYZE_SOURCE="${REPO}/Code/analyze_AMBER_eventgen.cpp"
ANALYZE_EXECUTABLE="${REPO}/analyze_AMBER_eventgen"

PLOT_MACRO="${REPO}/Code/plot_analyze_AMBER_eventgen.cpp"

figurefolder="${figurebase}/${figurefoldername}"


# ============================================================
# Check input
# ============================================================

if [ ! -f "$inputfile" ]; then
    echo "ERROR: TGEANT input file does not exist:"
    echo "  $inputfile"
    exit 1
fi

if [ ! -f "$ANALYZE_SOURCE" ]; then
    echo "ERROR: Analyze source file does not exist:"
    echo "  $ANALYZE_SOURCE"
    exit 1
fi

if [ ! -f "$PLOT_MACRO" ]; then
    echo "ERROR: Plot macro does not exist:"
    echo "  $PLOT_MACRO"
    exit 1
fi


# ============================================================
# Create output folders
# ============================================================

mkdir -p "$(dirname "$histogramfile")"

if [ ! -d "$figurefolder" ]; then
    echo "Creating figure folder:"
    echo "  $figurefolder"

    mkdir -p "$figurefolder"
else
    echo "Figure folder already exists:"
    echo "  $figurefolder"
fi


# ============================================================
# Show settings
# ============================================================

echo
echo "=============================================="
echo "Analyze TGEANT"
echo "=============================================="
echo "Input file:        $inputfile"
echo "Histogram output:  $histogramfile"
echo "Figure folder:     $figurefolder"
echo "=============================================="
echo


# ============================================================
# Load ROOT environment
# ============================================================

module load gcc/13.2.0

source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh


# ============================================================
# Compile analyzer
# ============================================================

echo "Compiling analyze_AMBER_eventgen.cpp ..."

g++ "$ANALYZE_SOURCE" \
    $(root-config --cflags --libs) \
    -lz \
    -o "$ANALYZE_EXECUTABLE"

echo "Compilation finished."
echo


# ============================================================
# Run analyzer
# ============================================================

echo "Running TGEANT analysis ..."

"$ANALYZE_EXECUTABLE" \
    "$inputfile" \
    "$histogramfile"

echo
echo "Analysis finished."
echo


# ============================================================
# Plot histograms
# ============================================================

echo "Plotting histograms ..."

root -l -b -q \
    "${PLOT_MACRO}(\"${histogramfile}\",\"${figurefolder}\")"

echo
echo "=============================================="
echo "Finished"
echo "=============================================="
echo "Histogram file:"
echo "  $histogramfile"
echo
echo "Figures:"
echo "  $figurefolder"
echo "=============================================="
