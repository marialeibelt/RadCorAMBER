#!/bin/bash

set -e


# ============================================================
# Read input configuration
# ============================================================

CONFIG="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_tgeant_amber.conf"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file does not exist:"
    echo "  $CONFIG"
    exit 1
fi

source "$CONFIG"


# ============================================================
# Check configuration
# ============================================================

: "${TGEANTrunsfolder:?TGEANTrunsfolder is not defined}"
: "${Eventnumber:?Eventnumber is not defined}"
: "${inputfile:?inputfile is not defined}"
: "${outputfolder:?outputfolder is not defined}"
: "${outputname:?outputname is not defined}"


# ============================================================
# Run-specific names
# ============================================================

# Remove .root from output name
RUNNAME="${outputname%.root}"

# Log folder
LOGFOLDER="${TGEANTrunsfolder}/log"
mkdir -p "$LOGFOLDER"

# Run-specific log
LOGFILE="${LOGFOLDER}/${RUNNAME}.log"

# Temporary settings folder
TEMPFOLDER="${TGEANTrunsfolder}/tmp"
mkdir -p "$TEMPFOLDER"

# Unique settings file for this job
RUNID="${JOB_ID:-$$}"
TEMP_SETTINGS="${TEMPFOLDER}/settings_${RUNNAME}_${RUNID}.xml"

# Remove temporary settings file when job ends
trap 'rm -f "$TEMP_SETTINGS"' EXIT


# ============================================================
# Redirect output to run-specific log
# ============================================================

exec > "$LOGFILE" 2>&1


# ============================================================
# Files
# ============================================================

TEMPLATE="${TGEANTrunsfolder}/settings_Mary_root.xml"


# ============================================================
# Check files and folders
# ============================================================

if [ ! -f "$TEMPLATE" ]; then
    echo "ERROR: Settings template does not exist:"
    echo "  $TEMPLATE"
    exit 1
fi

if [ ! -f "$inputfile" ]; then
    echo "ERROR: Input ROOT file does not exist:"
    echo "  $inputfile"
    exit 1
fi

mkdir -p "$outputfolder"


# ============================================================
# Create temporary settings file
# ============================================================

cp "$TEMPLATE" "$TEMP_SETTINGS"

sed -i \
    "s|<numParticles>.*</numParticles>|<numParticles>${Eventnumber}</numParticles>|" \
    "$TEMP_SETTINGS"

sed -i \
    "s|<path>.*</path>|<path>${inputfile}</path>|" \
    "$TEMP_SETTINGS"

sed -i \
    "s|<outputPath>.*</outputPath>|<outputPath>${outputfolder}</outputPath>|" \
    "$TEMP_SETTINGS"

sed -i \
    "s|<runName>.*</runName>|<runName>${outputname}</runName>|" \
    "$TEMP_SETTINGS"


# ============================================================
# Show used settings
# ============================================================

echo
echo "=============================================="
echo "TGEANT run"
echo "=============================================="
echo "Number of events:  $Eventnumber"
echo "Input ROOT file:   $inputfile"
echo "Output folder:     $outputfolder"
echo "Output name:       ${RUNNAME}_runXXX.tgeant.gz"
echo "Log file:          $LOGFILE"
echo "Settings file:     $TEMP_SETTINGS"
echo "SGE Job ID:        ${JOB_ID:-local}"
echo "=============================================="
echo


# ============================================================
# Set up TGEANT environment
# ============================================================

module unload gcc/13.2.0

source /cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/setup.sh

export GEANT4_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64/cmake/Geant4

export CLHEP_DIR=/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib/CLHEP-2.4.6.4

export TGEANT=/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/TGEANT/build

export LD_LIBRARY_PATH=$TGEANT/lib:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib64:/cvmfs/sft.cern.ch/lcg/views/LCG_104a/x86_64-el9-gcc12-opt/lib:$LD_LIBRARY_PATH

source "$TGEANT/thisgeant.sh"


# ============================================================
# Run TGEANT
# ============================================================

/usr/bin/time -v "$TGEANT/bin/TGEANT" "$TEMP_SETTINGS"