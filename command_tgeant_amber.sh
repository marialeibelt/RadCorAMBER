#!/bin/bash

set -e


# ============================================================
# Read input configuration
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/input_tgeant_amber.conf"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Configuration file does not exist:"
    echo "  $CONFIG"
    exit 1
fi

source "$CONFIG"


# ============================================================
# Check configuration
# ============================================================

if [ -z "$TGEANTrunsfolder" ]; then
    echo "ERROR: TGEANTrunsfolder is not defined in input_tgeant_amber.conf"
    exit 1
fi

if [ -z "$Eventnumber" ]; then
    echo "ERROR: Eventnumber is not defined in input_tgeant_amber.conf"
    exit 1
fi

if [ -z "$inputfile" ]; then
    echo "ERROR: inputfile is not defined in input_tgeant_amber.conf"
    exit 1
fi

if [ -z "$outputfolder" ]; then
    echo "ERROR: outputfolder is not defined in input_tgeant_amber.conf"
    exit 1
fi

if [ -z "$outputname" ]; then
    echo "ERROR: outputname is not defined in input_tgeant_amber.conf"
    exit 1
fi


# ============================================================
# Files
# ============================================================

TEMPLATE="${TGEANTrunsfolder}/settings_Mary_root.xml"
TEMP_SETTINGS="${TGEANTrunsfolder}/settings_Mary_root_tmp.xml"


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

if [ ! -d "$outputfolder" ]; then
    echo "Output folder does not exist."
    echo "Creating:"
    echo "  $outputfolder"

    mkdir -p "$outputfolder"
fi


# ============================================================
# Create temporary settings file
# ============================================================

cp "$TEMPLATE" "$TEMP_SETTINGS"


# Replace number of events
sed -i \
    "s|<numParticles>.*</numParticles>|<numParticles>${Eventnumber}</numParticles>|" \
    "$TEMP_SETTINGS"


# Replace ROOT input file
sed -i \
    "s|<path>.*</path>|<path>${inputfile}</path>|" \
    "$TEMP_SETTINGS"


# Replace output folder
sed -i \
    "s|<outputPath>.*</outputPath>|<outputPath>${outputfolder}</outputPath>|" \
    "$TEMP_SETTINGS"


# Replace output name
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
echo "Output name:       ${outputname}_runXXX.tgeant.gz"
echo "Settings file:     $TEMP_SETTINGS"
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

"$TGEANT/bin/TGEANT" "$TEMP_SETTINGS"
