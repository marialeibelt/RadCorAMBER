#!/bin/bash

set -e

# Input file
CONF="{/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_getcs.conf}"

if [[ ! -f "$CONF" ]]; then
    echo "ERROR: Config file not found:"
    echo "$CONF"
    exit 1
fi

source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-environment/bin/activate

# Read input variables
source "$CONF"

# Check required variables
: "${HOMEDIR:?HOMEDIR not defined in config}"
: "${OUTFOLDER:?OUTFOLDER not defined in config}"
: "${SAVENAME:?SAVENAME not defined in config}"

echo "========================================"
echo "Cross section calculation"
echo "========================================"
echo "HOMEDIR:   $HOMEDIR"
echo "OUTFOLDER: $OUTFOLDER"
echo "SAVENAME:  $SAVENAME"
echo "========================================"

python3 "$CODEDIR/getCS.py"
    --homedir "$HOMEDIR" \
    --outfolder "$OUTFOLDER" \
    --savename "$SAVENAME"
