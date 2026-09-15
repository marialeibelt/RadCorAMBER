#!/bin/bash

set -euo pipefail


# ============================================================
# FILE LOCATIONS
# ============================================================

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

CONFIG="${SCRIPT_DIR}/../inputs/input_tgeant_amber.conf"
JOBFILE="${SCRIPT_DIR}/amber_tgeant.job"


# ============================================================
# READ CONFIG
# ============================================================

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found:"
    echo "  $CONFIG"
    exit 1
fi

source "$CONFIG"

: "${INPUTDIR:?INPUTDIR not defined in config}"
: "${FILE_GLOB:?FILE_GLOB not defined in config}"
: "${TGEANTrunsfolder:?TGEANTrunsfolder not defined in config}"
: "${OUTPUTDIR:?OUTPUTDIR not defined in config}"
: "${SETTINGS_TEMPLATE:?SETTINGS_TEMPLATE not defined in config}"
: "${MAX_JOBS:?MAX_JOBS not defined in config}"
: "${QUEUE:?QUEUE not defined in config}"
: "${MEMORY:?MEMORY not defined in config}"

if [ ! -f "$JOBFILE" ]; then
    echo "ERROR: Job file not found:"
    echo "  $JOBFILE"
    exit 1
fi

if [ ! -d "$INPUTDIR" ]; then
    echo "ERROR: Input directory does not exist:"
    echo "  $INPUTDIR"
    exit 1
fi

if [ ! -f "$SETTINGS_TEMPLATE" ]; then
    echo "ERROR: Settings template does not exist:"
    echo "  $SETTINGS_TEMPLATE"
    exit 1
fi

if ! [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_JOBS must be a positive integer."
    echo "Current value: $MAX_JOBS"
    exit 1
fi


# ============================================================
# COUNT ROOT BATCH FILES
# ============================================================

mapfile -t INPUTFILES < <(
    find "$INPUTDIR" -maxdepth 1 -type f -name "$FILE_GLOB" -print | sort -V
)

NFILES="${#INPUTFILES[@]}"

if [ "$NFILES" -eq 0 ]; then
    echo "ERROR: No matching ROOT files found."
    echo "Directory: $INPUTDIR"
    echo "Pattern:   $FILE_GLOB"
    exit 1
fi

# No reason to allow more simultaneous jobs than files.
TC="$MAX_JOBS"
if [ "$TC" -gt "$NFILES" ]; then
    TC="$NFILES"
fi


# ============================================================
# SAVE CONFIG SNAPSHOT FOR THIS SUBMISSION
# ============================================================
#
# Every array job reads this frozen copy instead of the live config.
# You can edit input_tgeant_amber.conf for a later run without
# changing jobs that are already queued/running.

SNAPSHOT_DIR="${TGEANTrunsfolder}/submission_configs"
mkdir -p "$SNAPSHOT_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
CONFIG_SNAPSHOT="${SNAPSHOT_DIR}/input_tgeant_amber_${STAMP}_$$.conf"

cp "$CONFIG" "$CONFIG_SNAPSHOT"


# ============================================================
# SUBMISSION SUMMARY
# ============================================================

echo "============================================================"
echo "Submit AMBER TGEANT"
echo "============================================================"
echo "Input directory:   $INPUTDIR"
echo "File pattern:      $FILE_GLOB"
echo "ROOT files found:  $NFILES"
echo "Array tasks:       1-$NFILES"
echo "Max simultaneous:  $TC"
echo "Queue/class:       $QUEUE"
echo "Memory per task:   $MEMORY"
echo "Output directory:  $OUTPUTDIR"
echo "Settings template: $SETTINGS_TEMPLATE"
echo "Config snapshot:   $CONFIG_SNAPSHOT"
echo "============================================================"


# ============================================================
# SUBMIT
# ============================================================
#
# -t 1-NFILES:
#   exactly one array task per existing ROOT file
#
# -tc TC:
#   limits how many of those tasks may run simultaneously
#
# Queue/class:
#   controlled ONLY by QUEUE="..." in input_tgeant_amber.conf
#
# Memory:
#   h_pmem is controlled by MEMORY="..." in the config.

qsub \
    -N amber_tgeant \
    -S /bin/bash \
    -wd "$TGEANTrunsfolder" \
    -j y \
    -o /dev/null \
    -t "1-${NFILES}" \
    -tc "$TC" \
    -l "${QUEUE}=TRUE,h_pmem=${MEMORY}" \
    -v "CONFIG_SNAPSHOT=${CONFIG_SNAPSHOT}" \
    "$JOBFILE"
