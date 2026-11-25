#!/bin/bash

set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${SCRIPT_DIR}/../inputs/input_coral_phast.conf"
JOBFILE="${SCRIPT_DIR}/coral_phast.job"

if [ ! -f "$CONFIG" ]; then
    echo "ERROR: Config file not found: $CONFIG"
    exit 1
fi

source "$CONFIG"

for v in INPUTDIR FILE_GLOB CORALPHAST_RUNS OUTPUTDIR MDSTDIR HISTDIR OPTDIR FINALDIR LOGDIR \
         CORALDIR PHASTDIR TGEANTDIR PHASTCORAL_EXE PHAST_EXE PHAST_USER \
         LCG_VIEW TGEANT_SETUP COMPASS_FILES \
         PRM_PRODUCTION_DIR CORALOPTORIGINAL TRAFDICOPT DETECTORSDAT DETECTORSROOT \
         KEEP_MDST SAVE_MDST_ON_PHAST_FAILURE MAX_JOBS QUEUE MEMORY; do
    if [ -z "${!v:-}" ]; then
        echo "ERROR: $v not defined in config"
        exit 1
    fi
done

if [ ! -f "$JOBFILE" ]; then
    echo "ERROR: Job file not found: $JOBFILE"
    exit 1
fi

if [ ! -d "$INPUTDIR" ]; then
    echo "ERROR: Input directory does not exist: $INPUTDIR"
    exit 1
fi

for f in "$CORALOPTORIGINAL" "$TRAFDICOPT" "$DETECTORSDAT" "$DETECTORSROOT" \
         "$LCG_VIEW" "$TGEANT_SETUP"; do
    if [ ! -f "$f" ]; then
        echo "ERROR: Required production file does not exist: $f"
        exit 1
    fi
done

if [ ! -x "$PHASTCORAL_EXE" ]; then
    echo "ERROR: Tested CORAL executable missing or not executable:"
    echo "  $PHASTCORAL_EXE"
    exit 1
fi

if [ ! -x "$PHAST_EXE" ]; then
    echo "ERROR: Standalone PHAST executable missing or not executable:"
    echo "  $PHAST_EXE"
    exit 1
fi

if ! [[ "$PHAST_USER" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: PHAST_USER must be a positive integer (current: $PHAST_USER)"
    exit 1
fi

for flag in KEEP_MDST SAVE_MDST_ON_PHAST_FAILURE; do
    if ! [[ "${!flag}" =~ ^[01]$ ]]; then
        echo "ERROR: $flag must be 0 or 1 (current: ${!flag})"
        exit 1
    fi
done

if ! [[ "$MAX_JOBS" =~ ^[1-9][0-9]*$ ]]; then
    echo "ERROR: MAX_JOBS must be a positive integer (current: $MAX_JOBS)"
    exit 1
fi

mapfile -t INPUTFILES < <(
    find "$INPUTDIR" -maxdepth 1 -type f -name "$FILE_GLOB" -print | sort -V
)

NFILES="${#INPUTFILES[@]}"
if [ "$NFILES" -eq 0 ]; then
    echo "ERROR: No matching TGEANT files found."
    echo "Directory: $INPUTDIR"
    echo "Pattern:   $FILE_GLOB"
    exit 1
fi

DRY_RUN="${DRY_RUN:-1}"
EXPECTED_FILES="${EXPECTED_FILES:-0}"
for x in DRY_RUN EXPECTED_FILES; do
    if ! [[ "${!x}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $x must be a nonnegative integer"
        exit 1
    fi
done
if [[ "$DRY_RUN" != 0 && "$DRY_RUN" != 1 ]]; then
    echo "ERROR: DRY_RUN must be 0 or 1"
    exit 1
fi

MAX_FILES="${MAX_FILES:-0}"
if ! [[ "$MAX_FILES" =~ ^[0-9]+$ ]]; then
    echo "ERROR: MAX_FILES must be 0 or a positive integer (current: $MAX_FILES)"
    exit 1
fi
if [ "$MAX_FILES" -eq 0 ] && [ "$EXPECTED_FILES" -gt 0 ] && [ "$NFILES" -ne "$EXPECTED_FILES" ]; then
    echo "ERROR: Expected $EXPECTED_FILES TGEANT files for the complete run, found $NFILES"
    echo "Verify all TGEANT jobs and the FILE_GLOB before submitting."
    exit 1
fi

if [ "$MAX_FILES" -gt 0 ] && [ "$MAX_FILES" -lt "$NFILES" ]; then
    INPUTFILES=("${INPUTFILES[@]:0:$MAX_FILES}")
    NFILES="${#INPUTFILES[@]}"
fi

TC="$MAX_JOBS"
if [ "$TC" -gt "$NFILES" ]; then
    TC="$NFILES"
fi

mkdir -p "$MDSTDIR" "$HISTDIR" "$OPTDIR" "$FINALDIR" "$LOGDIR"

SNAPSHOT_DIR="${CORALPHAST_RUNS}/submission_configs"
mkdir -p "$SNAPSHOT_DIR"

STAMP="$(date +%Y%m%d_%H%M%S)"
CONFIG_SNAPSHOT="${SNAPSHOT_DIR}/input_coral_phast_${STAMP}_$$.conf"
INPUT_MANIFEST="${SNAPSHOT_DIR}/input_manifest_${STAMP}_$$.txt"

cp "$CONFIG" "$CONFIG_SNAPSHOT"
printf '%s\n' "${INPUTFILES[@]}" > "$INPUT_MANIFEST"

# Sanity-check the PRM settings we intentionally rely on.
if grep -Eq '^[[:space:]]*TraF[[:space:]]+Det2Ignore.*(^|[[:space:]])(AT|FT)([[:space:]]|$)' "$TRAFDICOPT"; then
    echo "ERROR: AT or FT is still present in active TraF Det2Ignore in:"
    echo "  $TRAFDICOPT"
    exit 1
fi

if ! grep -Eq '^[[:space:]]*define[[:space:]]+zone[[:space:]]+-7700[[:space:]]+-1700' "$CORALOPTORIGINAL"; then
    echo "ERROR: Expected upstream PRM zone (-7700 -1700) not found in base option file."
    exit 1
fi

cat <<SUMMARY
============================================================
Submit CORAL -> standalone PHAST UserEvent${PHAST_USER}
============================================================
Input directory:       $INPUTDIR
File pattern:          $FILE_GLOB
TGEANT files:          $NFILES
Array tasks:           1-$NFILES
Max simultaneous:      $TC
First input:           ${INPUTFILES[0]}
Last input:            ${INPUTFILES[$((NFILES-1))]}
Dry run (no qsub):     $DRY_RUN
Queue/class:           $QUEUE
Memory per task:       $MEMORY
CORAL hist output:     $HISTDIR
Option snapshots:      $OPTDIR
Final PHAST output:    $FINALDIR
mDST working/recovery: $MDSTDIR
Keep successful mDST:  $KEEP_MDST
Save mDST on failure:  $SAVE_MDST_ON_PHAST_FAILURE
Logs:                  $LOGDIR
Config snapshot:       $CONFIG_SNAPSHOT
Input manifest:        $INPUT_MANIFEST
CORAL executable:      $PHASTCORAL_EXE
Standalone PHAST:      $PHAST_EXE -u${PHAST_USER}
============================================================
SUMMARY

if [ "$DRY_RUN" -eq 1 ]; then
    echo "DRY_RUN=1: input preview completed; NO jobs submitted."
    echo "Set DRY_RUN=0 in the config when ready."
    exit 0
fi

qsub \
    -N coral_phast \
    -S /bin/bash \
    -wd "$CORALPHAST_RUNS" \
    -j y \
    -o /dev/null \
    -t "1-${NFILES}" \
    -tc "$TC" \
    -q "$SGE_QUEUE" \
    -l "${QUEUE}_queue=TRUE,h_pmem=${MEMORY}" \
    -v "CONFIG_SNAPSHOT=${CONFIG_SNAPSHOT},INPUT_MANIFEST=${INPUT_MANIFEST}" \
    "$JOBFILE"
