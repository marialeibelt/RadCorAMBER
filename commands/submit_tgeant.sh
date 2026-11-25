#!/usr/bin/env bash
# Submit sequential ROOT batches sorted by part and numeric batch index.
# Uses the existing amber_tgeant.job filelist/snapshot contract.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
CONFIG="${TGEANT_CONFIG:-${SCRIPT_DIR}/../inputs/input_tgeant_amber.conf}"
JOBFILE="${TGEANT_JOBFILE:-${SCRIPT_DIR}/amber_tgeant.job}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file missing: $CONFIG" >&2
    exit 1
fi
source "$CONFIG"
for key in INPUTDIR BATCH_SIZE TARGET_EVENTS START_BATCH DRY_RUN OUTPUTNAME TGEANTrunsfolder OUTPUTDIR SETTINGS_TEMPLATE MAX_JOBS QUEUE MEMORY EST_MB_PER_1000; do
    if [[ -z "${!key:-}" ]]; then
        echo "ERROR: Missing configuration value: $key" >&2
        exit 1
    fi
done

for key in BATCH_SIZE START_BATCH MAX_JOBS EST_MB_PER_1000; do
    if ! [[ "${!key}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $key must be a nonnegative integer (got '${!key}')." >&2
        exit 1
    fi
done
if (( BATCH_SIZE < 1 || MAX_JOBS < 1 )); then
    echo "ERROR: BATCH_SIZE and MAX_JOBS must be positive." >&2
    exit 1
fi
if [[ "$TARGET_EVENTS" != "ALL_REMAINING" ]]; then
    if ! [[ "$TARGET_EVENTS" =~ ^[1-9][0-9]*$ ]]; then
        echo "ERROR: TARGET_EVENTS must be a positive integer or ALL_REMAINING." >&2
        exit 1
    fi
    if (( TARGET_EVENTS % BATCH_SIZE != 0 )); then
        echo "ERROR: Numeric TARGET_EVENTS must be a multiple of BATCH_SIZE=$BATCH_SIZE." >&2
        exit 1
    fi
fi
if [[ "$DRY_RUN" != 0 && "$DRY_RUN" != 1 ]]; then
    echo "ERROR: DRY_RUN must be 0 or 1." >&2
    exit 1
fi
[[ -d "$INPUTDIR" ]] || { echo "ERROR: Input folder does not exist: $INPUTDIR" >&2; exit 1; }
[[ -f "$JOBFILE" ]] || { echo "ERROR: Job script does not exist: $JOBFILE" >&2; exit 1; }
[[ -f "$SETTINGS_TEMPLATE" ]] || { echo "ERROR: TGEANT settings XML missing: $SETTINGS_TEMPLATE" >&2; exit 1; }

# The snapshot makes queued/running jobs independent of subsequent config edits.
SNAPSHOT_DIR="${TGEANTrunsfolder}/submission_configs"
mkdir -p "$SNAPSHOT_DIR"
STAMP="$(date +%Y%m%d_%H%M%S)_$$"
CONFIG_SNAPSHOT="$SNAPSHOT_DIR/input_tgeant_amber_${STAMP}.conf"
FILELIST="$SNAPSHOT_DIR/inputfiles_${STAMP}.txt"
MANIFEST="$SNAPSHOT_DIR/inputfiles_${STAMP}.tsv"
cp -- "$CONFIG" "$CONFIG_SNAPSHOT"

# A numeric TARGET_EVENTS selects complete 1000-entry batches, for example
# day 1: START_BATCH=0 TARGET_EVENTS=7000000.
# ALL_REMAINING selects all unselected complete batches starting at START_BATCH,
# followed by ALL partial tail batches. For day 2, START_BATCH=7000.
# This preserves stage-1 selection order and includes every remaining event.
python3 - "$INPUTDIR" "$BATCH_SIZE" "$START_BATCH" "$TARGET_EVENTS" "$FILELIST" "$MANIFEST" <<'PYTHON_SELECT'
import glob
import os
import re
import sys
from collections import Counter

inputdir, batchsize, start, target, outlist, manifest = sys.argv[1:]
batchsize, start = int(batchsize), int(start)
pattern = re.compile(r'^part_(\d+)_.+_batch(\d+)_(\d+)\.root$')
full, partial, ignored = [], [], []
for path in glob.glob(os.path.join(inputdir, 'part_*_*.root')):
    basename = os.path.basename(path)
    match = pattern.fullmatch(basename)
    if not match:
        ignored.append(basename)
        continue
    part, batch, n = map(int, match.groups())
    item = (part, batch, n, path)
    if n == batchsize:
        full.append(item)
    elif 0 < n < batchsize:
        partial.append(item)
    else:
        raise SystemExit(f'ERROR: unexpected batch size {n} in {basename}')
if not full:
    raise SystemExit('ERROR: no full-size ROOT batches found; check INPUTDIR and filenames')
full.sort(key=lambda item: (item[0], item[1], item[3]))
partial.sort(key=lambda item: (item[0], item[1], item[3]))
identifiers = [(part, batch) for part, batch, _, _ in full + partial]
if len(set(identifiers)) != len(identifiers):
    raise SystemExit('ERROR: duplicate (part, batch) numbering. Refusing ambiguous selection.')
if start >= len(full):
    raise SystemExit(f'ERROR: START_BATCH={start} exceeds complete-batch count {len(full)}.')
if target == 'ALL_REMAINING':
    selected = full[start:] + partial
    print('ALL_REMAINING mode: all remaining full batches plus every partial tail file')
else:
    need = int(target) // batchsize
    end = start + need
    if end > len(full):
        raise SystemExit(f'ERROR: requested positions {start}..{end-1}, '
                         f'but there are only {len(full)} complete {batchsize}-event files.')
    selected = full[start:end]
    print(f'Exact event target: selected full-batch positions [{start}:{end})')
if not selected:
    raise SystemExit('ERROR: no files selected')
with open(outlist, 'w', encoding='utf-8') as fh:
    for _, _, _, path in selected:
        fh.write(os.path.abspath(path) + '\n')
with open(manifest, 'w', encoding='utf-8') as fh:
    fh.write('task_id\tpart\tbatch\tevents\tinput_root\n')
    for i, (part, batch, n, path) in enumerate(selected, 1):
        fh.write(f'{i}\t{part}\t{batch}\t{n}\t{os.path.abspath(path)}\n')
counts = Counter(part for part, _, _, _ in selected)
n_events = sum(n for _, _, n, _ in selected)
print(f'Full batches available: {len(full):,}; partial batches available: {len(partial)}')
print(f'Selected: {len(selected):,} ROOT files, {n_events:,} events')
print(f'First file: {os.path.basename(selected[0][3])}')
print(f'Last file:  {os.path.basename(selected[-1][3])}')
print('Files by part:', ', '.join(f'part_{p}: {n}' for p, n in sorted(counts.items())))
if ignored:
    print(f'WARNING: ignored {len(ignored)} files with unexpected names')
PYTHON_SELECT

NFILES=$(wc -l < "$FILELIST")
ACTUAL_EVENTS=$(awk -F '\t' 'NR>1 {sum += $4} END {printf "%.0f", sum}' "$MANIFEST")
if (( NFILES < 1 || ACTUAL_EVENTS < 1 )); then
    echo "ERROR: no events selected." >&2
    exit 1
fi
TC="$MAX_JOBS"
if (( TC > NFILES )); then TC="$NFILES"; fi

# The MB estimate is decimal. 127MB * 7000 = 889GB (~828GiB).
EST_MB=$((ACTUAL_EVENTS * EST_MB_PER_1000 / 1000))
EST_GB=$(awk -v mb="$EST_MB" 'BEGIN {printf "%.1f", mb/1000}')
EST_GIB=$(awk -v mb="$EST_MB" 'BEGIN {printf "%.1f", mb*1000000/1073741824}')

cat <<SUMMARY
============================================================
AMBER TGEANT — NLO-CU sequential selection
============================================================
Input folder:             $INPUTDIR
Output folder:            $OUTPUTDIR
Requested events:         $TARGET_EVENTS
Actual selected events:   $ACTUAL_EVENTS
Selected ROOT files:      $NFILES
Start at full-batch index:$START_BATCH
Concurrent array tasks:  $TC
Queue / memory:           $QUEUE / $MEMORY per task
Estimated new output:     ~${EST_GB} GB (~${EST_GIB} GiB), HIGHLY APPROXIMATE
Snapshot:                 $CONFIG_SNAPSHOT
Selected input file list: $FILELIST
Detailed manifest:        $MANIFEST
DRY_RUN:                  $DRY_RUN
============================================================
SUMMARY

echo "NOTE: The storage estimate comes from an older run and is not a capacity check.
NOTE: Keep the ROOT input set unchanged between stages so the sorted indices remain stable."
echo "NOTE: Check quota (e.g. quota -s) and leave room for existing files/logs."

if [[ "$DRY_RUN" == 1 ]]; then
    echo "PREVIEW ONLY: No SGE jobs submitted. Set DRY_RUN=0 after checking storage."
    exit 0
fi

# SERIALIZE SUBMISSIONS so that two submissions cannot reserve the same
# input ROOT files before the first TGEANT output has been written.
# This protects submissions made with THIS version of the script.
REGISTRY="${TGEANTrunsfolder}/submission_configs/tgeant_reserved_inputs"
mkdir -p "$REGISTRY"
exec 9> "$REGISTRY/.selection.lock"
flock -x 9

# Compare with previously submitted file lists (not previews). Reject ANY
# overlap even if the earlier array is still queued and has produced no output.
python3 - "$REGISTRY" "$FILELIST" <<'PYTHON_CHECK'
from pathlib import Path
import sys
registry, selected = map(Path, sys.argv[1:])
files = set(selected.read_text().splitlines())
for prev in sorted(registry.glob('submitted_*.txt')):
    overlap = files.intersection(prev.read_text().splitlines())
    if overlap:
        examples = "\n  ".join(sorted(overlap)[:3])
        raise SystemExit(f'ERROR: {len(overlap)} ROOT batches already reserved by {prev.name}.\n  {examples}')
print('Overlap check passed: none of these ROOT inputs was submitted previously with this script.')
PYTHON_CHECK

# Guard against accidentally reusing an output run name with completed
# *selected* files (non-overlapping stages remain possible).
if [[ -d "$OUTPUTDIR" ]]; then
    collisions=0
    while IFS= read -r infile; do
        runname="${infile##*/}"
        runname="${runname%.root}"
        shopt -s nullglob
        existing=("${OUTPUTDIR}/${runname}"*_run*.tgeant.gz)
        shopt -u nullglob
        if (( ${#existing[@]} > 0 )); then
            echo "ERROR: Output for selected input already exists: ${existing[0]}" >&2
            collisions=$((collisions+1))
            if (( collisions >= 3 )); then break; fi
        fi
    done < "$FILELIST"
    if (( collisions > 0 )); then
        echo "ERROR: refusing to submit possible duplicate TGEANT runs." >&2
        exit 1
    fi
fi


# SGE stdout/stderr are kept in task logs by amber_tgeant.job.
SGE_LOG_DIR="${TGEANTrunsfolder}/log/sge"
mkdir -p "$SGE_LOG_DIR" "$OUTPUTDIR"
# Copy selection to a temporary reservation before qsub, under the lock.
RESERVATION_TMP="${REGISTRY}/.pending_${STAMP}.txt"
RESERVATION="${REGISTRY}/submitted_${STAMP}.txt"
cp -- "$FILELIST" "$RESERVATION_TMP"
if qsub \
    -N amber_tgeant_nlo \
    -S /bin/bash \
    -wd "$TGEANTrunsfolder" \
    -j y \
    -o "$SGE_LOG_DIR" \
    -t "1-${NFILES}" \
    -tc "$TC" \
    -l "${QUEUE}=TRUE,h_pmem=${MEMORY}" \
    -v "CONFIG_SNAPSHOT=${CONFIG_SNAPSHOT},INPUT_FILELIST=${FILELIST}" \
    "$JOBFILE"; then
    mv -- "$RESERVATION_TMP" "$RESERVATION"
    echo "Submission registered in: $RESERVATION"
else
    rm -f -- "$RESERVATION_TMP"
    echo "ERROR: qsub failed; no files were reserved." >&2
    exit 1
fi
