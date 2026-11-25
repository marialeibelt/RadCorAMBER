#!/usr/bin/env bash
# Compile the existing C++ converter once, then parse every completed CU job.
# Install as: AMBER_RadCor/commands/mcmule_eventgen/command_lhe_to_root_MCMULE_unweighted_NLO_batch.sh
# Usage: bash command_lhe_to_root_MCMULE_unweighted_NLO_batch.sh
set -euo pipefail

REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONFIG="${LHE2ROOT_CONFIG:-$REPO/inputs/input_mcmule_lhe_to_root.conf}"
[[ -f "$CONFIG" ]] || { echo "ERROR: Config missing: $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

: "${codefolder:?Missing codefolder}"
: "${basefolder:?Missing basefolder}"
: "${batch_folder:?Missing batch_folder}"
: "${n_runs:?Missing n_runs}"
: "${base_seed:?Missing base_seed}"
: "${xi:?Missing xi}"
: "${output_prefix:?Missing output_prefix}"
[[ "$n_runs" =~ ^[1-9][0-9]*$ && "$base_seed" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: n_runs and base_seed must be positive integers" >&2; exit 1;
}
case "$batch_folder" in
    /*|*..*) echo "ERROR: batch_folder must be a relative path without '..'" >&2; exit 1 ;;
esac

RUN_ROOT="$basefolder/$batch_folder"
[[ -d "$RUN_ROOT" ]] || { echo "ERROR: Batch directory missing: $RUN_ROOT" >&2; exit 1; }

# Same compiler and ROOT environment as the user's old converter script.
# Skip only for offline tests (normal cluster runs must initialize ROOT).
if [[ "${LHE2ROOT_SKIP_ENV:-0}" != 1 ]]; then
    module load gcc/13.2.0

    set +u
    source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
    set -u
fi

SRC="$codefolder/lhe_to_root.cpp"
BIN="${LHE2ROOT_BIN:-$codefolder/lhe_to_root_nlo_batch}"
if [[ -z "${LHE2ROOT_BIN:-}" ]]; then
    [[ -f "$SRC" ]] || { echo "ERROR: Converter source missing: $SRC" >&2; exit 1; }
    echo "Compiling $SRC once ..."
    g++ "$SRC" $(root-config --cflags --libs) -o "$BIN"
fi
[[ -x "$BIN" ]] || { echo "ERROR: Converter executable missing: $BIN" >&2; exit 1; }

converted=0
skipped_not_ready=0
skipped_existing=0
for ((i=1; i<=n_runs; i++)); do
    seed=$((base_seed + i - 1))
    part="$RUN_ROOT/part_$i"
    lhe="$part/out/gen-CU-1-${xi}-${seed}.lhe"
    root="$part/out/${output_prefix}_part_${i}_seed${seed}.root"
    tmp_root="$part/out/.${output_prefix}_part_${i}_seed${seed}.partial.root"
    logdir="$part/log/lhe_to_root"
    logfile="$logdir/parse-CU-1-${xi}-${seed}.log"
    complete="$part/log/COMPLETE_CU_${seed}.ok"
    parsed_marker="$root.parse_complete.ok"

    printf '\n=== part %d/%d (seed %d) ===\n' "$i" "$n_runs" "$seed"
    if [[ ! -f "$complete" ]]; then
        echo "SKIP: generation not confirmed complete: $complete"
        ((skipped_not_ready+=1))
        continue
    fi
    [[ -s "$lhe" ]] || { echo "ERROR: Complete marker exists but CU LHE missing/empty: $lhe" >&2; exit 1; }

    if [[ -e "$root" ]]; then
        if [[ -f "$parsed_marker" ]]; then
            echo "SKIP: already parsed: $root"
            ((skipped_existing+=1))
            continue
        else
            echo "ERROR: ROOT exists without success marker; refusing overwrite: $root" >&2
            exit 1
        fi
    fi
    [[ ! -e "$tmp_root" ]] || { echo "ERROR: stale partial file, inspect/remove first: $tmp_root" >&2; exit 1; }

    mkdir -p "$logdir"
    echo "INPUT:  $lhe"
    echo "OUTPUT: $root"
    echo "LOG:    $logfile"
    # pipefail is set: a nonzero converter exit status is not hidden by tee.
    if "$BIN" "$lhe" "$tmp_root" 2>&1 | tee "$logfile"; then
        if [[ ! -s "$tmp_root" ]]; then
            echo "ERROR: Converter exited 0 but produced no ROOT file: $tmp_root" >&2
            exit 1
        fi
        mv -- "$tmp_root" "$root"
        touch "$parsed_marker"
        echo "DONE: $root"
        ((converted+=1))
    else
        echo "ERROR: Conversion failed, log: $logfile" >&2
        echo "       Partial output left for inspection (if created): $tmp_root" >&2
        exit 1
    fi
done

printf '\nSummary: converted=%d, already_parsed=%d, generator_not_finished=%d, total=%d\n' \
    "$converted" "$skipped_existing" "$skipped_not_ready" "$n_runs"
if [[ "$skipped_not_ready" -gt 0 ]]; then
    echo "Some jobs are still running or failed. Rerun this script later to parse those parts."
fi
