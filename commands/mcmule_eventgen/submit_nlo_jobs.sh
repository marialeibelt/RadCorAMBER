#!/usr/bin/env bash
# Submit independent SGE jobs with unique seeds, LHE outputs and logs.
# Usage: bash submit_nlo_jobs.sh [CU|NCU|N] [nenter per job] [itmx] [first seed|auto] [number of jobs]
# Example: SGE_QUEUE=long_epyc bash submit_nlo_jobs.sh CU 5000 20 auto 10
#   CU  : only CU (resampling + unweighting) per job   <- default
#   NCU : weighted N, then CU per job (N is NOT required for CU)
#   N   : only the weighted N stage
#   auto: first seed = highest seed in inputs/used_seeds_evtgen.txt + 1
# Output: $basefolder/$anafolder/<BATCH_LABEL>/part_<i>, i.e. next to the
#         integration given in inputs/input_evtgen_mcmule.conf.
# 10 x 10000 means ten runs with nenter=10000 per run,
# NOT a guarantee of 100000 written LHE events.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT="$SCRIPT_DIR/submit_nlo_sge.sh"
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONFIG="$REPO/inputs/input_evtgen_mcmule.conf"
SEED_REGISTRY="${SEED_REGISTRY:-$REPO/inputs/used_seeds_evtgen.txt}"
[[ -f "$SUBMIT" ]] || { echo "Missing $SUBMIT" >&2; exit 1; }
[[ -f "$SEED_REGISTRY" ]] || { echo "Missing $SEED_REGISTRY" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

MODE="${1:-CU}"
NENTER="${2:-10000}"
ITMX="${3:-20}"
BASE_SEED="${4:-auto}"
N_RUNS="${5:-${N_RUNS:-10}}"

[[ "$MODE" == "CU" || "$MODE" == "NCU" || "$MODE" == "N" ]] || { echo "MODE must be CU, NCU or N" >&2; exit 1; }
if [[ "$BASE_SEED" == "auto" ]]; then
    # Only regular seeds (< 1000000); offset N seeds of NCU runs are ignored.
    last="$(awk '!/^#/ && $1 < 1000000 && $1 > m {m=$1} END {print m+0}' "$SEED_REGISTRY")"
    BASE_SEED=$((last + 1))
fi
for n in "$NENTER" "$ITMX" "$BASE_SEED" "$N_RUNS"; do
    [[ "$n" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid positive integer: $n" >&2; exit 1; }
done

# Check ALL seeds first, so that we never submit only half of a batch.
for ((i=0; i<N_RUNS; i++)); do
    s=$((BASE_SEED + i))
    for t in "$s" $([[ "$MODE" == "NCU" ]] && echo $((s + 1000000))); do
        if awk -v s="$t" '!/^#/ && $1 == s {f=1} END {exit !f}' "$SEED_REGISTRY"; then
            echo "ERROR: Seed $t already used (see $SEED_REGISTRY). Nothing submitted." >&2
            exit 1
        fi
    done
done

BATCH_LABEL="${BATCH_LABEL:-evtgen_${MODE}_${N_RUNS}x${NENTER}_seed${BASE_SEED}}"
case "$BATCH_LABEL" in
    /*|*..*) echo "Invalid BATCH_LABEL (must be a relative path without '..')" >&2; exit 1 ;;
esac

echo "Integration folder: $basefolder/$anafolder (xi=$xi)"
echo "Seeds: $BASE_SEED .. $((BASE_SEED + N_RUNS - 1))"
for ((i=1; i<=N_RUNS; i++)); do
    seed=$((BASE_SEED + i - 1))
    folder="${anafolder}/${BATCH_LABEL}/part_${i}"
    echo "========================================================"
    echo "Submitting part $i/$N_RUNS: nenter=$NENTER itmx=$ITMX seed=$seed"
    ANA_FOLDER="$folder" bash "$SUBMIT" "$MODE" "$NENTER" "$ITMX" "$seed"
done

echo "========================================================"
echo "All $N_RUNS submissions requested. Check qstat -u \"\$USER\"."
echo "Output: $basefolder/$anafolder/$BATCH_LABEL/part_*"
echo "NOTE: Output event counts may vary; check each LHE and preserve its weight normalization."
