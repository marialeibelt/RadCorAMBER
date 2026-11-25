#!/usr/bin/env bash
# Submit independent SGE jobs; each runs N then CU (sequentially), with unique seeds,
# LHE outputs and persistent stdout/stderr logs.
# Usage: bash submit_nlo_jobs.sh [CU|N] [nenter per job] [itmx] [first CU seed] [number of jobs]
# Example: SGE_QUEUE=long_epyc bash submit_nlo_jobs.sh CU 5000 20 10000 10
# CU: N -> wait -> CU, in each job. N: only the weighted N stage.
# 10 x 10000 means ten runs with nenter=10000 per run,
# NOT a guarantee of 100000 written LHE events.
set -euo pipefail

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
SUBMIT="$SCRIPT_DIR/submit_nlo_sge.sh"
[[ -f "$SUBMIT" ]] || { echo "Missing $SUBMIT" >&2; exit 1; }

MODE="${1:-CU}"
NENTER="${2:-10000}"
ITMX="${3:-20}"
BASE_SEED="${4:-35421}"
N_RUNS="${5:-${N_RUNS:-10}}"

[[ "$MODE" == "N" || "$MODE" == "CU" ]] || { echo "MODE must be N or CU" >&2; exit 1; }
for n in "$NENTER" "$ITMX" "$BASE_SEED" "$N_RUNS"; do
    [[ "$n" =~ ^[1-9][0-9]*$ ]] || { echo "Invalid positive integer: $n" >&2; exit 1; }
done

# Keep all jobs in separate folders; the parent is unique by user label.
BATCH_LABEL="${BATCH_LABEL:-evtgen_nlo_${MODE}_NthenCU_${N_RUNS}x${NENTER}_seed${BASE_SEED}}"
case "$BATCH_LABEL" in
    /*|*..*) echo "Invalid BATCH_LABEL (must be a relative path without '..')" >&2; exit 1 ;;
esac

for ((i=1; i<=N_RUNS; i++)); do
    seed=$((BASE_SEED + i - 1))
    folder="evtgen_08_11/${BATCH_LABEL}/part_${i}"
    echo "========================================================"
    echo "Submitting part $i/$N_RUNS: nenter=$NENTER itmx=$ITMX seed=$seed"
    ANA_FOLDER="$folder" bash "$SUBMIT" "$MODE" "$NENTER" "$ITMX" "$seed"
done

echo "========================================================"
echo "All $N_RUNS submissions requested. Check qstat -u \"\$USER\"."
echo "Each part has its own log/gen-...log and log/sge/ scheduler log."
echo "NOTE: Output event counts may vary; check each LHE and preserve its weight normalization."
