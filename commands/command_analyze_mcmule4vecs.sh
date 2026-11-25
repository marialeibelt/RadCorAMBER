#!/bin/bash
set -euo pipefail

source "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_mcmule4vecs.conf"

mkdir -p "$OUTDIR" "$OUTDIR_VALS"

python3 "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/analyze_mcmule4vecs.py" \
    --momosdir "$MOMOSDIR" \
    --outdir "$OUTDIR" \
    --outdir-vals "$OUTDIR_VALS" \
    --run "$RUN" \
    --rundir "$RUNDIR" \
    --savename "$SAVENAME" \
    --nbins "$NBINS" \
    --pieces "${PIECES[@]}"

