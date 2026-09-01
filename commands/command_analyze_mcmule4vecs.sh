#!/bin/bash

source "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_analyze_mcmule4vecs.config"

python3 "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/Code/analyze_mcmule4vecs.py" \
    --homedir "$HOMEDIR" \
    --outdir "$OUTDIR" \
    --outdir-vals "$OUTDIR_VALS" \
    --run "$RUN" \
    --savename "$SAVENAME" \
    --nbins "$NBINS"
