#!/usr/bin/env bash
set -euo pipefail
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
source "${REPO}/inputs/input_analyze_mcmule4vecs.conf"
mkdir -p "$OUTDIR" "$OUTDIR_VALS"
python3 "${REPO}/Code/analyze_mcmule_meeting.py" \
  --rundir "$RUNDIR" --run "$RUN" \
  --outdir "$OUTDIR" --outdir-vals "$OUTDIR_VALS" \
  --savename "$SAVENAME" --pieces "${PIECES[@]}" "$@"