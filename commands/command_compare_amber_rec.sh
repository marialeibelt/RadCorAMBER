#!/usr/bin/env bash
set -euo pipefail
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONF="${REPO}/inputs/input_compare_amber_rec.conf"
[[ -f "$CONF" ]] || { echo "ERROR: Config missing: $CONF" >&2; exit 1; }
source "$CONF"
: "${BEFORE_ROOT_FILE:?}"
: "${AFTER_ROOT_FILE:?}"
: "${FIGUREBASE:?}"
: "${FIGUREFOLDER:?}"
: "${ROOT_SETUP:?}"
for file in "$BEFORE_ROOT_FILE" "$AFTER_ROOT_FILE" "$ROOT_SETUP"; do
    [[ -f "$file" ]] || { echo "ERROR: Missing $file" >&2; exit 1; }
done
set +u
source "$ROOT_SETUP"
set -u
export MPLBACKEND=Agg
OUTDIR="${FIGUREBASE}/${FIGUREFOLDER}"
mkdir -p "$OUTDIR"
python3 "${REPO}/Code/compare_amber_rec.py" \
  --before "$BEFORE_ROOT_FILE" --after "$AFTER_ROOT_FILE" --outdir "$OUTDIR" \
  --before-tree "${BEFORE_TREE:-TGEANTTruth}" --after-tree "${AFTER_TREE:-UserEvent300/PRM}" \
  --q2-log-min "${Q2_LOG_MIN:-1e-10}" \
  --q2-log-max "${Q2_LOG_MAX:-0.05}" \
  --q2-log-bins "${Q2_LOG_BINS:-450}" \
  --photon-diagnostic "${PHOTON_DIAGNOSTIC:-1}"
echo "Comparison plots and counts: $OUTDIR"
