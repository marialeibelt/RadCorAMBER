#!/usr/bin/env bash
set -euo pipefail
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONF="${REPO}/inputs/input_analyze_coralphast.conf"
[[ -f "$CONF" ]] || { echo "ERROR: Cannot find $CONF" >&2; exit 1; }
source "$CONF"
: "${PHAST_ROOT_FILE:?}"
: "${ROOT_SETUP:?}"
: "${FIGUREBASE:?}"
: "${FIGUREFOLDER:?}"
DISTRIBUTION_Y_SCALE="${DISTRIBUTION_Y_SCALE:-log}"
CUTFLOW_Y_SCALE="${CUTFLOW_Y_SCALE:-linear}"
[[ -f "$PHAST_ROOT_FILE" ]] || { echo "ERROR: ROOT input missing: $PHAST_ROOT_FILE" >&2; exit 1; }
[[ -f "$ROOT_SETUP" ]] || { echo "ERROR: ROOT environment missing: $ROOT_SETUP" >&2; exit 1; }
set +u
source "$ROOT_SETUP"
set -u
export MPLBACKEND=Agg
OUTDIR="${FIGUREBASE}/${FIGUREFOLDER}"
mkdir -p "$OUTDIR"
echo "INFO: Input config: $CONF"
echo "INFO: Distribution y-axis: $DISTRIBUTION_Y_SCALE (cutflow: $CUTFLOW_Y_SCALE)"
echo "INFO: Plotter: ${REPO}/Code/analyze_coralphast.py"
python3 "${REPO}/Code/analyze_coralphast.py" \
    --input "$PHAST_ROOT_FILE" --outdir "$OUTDIR" \
    --distribution-yscale "$DISTRIBUTION_Y_SCALE" \
    --cutflow-yscale "$CUTFLOW_Y_SCALE"
echo "Plots and event_statistics.txt: $OUTDIR"
