#!/usr/bin/env bash
set -euo pipefail
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONF="${REPO}/inputs/input_find_first_mcmule_photon_bin.conf"
[[ -f "$CONF" ]] || { echo "ERROR: missing config $CONF" >&2; exit 1; }
source "$CONF"
: "${INPUT_ROOT_GLOB:?}"
: "${ROOT_SETUP:?}"
: "${FIGUREBASE:?}"
: "${FIGUREFOLDER:?}"
[[ -f "$ROOT_SETUP" ]] || { echo "ERROR: ROOT setup missing: $ROOT_SETUP" >&2; exit 1; }
set +u
source "$ROOT_SETUP"
set -u
python3 "${REPO}/Code/find_first_mcmule_photon_bin.py" \
  --input-glob "$INPUT_ROOT_GLOB" \
  --outdir "${FIGUREBASE}/${FIGUREFOLDER}" \
  --max-mev "${MAX_PHOTON_ENERGY_MEV:-50}" \
  --bin-width-mev "${BIN_WIDTH_MEV:-1e-10}"
