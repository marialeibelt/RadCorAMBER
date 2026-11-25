#!/usr/bin/env bash
set -euo pipefail
REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONF="${REPO}/inputs/input_merge_coralphast_files.conf"
[[ -f "$CONF" ]] || { echo "ERROR: Missing configuration: $CONF" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONF"
: "${ROOT_SETUP:?}"
: "${MERGED_DIR:?}"
: "${MERGED_FILE:?}"
[[ -f "$ROOT_SETUP" ]] || { echo "ERROR: ROOT setup missing: $ROOT_SETUP" >&2; exit 1; }
# LCG setup uses optional environment variables; preserve strictness after sourcing.
set +u
source "$ROOT_SETUP"
set -u
[[ -f "${REPO}/Code/merge_coralphast_files.py" ]] || {
  echo "ERROR: Code/merge_coralphast_files.py is missing" >&2; exit 1;
}
mkdir -p "$MERGED_DIR"
export MERGE_FILES_INPUT_LIST="${INPUT_FILES[*]:-}"  # only for diagnostics; actual list via tempfile
# Serialize exact bash array values without shell re-evaluation/globbing.
MANIFEST_INPUT="$(mktemp "${MERGED_DIR}/.merge_selected_input_XXXXXX")"
trap 'rm -f "$MANIFEST_INPUT"' EXIT
for f in "${INPUT_FILES[@]:-}"; do
  [[ -n "$f" ]] && printf '%s\n' "$f" >> "$MANIFEST_INPUT"
done
python3 "${REPO}/Code/merge_coralphast_files.py" \
  --input-list "$MANIFEST_INPUT" \
  --glob-list "$(printf '%s\n' "${INPUT_GLOBS[@]:-}")" \
  --extra-list "${LIST_FILE:-}" \
  --out "${MERGED_DIR}/${MERGED_FILE}" \
  --tree "${TREE_NAME:-UserEvent300/PRM}" \
  --chunk-size "${MERGE_CHUNK_SIZE:-30}" \
  --overwrite "${OVERWRITE:-0}" \
  --keep-chunks "${KEEP_TEMP_CHUNKS:-0}" \
  --expected-files "${EXPECTED_FILES:-0}" \
  --expected-entries "${EXPECTED_TREE_ENTRIES:-0}"
