#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG="${1:-${REPO}/inputs/input_merge_coralphast.conf}"
[[ -f "$CONFIG" ]] || { echo "ERROR: Missing config $CONFIG" >&2; exit 1; }
source "$CONFIG"
for var in INPUT_DIR INPUT_PATTERN EXPECTED_FILES TREE_NAME TGEANT_DIR TGEANT_PATTERN COUNT_REPORT MERGED_DIR MERGED_FILE MERGE_CHUNK_SIZE USE_EXISTING_MERGE OVERWRITE DELETE_TGEANT_AFTER_SUCCESS; do
  [[ -n "${!var:-}" ]] || { echo "ERROR: Required config value '$var' is missing" >&2; exit 1; }
done
if [[ "$DELETE_TGEANT_AFTER_SUCCESS" != 0 ]]; then
  echo "ERROR: Automatic TGEANT deletion is intentionally disabled. Set DELETE_TGEANT_AFTER_SUCCESS=0." >&2
  exit 1
fi
if ! command -v root >/dev/null 2>&1 || ! command -v hadd >/dev/null 2>&1; then
  [[ -f "${ROOT_SETUP:-}" ]] || { echo "ERROR: ROOT_SETUP missing or not found" >&2; exit 1; }
  echo "Loading ROOT from $ROOT_SETUP"
  set +u
  source "$ROOT_SETUP"
  set -u
fi
command -v root >/dev/null || { echo "ERROR: root command unavailable" >&2; exit 1; }
command -v hadd >/dev/null || { echo "ERROR: hadd command unavailable" >&2; exit 1; }
export INPUT_DIR INPUT_PATTERN EXPECTED_FILES TREE_NAME TGEANT_DIR TGEANT_PATTERN COUNT_REPORT MERGED_DIR MERGED_FILE MERGE_CHUNK_SIZE USE_EXISTING_MERGE OVERWRITE DELETE_TGEANT_AFTER_SUCCESS
python3 "${REPO}/Code/coralphast_merge_tools.py" merge
