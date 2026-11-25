#!/usr/bin/env bash
set -euo pipefail
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
REPO="$(cd "$SCRIPT_DIR/.." && pwd)"
CONFIG="${1:-${REPO}/inputs/input_count_phast_expected_events.conf}"
[[ -f "$CONFIG" ]] || { echo "ERROR: Missing config $CONFIG" >&2; exit 1; }
source "$CONFIG"
: "${TGEANT_DIR:?}" "${TGEANT_PATTERN:?}" "${EXPECTED_FILES:?}" "${COUNT_REPORT:?}"
export TGEANT_DIR TGEANT_PATTERN EXPECTED_FILES COUNT_REPORT
python3 "${REPO}/Code/coralphast_merge_tools.py" count
