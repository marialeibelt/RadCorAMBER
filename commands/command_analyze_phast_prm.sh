#!/usr/bin/env bash
set -euo pipefail

# ============================================================
# Read input configuration
# ============================================================
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
PACKAGE_DIR="$(cd "${SCRIPT_DIR}/.." && pwd)"
CONFIG="${1:-${PACKAGE_DIR}/inputs/input_analyze_phast_prm.conf}"

if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Config file not found: $CONFIG" >&2
    exit 1
fi
source "$CONFIG"
: "${phast_histogramfile:?Set phast_histogramfile to a PHAST histograms ROOT file}"
: "${figurebase:?Missing figurebase}"
: "${figurefoldername:?Missing figurefoldername}"

# ============================================================
# Paths / checks
# ============================================================
PLOT_MACRO="${PACKAGE_DIR}/Code/plot_analyze_phast_prm.C"
figurefolder="${figurebase}/${figurefoldername}"
if [[ ! -f "$phast_histogramfile" ]]; then
    echo "ERROR: ROOT input file does not exist: $phast_histogramfile" >&2
    exit 1
fi
if [[ ! -f "$PLOT_MACRO" ]]; then
    echo "ERROR: Plot macro does not exist: $PLOT_MACRO" >&2
    exit 1
fi
mkdir -p "$figurefolder"

# ============================================================
# Environment
# ============================================================
if type module >/dev/null 2>&1; then module load gcc/13.2.0; fi
if [[ -n "${root_setup:-}" ]]; then
    if [[ ! -f "$root_setup" ]]; then
        echo "ERROR: ROOT environment missing: $root_setup" >&2
        exit 1
    fi
    source "$root_setup"
fi
if ! command -v root >/dev/null; then
    echo "ERROR: ROOT is not available" >&2
    exit 1
fi

# ============================================================
# Plot PHAST histograms
# ============================================================
echo "=============================================="
echo "Analyze PHAST / CORAL"
echo "=============================================="
echo "Histogram input: $phast_histogramfile"
echo "Figure folder:   $figurefolder"
echo "=============================================="
root -l -b -q "${PLOT_MACRO}(\"${phast_histogramfile}\",\"${figurefolder}\")"
echo "Finished: $figurefolder"
