#!/usr/bin/env bash
set -euo pipefail
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${scriptdir}/.." && pwd)"
config="${1:-${repo}/inputs/input_compare_tgeant_phast.conf}"
[[ -f "$config" ]] || { echo "No config: $config" >&2; exit 1; }
source "$config"
: "${mc_histogramfile:?Need MC histogram ROOT path}"
: "${phast_histogramfile:?Need PHAST histogram ROOT path}"
: "${figurebase:?Missing figurebase}"
: "${figurefoldername:?Missing figurefoldername}"
: "${root_setup:?Missing root_setup}"
[[ -f "$root_setup" ]] || { echo "ROOT setup not found: $root_setup" >&2; exit 1; }
[[ -f "$mc_histogramfile" ]] || { echo "Missing MC histogram file: $mc_histogramfile" >&2; exit 1; }
[[ -f "$phast_histogramfile" ]] || { echo "Missing PHAST histogram file: $phast_histogramfile" >&2; exit 1; }
if type module &>/dev/null; then module load gcc/13.2.0; fi
source "$root_setup"
figurefolder="${figurebase}/${figurefoldername}"
mkdir -p "$figurefolder"
root -l -b -q "${repo}/Code/plot_compare_tgeant_phast.C(\"${mc_histogramfile}\",\"${phast_histogramfile}\",\"${figurefolder}\")"
echo "Comparison figures: $figurefolder"
