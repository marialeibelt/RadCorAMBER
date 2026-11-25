#!/usr/bin/env bash
set -euo pipefail
scriptdir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
repo="$(cd "${scriptdir}/.." && pwd)"
config="${1:-${repo}/inputs/input_analyze_tgeant_mc_truth.conf}"
[[ -f "$config" ]] || { echo "No config: $config" >&2; exit 1; }
source "$config"
: "${input_glob:?Missing input_glob}"
: "${mc_histogramfile:?Missing mc_histogramfile}"
: "${root_setup:?Missing root_setup}"
[[ -f "$root_setup" ]] || { echo "ROOT setup not found: $root_setup" >&2; exit 1; }
if type module &>/dev/null; then module load gcc/13.2.0; fi
set +u
source "$root_setup"
set -u
python3 -c 'import ROOT; print("ROOT version:", ROOT.gROOT.GetVersion())' || {
 echo "PyROOT not available in this LCG environment" >&2; exit 1;
}
mapfile -t files < <(compgen -G "$input_glob" | sort -V)
((${#files[@]}>0)) || { echo "No files found for: $input_glob" >&2; exit 1; }
if ((${max_files:-0}>0)); then files=("${files[@]:0:${max_files}}"); fi
mkdir -p "$(dirname "$mc_histogramfile")"
echo "Processing ${#files[@]} TGEANT files -> $mc_histogramfile"
python3 "${repo}/Code/analyze_tgeant_mc_truth.py" \
  --output "$mc_histogramfile" \
  --max-events "${max_events:-0}" \
  --nominal-beam-energy "${nominal_beam_energy:-100}" \
  --muon-min-p "${muon_min_p:-70}" \
  "${files[@]}"
