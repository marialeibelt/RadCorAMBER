#!/usr/bin/env bash
# Run the new PyROOT-based analysis on ten CU ROOT files.
# Does NOT compile or modify analyze_root.cpp or merge ROOT files.
set -euo pipefail

REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONFIG="${ROOT_ANALYSIS_CONFIG:-$REPO/inputs/input_analyze_mcmule_evtgen.conf}"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: Missing config: $CONFIG" >&2
    exit 1
fi
# shellcheck source=/dev/null
source "$CONFIG"

: "${codefolder:?Missing codefolder}"
: "${basefolder:?Missing basefolder}"
: "${figure_folder:?Missing figure_folder}"
: "${beam_energy_gev:?Missing beam_energy_gev}"

# The LCG view references unset shell variables in its setup, e.g. COMPILER.
# Temporarily disable nounset when sourcing it, restore afterward.
module load gcc/13.2.0
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
set -u

# Fail early with a useful message if PyROOT or Python plotting deps are absent.
python3 -c 'import ROOT, numpy, matplotlib' || {
    echo "ERROR: LCG Python needs ROOT, numpy, matplotlib. Try python3 -c 'import ROOT, numpy, matplotlib'" >&2
    exit 1
}

script="$codefolder/analyze_mcmule_evtgen.py"
[[ -f "$script" ]] || { echo "ERROR: Missing analysis script: $script" >&2; exit 1; }
mkdir -p "$figure_folder"

args=(
    --outdir "$figure_folder"
    --beam-energy-gev "$beam_energy_gev"
    --input-frame "${input_frame:-auto}"
    --lo-frame "${lo_frame:-auto}"
    --threads "${threads:-1}"
)
if [[ -n "${nlo_root:-}" ]]; then
    args+=(--nlo-root "$nlo_root")
else
    : "${batch_folder:?Missing batch_folder}" "${n_runs:?Missing n_runs}" "${base_seed:?Missing base_seed}" "${xi:?Missing xi}" "${output_prefix:?Missing output_prefix}"
    args+=(--batch-root "$basefolder/$batch_folder" --n-runs "$n_runs" --base-seed "$base_seed" --xi "$xi" --prefix "$output_prefix")
fi
if [[ -n "${lo_root:-}" ]]; then
    args+=(--lo-root "$lo_root")
fi
if [[ "${allow_missing:-0}" == 1 ]]; then
    args+=(--allow-missing)
fi

echo "Analysing NLO ROOT events"
echo "Output: $figure_folder"
# Collect the full terminal output without putting any entry counts on the PNGs.
python3 "$script" "${args[@]}" 2>&1 | tee "$figure_folder/analyze_mcmule_evtgen.log"
