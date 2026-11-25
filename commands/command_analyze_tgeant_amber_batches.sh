#!/usr/bin/env bash
# Process AMBER TGEANT ASCII batches to ROOT, merge in resumable chunks,
# produce separate style-matched figures, and inventory the ROOT objects.
# This script NEVER removes the original .tgeant.gz files.
set -euo pipefail

REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONFIG="${REPO}/inputs/input_analyze_tgeant_amber_batches.conf"
ANALYZE_SOURCE="${REPO}/Code/analyze_AMBER_eventgen.cpp"
ANALYZE_EXECUTABLE="${REPO}/analyze_AMBER_eventgen"
PLOT_MACRO="${REPO}/Code/plot_analyze_AMBER_eventgen.cpp"
STYLED_PLOTTER="${REPO}/Code/plot_tgeant_amber_styled.py"

[[ -f "$CONFIG" ]] || { echo "ERROR: Missing config $CONFIG" >&2; exit 1; }
# shellcheck source=/dev/null
source "$CONFIG"

: "${inputdir:?missing inputdir}"
: "${file_glob:?missing file_glob}"
: "${histogramfile:?missing histogramfile}"
: "${figurebase:?missing figurebase}"
: "${figurefoldername:?missing figurefoldername}"
chunk_size="${chunk_size:-50}"
cache_dir="${cache_dir:-${histogramfile%.root}_analysis_chunks}"
resume="${resume:-1}"
make_styled_plots="${make_styled_plots:-1}"
make_legacy_plots="${make_legacy_plots:-0}"
expected_files="${expected_files:-0}"
expected_events="${expected_events:-0}"

[[ "$chunk_size" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: chunk_size must be positive integer" >&2; exit 1; }
[[ "$expected_files" =~ ^[0-9]+$ ]] || { echo "ERROR: expected_files must be integer" >&2; exit 1; }
[[ "$expected_events" =~ ^[0-9]+$ ]] || { echo "ERROR: expected_events must be integer" >&2; exit 1; }
[[ -d "$inputdir" ]] || { echo "ERROR: Missing input directory $inputdir" >&2; exit 1; }
[[ -f "$ANALYZE_SOURCE" ]] || { echo "ERROR: Missing analyzer source $ANALYZE_SOURCE" >&2; exit 1; }
[[ "$make_styled_plots" != 1 || -f "$STYLED_PLOTTER" ]] || { echo "ERROR: Missing plotter $STYLED_PLOTTER" >&2; exit 1; }
[[ "$make_legacy_plots" != 1 || -f "$PLOT_MACRO" ]] || { echo "ERROR: Missing macro $PLOT_MACRO" >&2; exit 1; }

figurefolder="${figurebase}/${figurefoldername}"
mkdir -p "$cache_dir" "$(dirname "$histogramfile")" "$figurefolder"

# Avoid two sessions writing the same cache/final ROOT file simultaneously.
exec 8>"${cache_dir}/.analysis.lock"
flock -n 8 || { echo "ERROR: another analysis uses this cache_dir: $cache_dir" >&2; exit 1; }

# Discover files in stable, numeric/natural order; no shell word-splitting of paths.
mapfile -d '' -t inputfiles < <(find "$inputdir" -maxdepth 1 -type f -name "$file_glob" -print0 | sort -zV)
nfiles=${#inputfiles[@]}
(( nfiles > 0 )) || { echo "ERROR: No TGEANT files match $file_glob in $inputdir" >&2; exit 1; }
if (( expected_files > 0 && nfiles != expected_files )); then
  echo "ERROR: expected $expected_files .tgeant.gz files, found $nfiles." >&2
  echo "Finish TGEANT before starting analysis, or correct expected_files." >&2
  exit 1
fi

# Preserve an exact input inventory. File sizes + modification times catch
# accidental input changes when resuming; no expensive hashing of gzip data.
manifest_tmp="$(mktemp "${cache_dir}/manifest.XXXXXX")"
for f in "${inputfiles[@]}"; do
  stat -c '%n	%s	%Y' "$f" >> "$manifest_tmp"
done
fingerprint="$( { printf 'chunk_size=%s\n' "$chunk_size"; sha256sum "$ANALYZE_SOURCE"; cat "$manifest_tmp"; } | sha256sum | awk '{print $1}')"

if [[ -f "${cache_dir}/input_fingerprint.sha256" ]]; then
  old_fingerprint="$(cat "${cache_dir}/input_fingerprint.sha256")"
  if [[ "$old_fingerprint" != "$fingerprint" ]]; then
    rm -f "$manifest_tmp"
    echo "ERROR: Input list, analyzer source or chunk size changed since cached chunks were generated." >&2
    echo "Use a new cache_dir, or deliberately clear the old cache and rerun." >&2
    exit 1
  fi
elif compgen -G "${cache_dir}/chunk_*.root" > /dev/null; then
  rm -f "$manifest_tmp"
  echo "ERROR: existing cached chunks have no fingerprint. Choose a fresh cache_dir." >&2
  exit 1
fi
mv -f "$manifest_tmp" "${cache_dir}/input_manifest.tsv"
printf '%s\n' "$fingerprint" > "${cache_dir}/input_fingerprint.sha256"

# Collect run log in the persistent cache.
exec > >(tee -a "${cache_dir}/analysis.log") 2>&1
printf '\n=== TGEANT analysis started %s ===\n' "$(date -Is)"
echo "Input directory: $inputdir"
echo "Matched files:   $nfiles"
echo "Output ROOT:     $histogramfile"
echo "Plot directory:  $figurefolder"
echo "Chunk size:      $chunk_size"
echo "Persistent cache:$cache_dir"
echo "Archive contains primary MC truth and histograms, NOT detector hits."
echo "Expected primary-event entries: $expected_events"

# ROOT environment: LCG setup can access unset variables; temporarily relax -u.
module load gcc/13.2.0
set +u
source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
set -u

command -v root-config >/dev/null || { echo "ERROR: root-config missing" >&2; exit 1; }
command -v hadd >/dev/null || { echo "ERROR: hadd missing" >&2; exit 1; }

echo "Compiling $ANALYZE_SOURCE ..."
g++ "$ANALYZE_SOURCE" $(root-config --cflags --libs) -lz -o "$ANALYZE_EXECUTABLE"

# Limit temporary space to one chunk. Completed merged chunks remain in the
# persistent cache and can be reused if a subsequent batch fails.
tmpdir=''
cleanup() { [[ -z "$tmpdir" || ! -d "$tmpdir" ]] || rm -rf -- "$tmpdir"; }
trap cleanup EXIT

chunk_roots=()
chunk_index=0
for ((start=0; start<nfiles; start+=chunk_size)); do
  stop=$((start + chunk_size))
  (( stop > nfiles )) && stop=$nfiles
  printf -v chunk_tag '%05d' "$chunk_index"
  chunkroot="${cache_dir}/chunk_${chunk_tag}.root"

  if [[ "$resume" == 1 && -s "$chunkroot" ]]; then
    echo "Chunk $chunk_index ($((${start}+1))-$stop / $nfiles): using existing $chunkroot"
    chunk_roots+=("$chunkroot")
    ((chunk_index+=1))
    continue
  fi

  echo "Chunk $chunk_index ($((${start}+1))-$stop / $nfiles): analyzing ..."
  tmpdir="$(mktemp -d "${TMPDIR:-/tmp}/tgeant_chunk_XXXXXX")"
  batch_histograms=()
  for ((j=start; j<stop; ++j)); do
    printf -v btag '%05d' "$((j+1))"
    batchroot="${tmpdir}/batch_${btag}.root"
    "$ANALYZE_EXECUTABLE" "${inputfiles[$j]}" "$batchroot"
    [[ -s "$batchroot" ]] || { echo "ERROR: Analyzer produced no output for ${inputfiles[$j]}" >&2; exit 1; }
    batch_histograms+=("$batchroot")
  done

  hadd -f "${tmpdir}/merged_chunk.root" "${batch_histograms[@]}"
  [[ -s "${tmpdir}/merged_chunk.root" ]] || { echo "ERROR: chunk merge failed" >&2; exit 1; }
  # Same-filesystem atomic publish: do not consider partially merged files done.
  staged="${cache_dir}/.chunk_${chunk_tag}.root.$$"
  cp "${tmpdir}/merged_chunk.root" "$staged"
  mv -f "$staged" "$chunkroot"
  rm -rf -- "$tmpdir"
  tmpdir=''

  chunk_roots+=("$chunkroot")
  ((chunk_index+=1))
done

echo "Merging ${#chunk_roots[@]} completed chunks into ONE ROOT file ..."
staged_final="${histogramfile%.root}.inprogress.$$.root"
trap 'cleanup; rm -f -- "${staged_final:-}"' EXIT
hadd -f "$staged_final" "${chunk_roots[@]}"
[[ -s "$staged_final" ]] || { echo "ERROR: final ROOT file is empty" >&2; exit 1; }
mv -f "$staged_final" "$histogramfile"

# Structural and completeness check on the merged TTree. Only parseable
# primary truth is archived, but we should at least have one row for every
# TGEANT #EVENT block. A missing event or failed merge is a hard error.
python3 - "$histogramfile" "$expected_events" <<'PY_ROOT_CHECK'
import sys
import ROOT
file_name, expected_s = sys.argv[1:]
expected = int(expected_s)
f = ROOT.TFile.Open(file_name, "READ")
if not f or f.IsZombie():
    raise SystemExit(f"ERROR: cannot open merged archive {file_name}")
t = f.Get("TGEANTTruth")
h = f.Get("h_event_summary")
if not t or not h:
    raise SystemExit("ERROR: missing TGEANTTruth or event summary in archive")
n = int(t.GetEntries())
nh = int(round(h.GetBinContent(1)))
print(f"VALIDATION: merged TGEANTTruth entries = {n:,}")
print(f"VALIDATION: merged #EVENT blocks      = {nh:,}")
if n != nh:
    raise SystemExit("ERROR: TTree entry count and event block count differ")
if expected and n != expected:
    raise SystemExit(f"ERROR: expected {expected:,} events, archived {n:,}")
if nh > 0:
    primary = int(round(h.GetBinContent(2)))
    valid_muon = int(round(h.GetBinContent(4)))
    print(f"VALIDATION: parsed primary block in {primary:,} events ({primary/nh:.1%})")
    print(f"VALIDATION: usable muon kinematics in {valid_muon:,} events ({valid_muon/nh:.1%})")
    if valid_muon == 0:
        raise SystemExit("ERROR: no usable muon kinematics; verify TGEANT particle format/units")
f.Close()
PY_ROOT_CHECK

# Save provenance alongside the archive; duplicate this manifest to the
# output directory so it remains available even if the cache is later purged.
cp "${cache_dir}/input_manifest.tsv" "${histogramfile%.root}_input_manifest.tsv"
cp "${cache_dir}/input_fingerprint.sha256" "${histogramfile%.root}_input_fingerprint.sha256"

if [[ "$make_styled_plots" == 1 ]]; then
  echo "Producing separate styled ROOT-histogram plots and ROOT-object inventory ..."
  python3 "$STYLED_PLOTTER" "$histogramfile" "$figurefolder"
fi
if [[ "$make_legacy_plots" == 1 ]]; then
  echo "Running original plot macro ..."
  mkdir -p "${figurefolder}/legacy"
  root -l -b -q "${PLOT_MACRO}(\"${histogramfile}\",\"${figurefolder}/legacy\")"
fi

echo "Final file: $histogramfile"
echo "Plot folder: $figurefolder"
echo "Analyzed inputs: $nfiles"
echo "NOTE: This archive contains parsed primary MC records, not detector hits."
echo "      Do NOT delete .tgeant.gz files until CORAL/PHAST and all validation are complete."
printf '=== Finished %s ===\n' "$(date -Is)"
