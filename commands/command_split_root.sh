#!/usr/bin/env bash
# Split 10 (or nparts) independently generated NLO-CU ROOT samples into
# 1000-entry files in one flat output directory for the next TGEANT step.
# Each file is uniquely prefixed by its parent part number.
set -euo pipefail

CONFIG="${SPLIT_ROOT_CONFIG:-/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_split_root.conf}"
if [[ ! -f "$CONFIG" ]]; then
    echo "ERROR: configuration not found: $CONFIG" >&2
    exit 1
fi
source "$CONFIG"

: "${codefolder:?Missing codefolder}"
: "${runfolder:?Missing runfolder}"
: "${root_prefix:?Missing root_prefix}"
: "${outputfolder:?Missing outputfolder}"
: "${batchsize:?Missing batchsize}"
: "${nparts:?Missing nparts}"
: "${first_seed:?Missing first_seed}"

for v in batchsize nparts first_seed; do
    if ! [[ "${!v}" =~ ^[0-9]+$ ]]; then
        echo "ERROR: $v must be an integer: ${!v}" >&2
        exit 1
    fi
done
if (( batchsize < 1 || nparts < 1 )); then
    echo "ERROR: batchsize and nparts must be > 0" >&2
    exit 1
fi

# CVMFS ROOT environment (LCG setup accesses unset COMPILER variables).
if [[ "${SPLIT_ROOT_SKIP_ENV:-0}" != 1 ]]; then
    module load gcc/13.2.0
    set +u
    source /cvmfs/sft.cern.ch/lcg/views/LCG_109/x86_64-el9-gcc13-opt/setup.sh
    set -u
fi

sourcefile="$codefolder/split_root.cpp"
executable="${SPLIT_ROOT_EXE:-$codefolder/split_root}"
if [[ "${SPLIT_ROOT_SKIP_COMPILE:-0}" != 1 ]]; then
    [[ -f "$sourcefile" ]] || { echo "ERROR: missing $sourcefile" >&2; exit 1; }
    echo "Compiling $sourcefile ..."
    g++ "$sourcefile" $(root-config --cflags --libs) -o "$executable"
fi
[[ -x "$executable" ]] || { echo "ERROR: not executable: $executable" >&2; exit 1; }

mkdir -p "$outputfolder/logs"

echo "============================================================"
echo "Split ROOT NLO-CU runs"
echo "Run folder: $runfolder"
echo "Output (flat): $outputfolder"
echo "Batch size: $batchsize entries"
echo "Parts: 1 to $nparts | Seed start: $first_seed"
echo "============================================================"

new_files=0
for (( part=1; part<=nparts; part++ )); do
    seed=$((first_seed + part - 1))
    inputfile="$runfolder/part_${part}/out/${root_prefix}_part_${part}_seed${seed}.root"
    tempfolder="$outputfolder/.working_part_${part}"
    marker="$outputfolder/.split_done_part_${part}"
    log="$outputfolder/logs/split_part_${part}.log"

    if [[ ! -f "$inputfile" ]]; then
        echo "ERROR: input for part $part not found: $inputfile" >&2
        exit 1
    fi
    if [[ -f "$marker" ]]; then
        echo "[part $part/$nparts] already split, skipping ($marker)"
        continue
    fi

    # Refuse ambiguous partial reruns; do not silently overwrite data.
    shopt -s nullglob
    existing=("$outputfolder"/part_${part}_*.root)
    leftovers=("$tempfolder"/*)
    shopt -u nullglob
    if (( ${#existing[@]} > 0 || ${#leftovers[@]} > 0 )); then
        echo "ERROR: part $part has leftover output from an incomplete run." >&2
        echo "Inspect $tempfolder and $outputfolder/part_${part}_*.root" >&2
        echo "No files were overwritten." >&2
        exit 1
    fi

    mkdir -p "$tempfolder"
    echo "[part $part/$nparts] $inputfile"
    echo "  -> splitting into $tempfolder; logfile: $log"
    if ! "$executable" "$inputfile" "$tempfolder" "$batchsize" > "$log" 2>&1; then
        echo "ERROR: split_root failed for part $part; see $log" >&2
        exit 1
    fi

    shopt -s nullglob
    generated=("$tempfolder"/*.root)
    shopt -u nullglob
    if (( ${#generated[@]} == 0 )); then
        echo "ERROR: split_root returned success but produced no ROOT files for part $part" >&2
        echo "See $log" >&2
        exit 1
    fi

    for f in "${generated[@]}"; do
        target="$outputfolder/part_${part}_$(basename "$f")"
        if [[ -e "$target" ]]; then
            echo "ERROR: refusing to overwrite $target" >&2
            exit 1
        fi
        mv -- "$f" "$target"
    done
    # Leave non-ROOT auxiliary files in tempfolder if splitter produced any.
    rmdir "$tempfolder" 2>/dev/null || true
    printf 'input=%s\nbatchsize=%d\nfiles=%d\n' \
        "$inputfile" "$batchsize" "${#generated[@]}" > "$marker"
    echo "  OK: ${#generated[@]} ROOT batch files"
    new_files=$((new_files + ${#generated[@]}))
done

echo "============================================================"
echo "Finished. Newly created files: $new_files"
echo "Total batch ROOT files: $(find "$outputfolder" -maxdepth 1 -type f -name '*.root' | wc -l)"
echo "Splitting logs: $outputfolder/logs"
echo "For the next TGEANT step use INPUTDIR=$outputfolder"
echo "============================================================"
