#!/usr/bin/env bash
# Submit one McMule NLO SGE job.
# Usage: bash submit_nlo_sge.sh [CU|NCU|N] [nenter] [itmx] [seed]
#   CU  : only CU (cellular resampling + unweighting)  <- default, what we need
#   NCU : weighted N run first, then CU (old behaviour; N is NOT needed for CU)
#   N   : only the weighted N run
# Example: bash submit_nlo_sge.sh CU 5000 20 10010
#
# Every seed is checked against and recorded in inputs/used_seeds_evtgen.txt,
# so a seed can never be used twice for event generation.
set -euo pipefail

REPO="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
CONFIG="$REPO/inputs/input_evtgen_mcmule.conf"
JOB_SCRIPT="$REPO/commands/mcmule_eventgen/run_nlo_sge.job"
GENERATOR="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator/build/src/mcmule"
SEED_REGISTRY="${SEED_REGISTRY:-$REPO/inputs/used_seeds_evtgen.txt}"

[[ -r "$CONFIG" ]] || { echo "ERROR: Missing $CONFIG" >&2; exit 1; }
[[ -r "$JOB_SCRIPT" ]] || { echo "ERROR: Missing $JOB_SCRIPT" >&2; exit 1; }
[[ -x "$GENERATOR" ]] || { echo "ERROR: Missing executable $GENERATOR" >&2; exit 1; }
[[ -f "$SEED_REGISTRY" ]] || { echo "ERROR: Missing seed registry $SEED_REGISTRY" >&2; exit 1; }
command -v qsub >/dev/null 2>&1 || { echo "ERROR: qsub not available" >&2; exit 1; }

# Existing configuration supplies the xi-dependent integration filenames and source folder.
# shellcheck source=/dev/null
source "$CONFIG"

MODE="${1:-CU}"
NENTER="${2:-50000}"
ITMX="${3:-20}"
SEED="${4:-}"
[[ "$MODE" == "CU" || "$MODE" == "NCU" || "$MODE" == "N" ]] || {
    echo "ERROR: Mode must be CU, NCU or N." >&2; exit 1;
}
[[ -n "$SEED" ]] || { echo "ERROR: Seed is required (see $SEED_REGISTRY for used seeds)." >&2; exit 1; }
[[ "$NENTER" =~ ^[1-9][0-9]*$ && "$ITMX" =~ ^[1-9][0-9]*$ && "$SEED" =~ ^[1-9][0-9]*$ ]] || {
    echo "ERROR: nenter, itmx and seed must be positive integers" >&2; exit 1;
}
[[ "$xi" == "01" || "$xi" == "10" ]] || { echo "ERROR: Unsupported xi='$xi'" >&2; exit 1; }

# Stages and their seeds. In NCU mode the N seed is offset to stay unique.
RUN_N=0; RUN_CU=0
case "$MODE" in
    CU)  RUN_CU=1; N_SEED=0 ;;
    NCU) RUN_N=1; RUN_CU=1; N_SEED=$((SEED + 1000000)) ;;
    N)   RUN_N=1; N_SEED="$SEED" ;;
esac

seed_used() { awk -v s="$1" '!/^#/ && $1 == s {found=1} END {exit !found}' "$SEED_REGISTRY"; }
check_seeds=()
[[ "$RUN_CU" -eq 1 ]] && check_seeds+=("$SEED")
[[ "$RUN_N" -eq 1 ]] && check_seeds+=("$N_SEED")
for s in "${check_seeds[@]}"; do
    if seed_used "$s"; then
        echo "ERROR: Seed $s was already used (see $SEED_REGISTRY). Choose a new seed." >&2
        exit 1
    fi
done

# Output folder: by default next to the integration it uses.
# Override using: ANA_FOLDER='my/new/path' bash submit_nlo_sge.sh CU 5000 20 10010
ANA_FOLDER="${ANA_FOLDER:-$anafolder/evtgen_${MODE}_${NENTER}_seed${SEED}}"
case "$ANA_FOLDER" in
    /*|*..*) echo "ERROR: ANA_FOLDER must be a relative path without '..'" >&2; exit 1 ;;
esac
SRC_DIR="$basefolder/$anafolder"
RUN_DIR="$basefolder/$ANA_FOLDER"

[[ -f "$SRC_DIR/user.so" ]] || { echo "ERROR: Missing source user.so: $SRC_DIR/user.so" >&2; exit 1; }
if [[ "$xi" == "10" ]]; then
    INPUT_NLO0="$mcmuleNLO0_xi10"
    INPUT_R="$mcmuleR_xi10"
else
    INPUT_NLO0="$mcmuleNLO0_xi01"
    INPUT_R="$mcmuleR_xi01"
fi

for f in "$INPUT_NLO0" "$INPUT_R"; do
    [[ -f "$SRC_DIR/input/$f" ]] || { echo "ERROR: Missing integration: $SRC_DIR/input/$f" >&2; exit 1; }
done

mkdir -p "$RUN_DIR/input" "$RUN_DIR/out" "$RUN_DIR/log/sge"
# Keep the exact compiled userfile and integration data as in the original run.
for f in "user.so" "input/$INPUT_NLO0" "input/$INPUT_R"; do
    if [[ -e "$RUN_DIR/$f" ]]; then
        cmp -s "$SRC_DIR/$f" "$RUN_DIR/$f" || {
            echo "ERROR: Existing $RUN_DIR/$f differs from the source (refusing to overwrite)." >&2; exit 1;
        }
    else
        cp -p "$SRC_DIR/$f" "$RUN_DIR/$f"
    fi
done

OUT_REL_N="out/gen-N-1-${xi}-${N_SEED}.lhe"
OUT_REL_CU="out/gen-CU-1-${xi}-${SEED}.lhe"
LOG_FILE_N="$RUN_DIR/log/gen-N-1-${xi}-${N_SEED}.log"
LOG_FILE_CU="$RUN_DIR/log/gen-CU-1-${xi}-${SEED}.log"
for existing in "$RUN_DIR/$OUT_REL_N" "$RUN_DIR/$OUT_REL_CU"; do
    if [[ -e "$existing" ]]; then
        echo "ERROR: Output LHE already exists: $existing" >&2
        echo "Use a new seed or ANA_FOLDER to protect existing data." >&2
        exit 1
    fi
done

# Optional scheduler choices (must match locally configured queue/parallel environment):
#   SGE_QUEUE=medium
#   SGE_PE=smp SGE_SLOTS=5  (check available PE using: qconf -spl)
#   SGE_MEMORY_REQUEST='mem_free=2G' (check resource names with cluster admins)
SLOTS="${SGE_SLOTS:-1}"
[[ "$SLOTS" =~ ^[1-9][0-9]*$ ]] || { echo "ERROR: SGE_SLOTS must be a positive integer" >&2; exit 1; }
qsub_args=(-cwd -V -j y -N "mcNLO_${MODE}_${SEED}" -o "$RUN_DIR/log/sge")
if [[ -n "${SGE_QUEUE:-}" ]]; then qsub_args+=(-q "$SGE_QUEUE"); fi
if [[ "$SLOTS" -gt 1 ]]; then
    [[ -n "${SGE_PE:-}" ]] || {
        echo "ERROR: For SGE_SLOTS>1, set a valid SGE_PE (check qconf -spl)." >&2; exit 1;
    }
    qsub_args+=(-pe "$SGE_PE" "$SLOTS")
fi
if [[ -n "${SGE_MEMORY_REQUEST:-}" ]]; then
    qsub_args+=(-l "$SGE_MEMORY_REQUEST")
fi

export RUN_DIR MODE RUN_N RUN_CU NENTER ITMX SEED N_SEED OUT_REL_N OUT_REL_CU LOG_FILE_N LOG_FILE_CU GENERATOR INPUT_NLO0 INPUT_R xi
# qsub -V passes these variables from the submission shell to the job.
echo "Submitting NLO job: mode=$MODE, nenter=$NENTER, itmx=$ITMX, seed=$SEED, xi=$xi"
echo "Integration: $SRC_DIR/input/{$INPUT_NLO0,$INPUT_R}"
[[ "$RUN_N" -eq 1 ]] && echo "STAGE N : seed=$N_SEED, LHE=$RUN_DIR/$OUT_REL_N"
[[ "$RUN_CU" -eq 1 ]] && echo "STAGE CU: seed=$SEED, LHE=$RUN_DIR/$OUT_REL_CU"
echo "Scheduler logs: $RUN_DIR/log/sge/"
echo "Requested SGE slots: $SLOTS"
qsub "${qsub_args[@]}" "$JOB_SCRIPT"

# Record the seeds only after a successful submission.
now="$(date -Is)"
[[ "$RUN_CU" -eq 1 ]] && printf '%s\tCU\t%s\t%s\n' "$SEED" "$now" "$RUN_DIR" >> "$SEED_REGISTRY"
[[ "$RUN_N" -eq 1 ]] && printf '%s\tN\t%s\t%s\n' "$N_SEED" "$now" "$RUN_DIR" >> "$SEED_REGISTRY"
exit 0
