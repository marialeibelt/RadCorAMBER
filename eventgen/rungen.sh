#!/bin/bash

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

source "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/inputs/input_evtgen_mcmule.conf"
source "/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/gen.sh"

export MCMULE_JOBS=5


# LO
generate 0 N "$xi" "$seed_LO" &

# NLO
generate 1 N "$xi" "$seed_NLO" &

wait

# NLO with CRES + unweighting
generate 1 CU "$xi" "$seed_NLO_CU" &
