#!/bin/bash
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/gen.sh

export MCMULE_JOBS=5

generate 0 N 10 29374 &
generate 1 N 01 75274 &
wait
generate 1 CU 01 67374 &
