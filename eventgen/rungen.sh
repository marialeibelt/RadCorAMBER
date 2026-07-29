#!/bin/bash
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/gen.sh

export MCMULE_JOBS=5

generate 0 N 10 29437 &
generate 1 N 01 75447 &
wait
generate 1 C 01 67697 &
