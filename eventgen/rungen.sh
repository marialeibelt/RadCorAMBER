#!/bin/bash
source /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/eventgen/gen.sh

export MCMULE_JOBS=5

generate 0 N 10 29517 &
#wait
#generate 1 N 01 75329 &
#wait
#generate 1 C 01 68990 &
