#!/bin/bash
source pymule-generate.sh

export MCMULE_JOBS=5

generate 0 N 10 29517 &
generate 1 N 01 75329 &
wait

generate 1 C 01 68990 &
wait
