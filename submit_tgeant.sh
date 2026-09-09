#!/bin/bash

#$ -S /bin/bash
#$ -N tgeant
#$ -l long=TRUE
#$ -l h_pmem=12G
#$ -l h_vmem=12G

#$ -wd /nfs/momos/user/mleibelt/TGEANT_runs

#$ -o /dev/null
#$ -j y

/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor/commands/command_tgeant_amber.sh

exit 0