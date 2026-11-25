#!/bin/bash
#$ -S /bin/bash
#$ -N McMule_Gen4vecs
#$ -q long_epyc,long_xeon2,long_xeonsb
#$ -l long_queue=TRUE
#$ -wd /nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor
#$ -V
#$ -j y

set -e

# Project directories
PROJECT="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/AMBER_RadCor"
MCMULE="/nfs/freenas/tuph/e18/project/prm/mleibelt/AMBER_Repo/mcmule-event-generator"

# Read input
source "$PROJECT/inputs/input_mcmule_gen4vecs.conf"

# Userfile to compile
USERFILE_PATH="$PROJECT/userfiles/$USERFILE"

# Check userfile
if [[ ! -f "$USERFILE_PATH" ]]; then
    echo "ERROR: Userfile not found: $USERFILE_PATH"
    exit 1
fi

# Create output and input directories
mkdir -p "$OUTPUTFOLDER/input"

# Go to individual output folder
cd "$OUTPUTFOLDER" || exit 1

# Compile selected userfile
gfortran -fPIC \
    -I "$MCMULE/build/src/libmcmule.so.p" \
    -shared \
    -o user.so \
    "$USERFILE_PATH"

# Start McMule
pymule batch shell menu-mp2mp.menu