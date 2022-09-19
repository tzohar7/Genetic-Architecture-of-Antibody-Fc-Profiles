#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of Tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --job-name=postGWAS
#SBATCH --array=0-202

###########################################################################################
# Setup and Parallelization ###############################################################
###########################################################################################

printf "Job starting...\n\n"

INPUTDIR=~/data/inputdir
FILES=($(ls $INPUTDIR/*.txt))

CURRENT_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}

###########################################################################################
# Change headers ##########################################################################
###########################################################################################

sed '1{ s/Chr/CHR/; s/ChrPos/BP/; s/PValue/P/; s/GenDist//;}' $CURRENT_FILE | awk '$1=$1' > $CURRENT_FILE.tmp

###########################################################################################
# Clumping ################################################################################
###########################################################################################

module load plink

WORDTOREMOVE="~/data/location_name"

PHENO="${CURRENT_FILE//".txt"/}"
PHENO="${PHENO//$WORDTOREMOVE/}"

OUTPUT="~/data/clumped/"
OUTPUT+=$PHENO

plink --bfile ~data/bedfiles/filtered/genome --clump $CURRENT_FILE --clump-p1 0.0001 --clump-p2 0.01 --clump-r2 0.50 --clump-kb 100 --out $OUTPUT

rm $OUTPUT.log
rm $OUTPUT.nosex
rm $CURRENT_FILE.tmp
