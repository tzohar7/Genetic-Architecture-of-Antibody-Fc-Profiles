#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of Tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --job-name=preproc
#SBATCH --array=0-22

printf "jobstarting...\n\n"

module load plink

# Change to your directories
INPUTDIR=~/data/vcffiles
WORKDIR=~/data/vcffiles
cd $WORKDIR

FILES=($(ls $INPUTDIR/*.dose.vcf.gz))

CURRENT_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}

echo current dataset is $CURRENT_FILE

# replace with extension string
WORDTOREMOVE1="extension"
# replace with directory location
WORDTOREMOVE2="/home/"

OUTPUT_PHRASE="${CURRENT_FILE//$WORDTOREMOVE1/}"
OUTPUT_PHRASE="${OUTPUT_PHRASE//$WORDTOREMOVE2/}"

printf "Output phrase is $OUTPUT_PHRASE\n"

#######################################
# 1) Convert to bed and quality control
plink --vcf $CURRENT_FILE --update-ids Pair.txt --keep Keep.txt --maf 0.1 --hwe 1E-8 --geno 0.1 --mind 0.1 --snps-only just-acgt --make-bed --out ~/data/bedfiles/$OUTPUT_PHRASE

printf "\nDone converting to bed format.\n"

#######################################
# 2) Remove duplicate snps due to multiallelic variants
cut -f 2 ~/data/bedfiles/$OUTPUT_PHRASE.bim | sort | uniq -d > ~/data/bedfiles/$OUTPUT_PHRASE.dups
plink --bfile ~/data/bedfiles/$OUTPUT_PHRASE --exclude ~/data/bedfiles/$OUTPUT_PHRASE.dups --make-bed --out ~/data/bedfiles/filtered/filtered_$OUTPUT_PHRASE
rm ~/data/bedfiles/$OUTPUT_PHRASE.dups

printf "\nDone removing duplicates.\n"

printf '\nFinished!\n\n'
