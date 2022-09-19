#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 4                      # Number of Tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --job-name=merge_chr_PCA

printf "jobstarting...\n\n"

module load plink

# Change your specific directories
INPUTDIR=~/data/inputdir
WORKDIR=~/data/workingdir
cd $WORKDIR

printf "working directory is $WORKDIR\n\n"

############################## Merge and QC #####################################
plink --merge-list all_my_files.txt --indiv-sort 0 --make-bed --out ~/data/workingdir/genome_temp

plink --bfile genome_temp --maf 0.1 --hwe 1E-8 --geno 0.1 --mind 0.1 --snps-only just-acgt --make-bed --out genome

rm *genome_temp*

############################## PCA #######################################
# Filter data (you can speed up things by adding --memory 119500 --threads 19 in PLINK)
plink --bfile genome --maf 0.01 --geno 0.01 --hwe 5e-6 --autosome --exclude exclusion_regions_hg19.txt --make-bed --out FILENAME_2
  
# Prune snps 
plink --bfile FILENAME_2 --indep-pairwise 1000 10 0.02 --autosome --out pruned_data

# Extract pruned SNPs and only these variants will be used for PC calculation
plink --bfile FILENAME_2 --extract pruned_data.prune.in --make-bed --out FILENAME_3
rm *FILENAME_2*
rm *pruned_data*

# Calculate/generate PCs based on pruned data set
plink --bfile FILENAME_3 --pca --out PCA
rm *FILENAME_3*
rm PCA.log
rm PCA.nosex

printf "\ndone!\n\n"
