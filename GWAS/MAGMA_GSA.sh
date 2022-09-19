#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 1                      # Number of Tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --job-name=magma_snp
#SBATCH --array=1-203 

###########################################################################################
# Setup and Parallelization ###############################################################
###########################################################################################

printf "Job starting...\n\n"

genome="~data/bedfiles/filtered/genome"

cd ~/data/software/magma

INPUTDIR=~data/inputdir
FILES=($(ls $INPUTDIR/*.txt))

CURRENT_FILE=${FILES[${SLURM_ARRAY_TASK_ID}]}

WORDTOREMOVE="inputdir-location-name"

PHENO="${CURRENT_FILE//".txt"/}"
PHENO="${PHENO//$WORDTOREMOVE/}"

###########################################################################################
# Annotation ##############################################################################
###########################################################################################

# Annotate or "map" snps to genes for next process
# This can be done before
./magma --annotate window=10 --snp-loc $genome.bim --gene-loc humandb/NCBI37.3.gene.loc --out anno
rm anno.log

###########################################################################################
# Gene GWAS ###############################################################################
###########################################################################################

# For all SNPs
awk -v OFS="\t" '{print $2,$5}' $CURRENT_FILE | sed '1d' > $INPUTDIR/$PHENO.tmp

gene_output="~data/magma_snpGWAS/$PHENO"

./magma --bfile $genome --gene-annot anno.genes.annot --gene-model snp-wise=top,1 --pval $INPUTDIR/$PHENO.tmp  N=498 --out $gene_output

rm $gene_output.log
rm $INPUTDIR/$PHENO.tmp

###########################################################################################
# Gene Set Analysis #######################################################################
###########################################################################################

i=gene_sets/GO.entrez.gmt

prefix="${i//.entrez.gmt/}"
prefix="${prefix//gene_sets/}"
prefix="${prefix//\//}"
geneset_output="~data/magma/snpGWAS/${PHENO}_${prefix}"
./magma --gene-results $gene_output.genes.raw --set-annot $i --out $geneset_output
rm $geneset_output.log
