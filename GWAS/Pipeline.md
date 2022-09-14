# Genetic analysis pipeline
\
Preprocessing and GWAS
1. Preprocessing.sh
- Keep.txt file was made to inform which samples are specific to our study. (Must be requested with genotypic data)
- Pair.txt file was made to update family information used in plink software. (Must be requested with genotypic data)  
- VCF chromosome files are converted to BED files and qc’ed in the process with the Preprocessing.sh file
  
2. Merge_QC_PCA.sh
- All BED chromosome files were merged, and population stratification was achieved using the merge_qc_PCA.sh file
- Covariates.txt file is made for processing with fast-lmm and includes Age and first 10 PCs of population stratification. (Must be requested with genotypic data)

3. GWAS.sh
- Run GWAS with slurm file GWAS.sh that uses the fast-lmm_assoc.py file

\
Post-GWAS
1. Clumping.sh
- Clumping with post-GWAS_array_clump.sh

2. Make results tables for both lead SNPs and LD SNPs
- Run post-GWAS-analysis to generate results table
- Get rsIDs
- Get RegulomeDB ranks  
- Annotation annovar
- Get gene info

3. MAGMA_GSA.sh
- Run gene-set analysis with MAGMA_GSA.sh to get genes and gene sets for all features