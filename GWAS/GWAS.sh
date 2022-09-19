#!/bin/bash
#SBATCH -N 1                      # Number of nodes. You must always set -N 1 unless you receive special instruction from the system admin
#SBATCH -n 16                     # Number of Tasks. Don't specify more than 16 unless approved by the system admin
#SBATCH --job-name=GWAS
#SBATCH --array=0-203             # Number of GWAS run in parallel 

printf "Program has begun...\n\n"

# Change to location of python app
~/data/software/python-3.8.12/bin/python3.8 fast-lmm_assoc.py ${SLURM_ARRAY_TASK_ID}

printf "\n\nPython Closed and Done!"

