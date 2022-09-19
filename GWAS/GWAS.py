import os
import seaborn as sns
import numpy as np
from scipy import stats
import pandas as pd
import timeit
import sys

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.lines import Line2D
import matplotlib.ticker as plticker
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize

from fastlmm.association import single_snp_scale, single_snp
from pysnptools.util.mapreduce1.runner import LocalMultiProc
from pysnptools.util.pheno import loadOnePhen
from pysnptools.snpreader import Bed,Pheno

# Panda Display Settings
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', lambda x: '%.2f' % x)

################################################################
# Set up data
################################################################
# Location and name of files
bed_fn = "~data/bedfiles/filtered/genome.bed"
pheno_fn = "~data/bedfiles/filtered/Phenotypes.txt"
cov_fn = "~data/bedfiles/filtered/Covariates.txt"

# Input for each parallel  process
pheno_number = int(sys.argv[1])
pheno1 = loadOnePhen(pheno_fn, i_pheno=pheno_number, vectorize=False)

# Reorder the samples from the files in the same order
bed_data = Bed(bed_fn)
cov_data = Pheno(cov_fn).read()

bed_ids = bed_data.iid[:, 1].tolist()
cov_ids = cov_data.iid[:, 1].tolist()
pheno_ids = pheno1['iid'][:, 1].tolist()

idxp = [pheno_ids.index(i) for i in bed_ids]
pheno1['iid'] = pheno1['iid'][idxp, :]
pheno1['iid'][:, 0] = [i.split('.')[0] for i in pheno1['iid'][:, 0]]
pheno1['vals'] = pheno1['vals'][idxp]

idxc = [cov_ids.index(i) for i in bed_ids]
cov_data._row = cov_data._row[idxc, :]
cov_data._val = cov_data._val[idxc, :]

phenoname = pheno1['header'][0]
output_filename = "~/data/outputdir/" + phenoname + ".txt"
print('\n\n', phenoname)

################################################################
# Run GWAS one phenotype at a time
################################################################

print('\n\nStarting GWAS...\n\n')
start = timeit.default_timer()
from pysnptools.util.mapreduce1.runner import LocalMultiProc
runner = LocalMultiProc(taskcount=16) # Run on 16 additional Python processes
single_snp(test_snps=bed_fn, pheno=pheno1, covar=cov_data, count_A1=False, runner=runner, output_file_name=output_filename)
stop = timeit.default_timer()
print('\n\nRuntime: ', stop - start)
print('\n\nDone!')
