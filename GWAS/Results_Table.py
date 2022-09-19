import os
import seaborn as sns
import numpy as np
from scipy import stats
import pandas as pd
import timeit
import sys
import re

import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
from matplotlib.lines import Line2D
import matplotlib.ticker as plticker
import matplotlib.ticker as mticker
from matplotlib.colors import Normalize

from fastlmm.association import single_snp_scale, single_snp
from pysnptools.util.mapreduce1.runner import LocalMultiProc
from pysnptools.util.pheno import loadOnePhen
import statsmodels.stats.multitest as multitest

# Panda Display Settings
pd.set_option('display.max_rows', 600)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', lambda x: '%.2f' % x)

################################################################
# FOR CLUMPED FILES
files = os.listdir("~data/clumped/")

# Only '*.txt' files
files = [x for x in files if ".clumped" in x]

os.chdir("~data/clumped/")

dflist = []
for i in range(len(files)):
	df = pd.read_csv(files[i], delim_whitespace=True)
	df = df.loc[:, ['CHR', 'SNP', 'BP', 'P']]
	df['Phenotype'] = files[i].split('.', 2)[0]
	dflist.append(df)

dfu = pd.concat(dflist)
dfu = dfu.sort_values(by=['P'])
dfl = dfu.loc[dfu.P < 5e-8]
dfl['P'] = -np.log10(dfl.P.values)

dfl['Feature'] = dfl.Phenotype.str.split('_',expand=True)[0]
dfl['Antigen'] = dfl.Phenotype.str.split('_',expand=True)[1]

dfl = dfl.reset_index(drop=True)

dfl = dfl.loc[dfl.Phenotype.isin(filter_pheno), :]

len(dfl['Phenotype'].unique())

dfl = dfl.sort_values(by=['CHR', 'BP', 'P'], ascending=[True, True, False])

#############################
# Getting the Alelle Information
#############################

bim = pd.read_csv("~data/bedfiles/filtered/genome.bim", delim_whitespace=True, header=None, usecols=[1, 4, 5])

bim.columns = ["SNP", "A1", "A2"]
bim["EA"] = bim.A1
bim["NEA"] = bim.A2
bim.index = bim.SNP
bim.index.name = None
bim = bim.loc[:, ["EA", "NEA"]]

# Storage Dataframe
dfl.insert(4,"BETA", 0)
dfl.insert(5,"SE", 0)
dfl.insert(6,"%Var",0)
dfl["EA"] = bim.loc[dfl.SNP,"EA"].values
dfl["NEA"] = bim.loc[dfl.SNP,"NEA"].values

#############################
# Get GWAS stastistics
#############################

allfiles = os.listdir("~data/results")

os.chdir("~data/results")
files = [i+'.txt' for i in filter_pheno]
files = [s for s in allfiles if any(xs in s for xs in files)]

for i in range(len(files)):
	df = pd.read_csv(files[i], sep="\t", nrows=1e4)
	p = files[i].split('.', 2)[0]
	df["Phenotype"] = p
	df = df.sort_values(by=['Chr', 'ChrPos', 'PValue'], ascending=[True, True, False])
	snpsi = dfl.loc[dfl.Phenotype == p, "SNP"].values
	dfl.loc[((dfl.Phenotype == p) & (dfl.SNP.isin(snpsi))), ["BETA", "SE", "%Var"]] = df.loc[df.SNP.isin(snpsi), ["SnpWeight", "SnpWeightSE", "SnpFractVarExpl"]].values

dfl.to_csv("~data/results/Results_table.txt", sep="\t", index=False)
