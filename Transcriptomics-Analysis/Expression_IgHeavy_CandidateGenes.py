import pandas as pd
import os
import seaborn as sns
import numpy as np
from scipy import stats
import matplotlib.pyplot as plt
import matplotlib.backends.backend_pdf
import statsmodels.api as sm
import statsmodels.formula.api as smf
from matplotlib.lines import Line2D
import random
import scanpy as sc
import biomart
import anndata

import matplotlib
matplotlib.use('Qt5Agg')

#######################################################################################################################
#######################################################################################################################
"""Importing"""

################################################################
# Set Figure Settings
################################################################
# Close all figures
plt.close('all')

sns.set(font="Arial")
sns.set_context("paper")
user_dpi = 100
matplotlib.rcParams['pdf.fonttype'] = 42
matplotlib.rcParams['ps.fonttype'] = 42

# Panda Display Settings
pd.set_option('display.max_rows', 500)
pd.set_option('display.max_columns', 500)
pd.set_option('display.width', 1000)
pd.set_option('display.float_format', lambda x: '%.2f' % x)

################################################################
# Load Data
################################################################
# Data needs to be downloaded from here: https://cellxgene.cziscience.com/collections/e5f58829-1a66-40b5-a624-9046778e74f5
fn = "local.h5ad"
adata_orig = sc.read_h5ad(fn)

adata_norm = adata_orig.raw.to_adata()
sc.pp.normalize_total(adata_norm, target_sum=1e6, inplace=True)
sc.pp.log1p(adata_norm)

################################################################
# Get Gene Type
################################################################
# Set up connection to server
server = biomart.BiomartServer('http://uswest.ensembl.org/biomart')
mart = server.datasets['hsapiens_gene_ensembl']

# List the types of data we want
attributes = ['gene_biotype', 'hgnc_symbol']

# Get the mapping between the attributes
response = mart.search({'attributes': attributes})
data = response.raw.data.decode('ascii')

genesymbol_to_biotype = {}
# Store the data in a dict
for line in data.splitlines():
    line = line.split('\t')
    # The entries are in the same order as in the `attributes` variable
    ensembl_gene = line[0]
    gene_symbol = line[1]

    genesymbol_to_biotype[gene_symbol] = ensembl_gene

################################################################
# My gene processing and organization
################################################################
genesin = pd.read_csv("rnaseq_genes.txt", delim_whitespace=True, header=None)
genesin.index = genesin[1]
genesin.index.name = None
genesin = genesin.iloc[:, 2:]
genesin = genesin[~genesin.index.str.contains('IgE|IgD')]
glyco_genes = pd.read_csv("glyco_genes.csv", sep=',')

dfiso = genesin[genesin.index.str.contains('isotype')]

set_list = []
set_names = []

for i in range(dfiso.shape[0]):
    si = dfiso.iloc[i, :]
    si.loc[si.notna()].tolist()
    set_list.append(set(si.loc[si.notna()]))
    set_names.append(si.name)

set_list = set_list[::-1]
set_names = set_names[::-1]

# LOWEST LEVEL
posigg = list(set_list[0].difference(set_list[2], set_list[3]))
posiga = list(set_list[2].difference(set_list[0], set_list[3]))
negiga = list(set_list[3].difference(set_list[0], set_list[2]))

# LEVEL REG ISOTYPE SPECIFIC
regigg = list(set_list[1].difference(set_list[0], set_list[2], set_list[3], set_list[4]))
regiga = list(set_list[4].difference(set_list[0], set_list[2], set_list[3], set_list[1]))

# LEVEL ISOTYPE SPECIFIC
isoigg = list(set_list[5].difference(*set_list[:4], set_list[6]))
isoiga = list(set_list[6].difference(*set_list[:4], set_list[5]))

# LEVEL SIGNED REG
pos = list(set_list[7].difference(*set_list[:6], set_list[8])) + list(set_list[0].intersection(set_list[2]))
neg = list(set_list[8].difference(*set_list[:6], set_list[7]))

# LEVEL REG
reg = list(set_list[9].difference(*set_list[:9]))

# LEVEL FINAL
iso = list(set_list[10].difference(*set_list[:10]))

all_levels = [iso, reg, neg, pos, isoiga, isoigg, regiga, negiga, posiga, regigg, posigg]
idx = set_names[::-1]
all_names = [[idx[i]]*len(all_levels[i]) for i in range(len(all_levels))]

for i in range(len(all_levels)):
    genesin.loc[genesin.index == idx[i],
                (~genesin.loc[genesin.index == idx[i], :].isin(all_levels[i])).values[0].tolist()] = np.nan

all_levels = [sorted(i) for i in all_levels]

gi = np.sum(all_levels)
gi_lut = dict(zip(gi, np.sum(all_names)))

gn = glyco_genes['Name'].tolist()
gn_lut = dict(zip(glyco_genes['Name'].tolist(), glyco_genes['Group'].tolist()))

################################################################
# Redefine genes
################################################################
csr_genes = pd.read_csv("CSR_genes.csv", sep=',')
csr_genes = csr_genes.set_index(csr_genes['Process']).drop(columns='Process')

csr_genes_all = csr_genes.unstack().reset_index()[0].unique().tolist()
csr_genes_all.remove(np.nan)

set_list = []
set_names = []

for i in range(csr_genes.shape[0]):
    si = csr_genes.iloc[i, :]
    set_list.append(set(si.loc[si.notna()]))
    set_names.append(si.name)

gi = [j for i in set_list for j in i]
gi = np.unique(gi).tolist()

dflut = csr_genes.unstack().reset_index().drop(columns='level_0')
dflut.columns = ["Process", "Gene"]
dflut = dflut.loc[dflut["Gene"].notna(),:]
dflut = dflut.sort_values('Gene')

################################################################
# Reduce gene set based on expression in plasma cells
################################################################

col_mask = adata_orig.var.loc[:, 'feature_name'].isin(gi + gn).tolist()
npdata = adata_norm.X[:, col_mask].toarray()
colnames = adata_orig.var.loc[:, 'feature_name'][col_mask].tolist()
indexnames = adata_orig.obs.loc[:, 'cell_type'].tolist()
df = pd.DataFrame(index=indexnames, columns=colnames, data=npdata)
df = df.loc[:, gi]
df = df.loc[df.index == 'B cell']

cols_keep = (df.apply(lambda x: (x > 0).sum()) > 0)
cols_keep = cols_keep[cols_keep].index.tolist()

gi = [i for i in gi if i in cols_keep]
gn = ['MGAT3', 'B4GALT1', 'ST6GAL1', 'FUT8']

gi_type = pd.Series(gi).map(gi_lut).tolist()
gn_type = pd.Series(gn).map(gn_lut).tolist()

mygenes = [['DLC1'], ['RNF141', 'IRAG1', 'IRAG1-AS1'], ['SPRY2', 'PTMAP5'], ['TRAV7', 'TRAV12-2', 'TRAV13-1'], ['FUT2'],
           ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHA1', 'IGHA2', 'IGHM']]

mygene_type = ['GWAS gene']*10 + ['IGH gene']*7

all_genes = gi + gn + np.sum(mygenes)
all_types = gi_type + gn_type + mygene_type

################################################################
# IgHeavy Chain Function
################################################################

# Example inputs
#celltypes = ['common myeloid progenitor', 'myeloid cell',
#             'monocyte', 'classical monocyte', 'intermediate monocyte', 'non-classical monocyte', 'macrophage']
#candidate_gene = 'DLC1'

def igheavv_corr(celltypes, candidate_gene):

    Donors = adata_orig.obs['donor'].unique().tolist()
    Donors.remove('TSP12')

    gwas_gene = [candidate_gene]
    gene_test = all_genes[-7:]
    gene_coc = gwas_gene + gene_test

    dons = np.repeat(Donors, len(gwas_gene)).tolist()
    gwas_gene_rep = np.tile(gwas_gene, len(Donors)).tolist()
    ind_names = [dons[i] + '_' + gwas_gene_rep[i] for i in range(len(dons))]
    ind_names = [[i + '_bcell', i + '_xcell'] for i in ind_names]
    ind_names = np.array(ind_names).ravel().tolist()

    col_mask = adata_orig.var.loc[:, 'feature_name'].isin(gene_coc).tolist()
    colnames = adata_orig.var.loc[:, 'feature_name'][col_mask].tolist()

    # Takes in multiple cell types or only one
    if isinstance(celltypes, list):
        cellsi = celltypes
    else:
        cellsi = [celltypes]

    base_cell = ['plasma cell', 'plasmablast']
    df_save = pd.DataFrame(index=ind_names, columns=colnames)

    place = 0

    for w in range(len(Donors)):

        row_mask = (adata_orig.obs['donor'] == Donors[w]).tolist()

        npdata = adata_norm.X[row_mask, :]
        npdata = npdata[:, col_mask].toarray()

        indexnames = adata_orig.obs.loc[row_mask, 'cell_type'].tolist()

        df = pd.DataFrame(index=indexnames, columns=colnames, data=npdata)

        df = df.reset_index()
        df = df.rename(columns={'index': 'Cells'})

        # Expression of Ig Heavy Chain in Base cell and Candidate Gene in other cell-type
        for k in range(2):

            # Expression of Ig Heavy Chain in Base cell
            if k == 0:
                dfb = df.loc[df.Cells.isin(base_cell), :]
                if not dfb.empty:
                    df_save.iloc[place, :] = dfb.mean()

                    total_ig = dfb.loc[:, dfb.columns.isin(mygenes[-1])].sum(axis=1).sum()

                    for g in mygenes[-1]:
                        df_save.iloc[place, df_save.columns.get_loc(g)] = (dfb.loc[:, g].sum()/total_ig)*100

                place = place + 1

            # Expression of Candidate Gene in other cell-type
            if k == 1:
                dfx = df.loc[df.Cells.isin(cellsi), gwas_gene[0]]
                print(df.loc[df.Cells.isin(cellsi), 'Cells'].unique().tolist())
                df_save.iloc[place, :] = 0
                df_save.iloc[place, df_save.columns.get_loc(gwas_gene[0])] = dfx[dfx>0].mean()

                place = place + 1

        if ((dfb.shape[0] < 10) | (dfx[dfx > 0].shape[0] < 10)):
            df_save.iloc[place-1, :] = np.nan
            df_save.iloc[place-2, :] = np.nan

    df_save = df_save.reset_index()
    df_save = df_save.loc[~df_save.isna().any(axis=1), :]

    df_save = df_save.rename(columns={'index': 'inds'})

    df_save['Donor'] = df_save.inds.str.split('_', expand=True)[0]
    df_save['GWAS_Gene'] = df_save.inds.str.split('_', expand=True)[1]
    df_save['Cell_Axis'] = df_save.inds.str.split('_', expand=True)[2]

    df_corr = pd.DataFrame(index=gwas_gene, columns=colnames)
    dfsig = df_corr.copy()
    dfis = []

    for i in range(df_corr.shape[0]):
        dfi = df_save.loc[df_save.GWAS_Gene == gwas_gene[i], :].reset_index(drop=True)
        dfi = pd.concat([dfi.loc[dfi.Cell_Axis == 'xcell', dfi.GWAS_Gene[0]].reset_index(drop=True),
                         dfi.loc[dfi.Cell_Axis == 'bcell', dfi.columns != dfi.GWAS_Gene[0]].reset_index(drop=True)],
                        axis=1)
        dfi = dfi.drop(columns=['inds', 'GWAS_Gene', 'Cell_Axis', 'Donor']).astype(float)

        dfi = dfi.loc[~dfi.isna().all(axis=1), :]
        dfi = dfi.loc[~dfi.isna().any(axis=1), :]

        dfis.append(dfi)

        # Calculate correlations
        dfcorr = dfi.corr(method='spearman')
        coef, p = stats.spearmanr(dfi)
        ps = p[0, :]

        dfcorr = dfcorr.iloc[0, :]

        idx = [dfcorr.index.get_loc(i) for i in df_corr.columns]

        df_corr.iloc[i, :] = dfcorr.loc[df_corr.columns]
        dfsig.iloc[i, :] = ps[idx]

    df_corr = df_corr.astype(float)
    df_corr = df_corr.loc[:, ]

    # Plot Heatmap
    sns.set_style("white")
    f, axs = plt.subplots(figsize=(3, 1), dpi=100, constrained_layout=True)

    sns.heatmap(df_corr.loc[:, gene_test], cmap='RdBu_r', ax=axs, robust=False, vmin=-1, vmax=1, square=False,
                linecolor='#DCE2F3', linewidth=0.5,
                fmt='',
                annot_kws={"size": 7, "ha": "center", "va": "center"},
                cbar_kws={"shrink": 1})

    axs.set_title(None, fontsize=10)
    axs.tick_params(axis='both', pad=0, length=0)
