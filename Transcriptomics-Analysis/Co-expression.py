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
import gseapy as gp

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
# Co-expression with candidate gene and based on cell type
################################################################

# Example inputs
#celltypes = ['common myeloid progenitor', 'myeloid cell',
#               'CD141-positive myeloid dendritic cell', 'CD1c-positive myeloid dendritic cell', 'dendritic cell',
#               'monocyte', 'classical monocyte', 'intermediate monocyte', 'non-classical monocyte', 'macrophage',
#               'granulocyte', 'neutrophil']
#
#candidate_gene = 'RNF141'


def coexpression(celltypes, candidate_gene):

    Donors = adata_orig.obs['donor'].unique().tolist()
    Donors.remove('TSP12')

    gwas_gene = [candidate_gene]
    cellsi_list = celltypes

    colnames_forfilter = pd.DataFrame(adata_orig.var.loc[:, 'feature_name'].tolist())
    colnames_forfilter['coltype'] = colnames_forfilter[0].map(genesymbol_to_biotype)
    colnames_forfilter.columns = ['gene', 'coltype']
    colnames_forfilter = colnames_forfilter.loc[colnames_forfilter.coltype.notna(), :]
    colnames_forfilter = colnames_forfilter.loc[colnames_forfilter.coltype.str.contains('protein_coding'), :]  # Only protein coding
    colnames_forfilter = colnames_forfilter.loc[~colnames_forfilter.gene.str.contains('MT-'), :]  # No Mito genes

    # Prune genes with very poor expression across cells
    colnames_forfilter = colnames_forfilter.gene.tolist()

    ##########################################################################
    # Method: Donor & Cell Bulking
    ##########################################################################

    row_mask = adata_orig.obs['cell_type'].isin(cellsi_list).tolist()

    # Prune genes
    col_mask = adata_orig.var.loc[:, 'feature_name'].isin(colnames_forfilter).tolist()
    colnames = adata_orig.var.loc[:, 'feature_name'][col_mask].tolist()

    npdata = adata_norm.X[row_mask, :]
    npdata = npdata[:, col_mask]
    indexnames = adata_orig.obs.loc[row_mask, 'cell_type'].tolist()
    dfi = pd.DataFrame(columns=colnames, data=npdata.toarray(), index=indexnames)
    dfi = dfi.loc[dfi[gwas_gene[0]] > 0, :]
    select = ((dfi != 0).sum(axis=0) > dfi.shape[0]*0.5) | (dfi.columns == gwas_gene[0])
    dfi = dfi.loc[:, select]

    res = dfi.loc[:, dfi.columns != gwas_gene[0]].apply(lambda x: list(stats.spearmanr(dfi[gwas_gene[0]], x)))

    bcorr = res.copy()
    p_threshold = (0.05 / bcorr.shape[1])
    bcorr.loc[0, bcorr.loc[1, :] > p_threshold] = 0
    bcorr = bcorr.loc[0, :]
    bcorr = bcorr.loc[bcorr.notna()]
    bcorr = bcorr.loc[bcorr != 0]
    bcorr = bcorr[(bcorr < -0.30) | (bcorr > 0.30)]
    vals = bcorr.sort_values()
    output_corr = vals

    ##########################################################################
    # Plotting Enrichr Results
    ##########################################################################

    gene_list = [vals[vals > 0].index.tolist(), vals[vals < 0].index.tolist()]
    gene_list = [i for i in gene_list if any(i)]

    for i in gene_list:
        # Run enrichr
        # if you are only intrested in dataframe that enrichr returned, please set no_plot=True
        # list, dataframe, series inputs are supported
        enr = gp.enrichr(gene_list=i,
                         gene_sets=['GO_Biological_Process_2021'],
                         organism='Human',
                         no_plot=True,
                         cutoff=0.5)

        enr_results = enr.results
        enr_results = enr_results.iloc[:10, :]
        n_genes = enr_results['Genes'].str.split(';', expand=True).notna().sum(axis=1)

        sizemax = max(n_genes)
        sizemin = 0

        x = -np.log10(enr_results['Adjusted P-value'])
        x[x > 16] = 16
        y = enr_results['Term']

        colorz = enr_results['Combined Score']
        colorz = (colorz - colorz.min()) / (colorz.max() - colorz.min())
        size = n_genes
        size_scaled = (size - sizemin) / (sizemax - sizemin)
        # size_scaled = size_scaled + 0.01

        n_colors = 256  # Use 256 colors for the diverging color palette
        palette = sns.color_palette('viridis_r', n_colors=n_colors)
        color_min, color_max = [0, 1]

        def value_to_color(val):
            val_position = float((val - color_min)) / (color_max - color_min)
            ind = int(val_position * (n_colors - 1))
            return palette[ind]

        sns.set_style("whitegrid")
        f, ax = plt.subplots(figsize=(7, 2), dpi=100, constrained_layout=True)

        # Mapping from column names to integer coordinates
        x_labels = [v for v in x.unique()]
        y_labels = [v for v in y.unique()]
        x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
        y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}

        size_scale = 75
        ax.scatter(
            x=x,  # Use mapping for x
            y=-y.map(y_to_num),  # Use mapping for y
            s=size_scaled * size_scale,  # Vector of square sizes, proportional to size parameter
            c=colorz.apply(value_to_color),
            marker='o',  # Use square as scatterplot marker
            linewidth=0.5,
            edgecolors='#6C6F70',
        )

        ax.set_xlabel('-log10 q-value', fontsize=8)

        ax.set_yticks([-y_to_num[v] for v in y_labels])
        ax.set_yticklabels(y_labels, ha="right", fontsize=7)

        ax.tick_params(axis='both', pad=2, colors='k', length=2, width=1, tickdir='out')
        [q.set_linewidth(0.5) for q in ax.spines.values()]
        [q.set_color('grey') for q in ax.spines.values()]

        ax.set_ylim([min([-v for v in y_to_num.values()]) - 0.5, 0.5])

        prots = enr_results['Genes'].str.replace(';', '|')[0]
        print(vals.loc[vals.index.str.contains(prots)])
        print([dfi.groupby(dfi.index).mean()[i].idxmax() for i in enr_results['Genes'].str.split(';')[0]])

    return output_corr
