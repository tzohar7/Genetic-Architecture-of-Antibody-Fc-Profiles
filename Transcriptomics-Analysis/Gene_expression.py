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
#
import matplotlib
matplotlib.use('Qt5Agg')

#######################################################################################################################
"""Importing"""

################################################################
# Set Figure Settings
################################################################
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

# Set PDF printing settings
out_pdf = r'Expression' + '.pdf'
pdf = matplotlib.backends.backend_pdf.PdfPages(out_pdf)

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
# Reduce gene set based on expression in B cells
################################################################

# METHOD IS TO LOOK FOR EXPESSION ACROSS CELL-TYPES AND DONORS
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

mygenes = [['DLC1'], ['RNF141'], ['SPRY2'], ['TRAV10', 'TRAV8-2', 'TRAV12-1'], ['FUT2'],
           ['IGHG1', 'IGHG2', 'IGHG3', 'IGHG4', 'IGHA1', 'IGHA2', 'IGHM']]

mygene_type = ['GWAS gene']*7 + ['IGH gene']*7

all_genes = gi + gn + np.sum(mygenes)
all_types = gi_type + gn_type + mygene_type

########################################################################################################################
########################################################################################################################
# Expression of Candidate Genes
########################################################################################################################
########################################################################################################################

col_mask = adata_orig.var.loc[:, 'feature_name'].isin(all_genes).tolist()
colnames = adata_orig.var.loc[:, 'feature_name'][col_mask].tolist()

dfs = []
sp_interaction = []
npdata = adata_norm.X[:, col_mask]
indexnames = adata_orig.obs.loc[:, 'cell_type'].tolist()
df = pd.DataFrame(index=indexnames, columns=colnames, data=npdata.toarray())
df = df.loc[:, np.sum(mygenes[:-1])]

flt = df.unstack().reset_index(drop=False)
flt.columns = ['Gene', 'Cells', 'Value']

cell_drop = 'muscle|eye|retina|tendon|pericyte|mesothelial|urothelial|adventitial|Schwann|Muller|pancreatic|sperm|glial|transit |' \
            'myometrial|ectodermal|platelet|fat|enteroendo|leuko'
flt = flt.loc[~flt.Cells.str.contains(cell_drop), :]

flt.loc[flt.Cells.str.contains('myeloid cell'), 'Cells'] = 'Myeloid cell'
flt.loc[flt.Cells.str.contains('common myeloid progenitor'), 'Cells'] = 'Common Myeloid Progenitor'

flt.loc[flt.Cells.str.contains('NK|Nk|killer'), 'Cells'] = 'NK cell'

flt.loc[flt.Cells.str.contains('plasmacytoid dendritic cell'), 'Cells'] = 'Plasmacytoid Dendritic cell'
flt.loc[flt.Cells.str.contains('CD141'), 'Cells'] = 'CD141+ Myeloid Dendritic cell'
flt.loc[flt.Cells.str.contains('CD1c'), 'Cells'] = 'CD1c+ Myeloid Dendritic cell'
flt.loc[flt.Cells == 'dendritic cell', 'Cells'] = 'Dendritic cell'

flt.loc[flt.Cells.str.contains('regulatory'), 'Cells'] = 'T reg'
flt.loc[flt.Cells.str.contains('helper'), 'Cells'] = 'Helper T cell'
flt.loc[flt.Cells.str.contains('CD4'), 'Cells'] = 'CD4+ T cell'
flt.loc[flt.Cells.str.contains('CD8'), 'Cells'] = 'CD8+ T cell'

flt.loc[flt.Cells.str.contains('plasma cell|plasmablast'), 'Cells'] = 'Plasma cell'
flt.loc[flt.Cells.str.contains('thymocyte|thymic'), 'Cells'] = 'Thymocyte'

flt.loc[flt.Cells.str.contains('ciliated'), 'Cells'] = 'Ciliated cell'

flt.loc[flt.Cells.str.contains('epithe|club|cholangiocyte'), 'Cells'] = 'Epithelial cell'
flt.loc[flt.Cells.str.contains('endothelial'), 'Cells'] = 'Endothelial cell'
flt.loc[flt.Cells.str.contains('fibroblast'), 'Cells'] = 'Fibroblast'
flt.loc[flt.Cells.str.contains('ionocyte'), 'Cells'] = 'Ionocyte'

flt.loc[flt.Cells.str.contains('tuft'), 'Cells'] = 'Tuft cell'
flt.loc[flt.Cells.str.contains('goblet'), 'Cells'] = 'Goblet cell'
flt.loc[flt.Cells.str.contains('stromal'), 'Cells'] = 'Stromal cell'
flt.loc[flt.Cells.str.contains('pneumocyte'), 'Cells'] = 'Pneumocyte'
flt.loc[flt.Cells.str.contains('erythro'), 'Cells'] = 'Erythrocyte'
flt.loc[flt.Cells.str.contains('intestinal crypt stem cell'), 'Cells'] = 'Intestinal crypt stem cell'
flt.loc[flt.Cells.str.contains('salivary gland'), 'Cells'] = 'Salivary gland cell'
flt.loc[flt.Cells.str.contains('kerat'), 'Cells'] = 'Keratocyte'
flt.loc[flt.Cells.str.contains('paneth'), 'Cells'] = 'Paneth cell'

######################################
# DOT PLOTS
######################################

dots = flt.groupby(['Cells', 'Gene']).apply(lambda x: (x > 0).sum()/len(x)*100)
dots = dots.reset_index(drop=False)
dots['Expression'] = flt.groupby(['Cells', 'Gene']).mean().reset_index(drop=False)['Value'].values
dots['Number_of_Cells'] = flt.groupby(['Cells', 'Gene']).apply(lambda x: (x > 0).sum()).reset_index(drop=False)['Value'].values
dots = dots.loc[~dots.Gene.isin(['TRAV10', 'TRAV12-1']), :]
dots.loc[dots['Number_of_Cells'] < 10, ['Expression', 'Value']] = [0, 0]

dots.loc[210:, 'Cells'] = dots.loc[210:, 'Cells'].str.capitalize()
dots.loc[dots.Cells == 'Memory b cell', 'Cells'] = 'Memory B cell'
dots.loc[dots.Cells == 'Naive b cell', 'Cells'] = 'Naive B cell'
dots = dots.loc[dots.Cells != 'Monocyte', :]
dots = dots.loc[dots.Cells != 'Langerhans cell', :]
dots = dots.loc[dots.Cells != 'Erythrocyte', :]
dots = dots.loc[dots.Cells != 'Myeloid dendritic cell', :]

sorter = ['Naive B cell', 'B cell', 'Memory B cell', 'Plasma ',
          'T cell', 'Helper', 'CD4', 'CD8', 'T reg', 'Thymocyte',
          'NK', 'Innate ',
          'Hematopoietic', 'Common', 'Myeloid', 'monocyte', 'Macrophage', 'Dendritic', 'Mast',
          'Gran', 'Neutrophil', 'Basophil',
          'Epith', 'Endo', 'Hepato', 'Mucus|Secretory', 'Pneumo|Salivary|Ion',
          'Ciliated|Duoden|Goblet|Tuft|Intest|Paneth']

dots['rank_ctype'] = len(sorter)+1
for i in range(len(sorter)):
    dots.loc[dots.Cells.str.contains(sorter[i]), 'rank_ctype'] = i

sorter = ['DLC1', 'RNF141', 'SPRY2', 'TRAV8-2', 'FUT2']

dots['rank_gtype'] = 0
for i in range(len(sorter)):
    dots.loc[dots.Gene.str.contains(sorter[i]), 'rank_gtype'] = i

dots = dots.sort_values(by=['rank_ctype', 'Cells', 'rank_gtype'], ascending=[True, True, True]).reset_index(drop=True)

###############################
# Plot
##############################

colormax = 3
colormin = 0
sizemax = 100
sizemin = 0

dotsi = dots.loc[:134, :]
x = dotsi['Gene']
y = dotsi['Cells']
colorz = dotsi['Expression']
colorz = (colorz - colormin) / (colormax - colormin)
size = dotsi['Value']
size_scaled = (size - sizemin) / (sizemax - sizemin)

n_colors = 256  # Use 256 colors for the diverging color palette
palette = sns.color_palette('magma_r', n_colors=n_colors)
# palette = sns.diverging_palette(20, 220, n=n_colors)  # Create the palette
# Range of values that will be mapped to the palette, i.e. min and max possible correlation
color_min, color_max = [0, 1]


def value_to_color(val):
    # position of value in the input range, relative to the length of the input range
    val_position = float((val - color_min)) / (color_max - color_min)
    ind = int(val_position * (n_colors - 1))  # target index in the color palette
    return palette[ind]

sns.set()
f = plt.figure(figsize=(3.5, 5), constrained_layout=True)

gs = f.add_gridspec(1, 15, hspace=0, wspace=0)
ax = f.add_subplot(gs[:, :])
ax.set_facecolor('w')

# Mapping from column names to integer coordinates
x_labels = [v for v in x.unique()]
y_labels = [v for v in y.unique()]
x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}

size_scale = 195
ax.scatter(
    x=x.map(x_to_num),  # Use mapping for x
    y=-y.map(y_to_num),  # Use mapping for y
    s=size_scaled * size_scale,  # Vector of square sizes, proportional to size parameter
    c=colorz.apply(value_to_color),
    marker='o',  # Use square as scatterplot marker
    linewidth=0.5,
    edgecolors='#6C6F70',
)

# Show column labels on the axes
ax.set_xticks([x_to_num[v] for v in x_labels])
ax.set_xticklabels(x_labels, rotation=45, ha="center", fontsize=10)

ax.set_yticks([-y_to_num[v] for v in y_labels])
ax.set_yticklabels(y_labels, ha="right", fontsize=10)

ax.grid(False, 'major')
ax.grid(True, 'minor')
ax.set_xticks([t + 0.5 for t in ax.get_xticks()], minor=True)
ax.set_yticks([t + 0.5 for t in ax.get_yticks()], minor=True)

ax.tick_params(axis='both', pad=0, colors='k', length=2, width=1, tickdir='out')
[q.set_linewidth(0.5) for q in ax.spines.values()]
[q.set_color('grey') for q in ax.spines.values()]

ax.set_xlim([-0.5, max([v for v in x_to_num.values()]) + 0.5])
ax.set_ylim([min([-v for v in y_to_num.values()]) - 0.5, 0.5])
pdf.savefig(f)

###############################
# Legend
###############################

# HEATMAP
f, ax = plt.subplots(figsize=(1.25, 0.90), constrained_layout=True)
cmap = sns.color_palette('magma_r', n_colors=n_colors, as_cmap=True)

norm = matplotlib.colors.Normalize(0, 2)
ax.grid(False)
# space = np.linspace(0, 2, 3, dtype=int)
space = [0, 1, 2]
cb = matplotlib.colorbar.ColorbarBase(ax, cmap=cmap, norm=norm, ticks=space,
                                      orientation='horizontal')

cb.set_label('Average Expression\n(log-CPM)', fontsize=10)
ax.set_xticks([0.10, 1, 1.90])
ax.set_xticklabels([0, 1, 2],)
ax.tick_params(axis='both', labelsize=10, pad=1)
ax.xaxis.set_ticks_position('top')
pdf.savefig(f)

# CIRCLES
dotsl = dotsi.loc[:5, :]
x = pd.Series(list(range(5)))
y = pd.Series([0]*5)
size = np.arange(20, 120, 20)
size_scaled = (size - sizemin) / (sizemax - sizemin)

sns.set_style("ticks")
f, ax = plt.subplots(figsize=(2.5, 1), constrained_layout=True)

ax.set_yticks([])

x_labels = [v for v in x.unique()]
y_labels = [v for v in y.unique()]
x_to_num = {p[1]: p[0] for p in enumerate(x_labels)}
y_to_num = {p[1]: p[0] for p in enumerate(y_labels)}

size_scale = 195
ax.scatter(
    x=x.map(x_to_num),  # Use mapping for x
    y=-y.map(y_to_num),  # Use mapping for y
    s=size_scaled * size_scale,  # Vector of square sizes, proportional to size parameter
    c=['k']*5,
    marker='o',  # Use square as scatterplot marker
    linewidth=0.5,
    edgecolors='#6C6F70',
)

ax.set_xticks(x)
ax.set_xticklabels(size.tolist(), fontsize=10)
ax.set_xlabel('Fraction of cells\nexpressing gene (%)', fontsize=10)
[q.set_linewidth(0) for q in ax.spines.values()]
pdf.savefig(f)

pdf.close()
