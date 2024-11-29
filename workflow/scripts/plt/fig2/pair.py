import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import mudata as mu
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scanpy as sc
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


pmdata = mu.read(sys.argv[1])
nmdata = mu.read(sys.argv[2])
df_ral = pd.read_csv(sys.argv[3])
df_qc = pd.read_csv(sys.argv[4])
df_nc = pd.read_csv(sys.argv[5])

palette = {
    'rna': 'purple',
    'atac': '#ffd700'
}

figs = []
fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
sc.pl.umap(
    pmdata,
    color='celltype',
    ax=ax,
    frameon=False,
    add_outline=True,
    title='',
)
figs.append(fig)
fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
sc.pl.umap(
    nmdata,
    color='celltype',
    ax=ax,
    frameon=False,
    add_outline=True,
    title='',
)
figs.append(fig)
fig, ax = plt.subplots(1, 1, figsize=(1.7, 3), dpi=150)
sns.barplot(
    data=df_nc,
    x='size',
    y='celltype',
    hue='type',
    ax=ax
)
ax.set_xlabel('Number of cells')
ax.set_ylabel('')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(1.7, 3), dpi=150)
sns.barplot(
    data=df_ral[df_ral['type'] == 'obs'],
    x='stat',
    y='name',
    hue='omic',
    palette=palette,
    ax=ax
)
ax.set_xlabel('Spearman\'s ρ')
ax.set_ylabel('')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), dpi=150, tight_layout=True)
sns.violinplot(
    data=df_qc,
    x='type',
    y='log1p_total_counts',
    hue='omic',
    density_norm='width',
    ax=ax,
    palette=palette
)
ax.set_xlabel('')
ax.set_ylabel('Total counts (log1p)')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(2.5, 2.5), dpi=150, tight_layout=True)
sns.violinplot(
    data=df_qc,
    x='type',
    y='log1p_n_genes_by_counts',
    hue='omic',
    density_norm='width',
    ax=ax,
    palette=palette
)
ax.set_xlabel('')
ax.set_ylabel('Total genes (log1p)')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

savefigs(figs, sys.argv[6])
