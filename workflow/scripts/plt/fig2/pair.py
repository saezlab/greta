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
df_oc = pd.read_csv(sys.argv[6])

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

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

cmap = {
    'rna': 'purple',
    'atac': '#ffd700'
}
fig, ax = plt.subplots(1, 1, figsize=(1.7, 3), dpi=150)
sns.barplot(
    data=df_ral[df_ral['type'] == 'obs'],
    x='stat',
    y='name',
    hue='omic',
    palette=cmap,
    ax=ax
)
ax.set_xlabel('Spearman\'s ρ')
ax.set_ylabel('')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(2, 1.7), dpi=150, tight_layout=True)
df_qc = df_qc.groupby(['celltype', 'omic', 'type'], as_index=False).mean(numeric_only=True)
sns.boxplot(
    data=df_qc,
    x='type',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    fill=None,
)
sns.stripplot(
    data=df_qc,
    x='type',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    dodge=True
)
ax.set_xlabel('')
ax.set_ylabel('Total counts (log1p)')
ax.legend().set_visible(False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(2, 1.7), dpi=150, tight_layout=True)
sns.boxplot(
    data=df_qc,
    x='type',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    fill=None,
)
sns.stripplot(
    data=df_qc,
    x='type',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    dodge=True
)
ax.set_xlabel('')
ax.set_ylabel('Total features (log1p)')
ax.legend().set_visible(False)
figs.append(fig)

# Format names and filter for o_methods and baselines
df_oc['mth'] = [n.split('.')[0] for n in df_oc['mth']]
msk = df_oc['mth'].str.startswith('o_') | df_oc['mth'].isin(baselines)
df_oc = df_oc.loc[msk]
df_oc['mth'] = [m.replace('o_', '') for m in df_oc['mth']]

def base_stability(df, col, mthds, baselines, palette, figs):
    fig, axes = plt.subplots(2, 1, figsize=(1.7, 3), dpi=150, gridspec_kw={'height_ratios': [len(mthds), len(baselines)]})
    ax = axes[0]
    sns.barplot(data=df[df['mth'].isin(mthds)], y='mth', x=col, hue='mth', ax=ax, palette=palette)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Methods')
    ax.set_xlim(-0.05, 1.05)
    
    ax = axes[1]
    sns.barplot(data=df[df['mth'].isin(baselines)], y='mth', x=col, hue='mth', ax=ax, palette=palette)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('Overlap coefficient')
    ax.set_ylabel('Baselines')
    ax.set_xlim(-0.05, 1.05)
    
    fig.subplots_adjust(wspace=0.05, hspace=0)
    figs.append(fig)

base_stability(df_oc, col='ocoef', mthds=mthds, baselines=baselines, palette=palette, figs=figs)

# Write
savefigs(figs, sys.argv[7])
