import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


# Read
df_knn = pd.read_csv(sys.argv[1]).sort_values('anchor', ascending=False)
df_knn = df_knn.groupby(['type', 'ctype', 'anchor'], as_index=False).mean(numeric_only=True)
df_cat = pd.read_csv(sys.argv[2])
df_cor = pd.read_csv(sys.argv[3]).sort_values('omic', ascending=False)
df_cor = df_cor.groupby(['type', 'ctype', 'omic'], as_index=False).mean(numeric_only=True)
df_oc = pd.read_csv(sys.argv[4])

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Plot
figs = []
fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150, tight_layout=True)
sns.boxplot(
    data=df_knn,
    x='anchor',
    y='k',
    hue='type',
    ax=ax,
    fill=None,
    fliersize=0,
)
sns.stripplot(
    data=df_knn,
    x='anchor',
    y='k',
    hue='type',
    ax=ax,
    dodge=True
)
ax.set_ylim(0, None)
ax.set_ylabel('K-Neighbor')
ax.set_xlabel('Anchor')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150, tight_layout=True)
mat_cat = df_cat.pivot(index='ctype_atac', columns='ctype_rna', values='prop')
sns.heatmap(mat_cat, cmap='Purples', annot=False, fmt='.2f', 
            square=True, ax=ax, cbar_kws={"shrink": 0.5}, vmin=0, vmax=1)
ax.set_xlabel('ATAC')
ax.set_ylabel('RNA')
ax.set_title('Proportion')
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150, tight_layout=True)
sns.boxplot(
    data=df_cor,
    x='omic',
    y='stat',
    hue='type',
    ax=ax,
    fill=None,
    fliersize=0,
)
sns.stripplot(
    data=df_cor,
    x='omic',
    y='stat',
    hue='type',
    ax=ax,
    dodge=True
)
ax.set_ylabel('Spearman\'s ρ')
ax.set_xlabel('')
ax.set_ylim(0, 1)
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

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

# Format names and filter for o_methods and baselines
df_oc['mth'] = [n.split('.')[0] for n in df_oc['mth']]
msk = df_oc['mth'].str.startswith('o_') | df_oc['mth'].isin(baselines)
df_oc = df_oc.loc[msk]
df_oc['mth'] = [m.replace('o_', '') for m in df_oc['mth']]
base_stability(df_oc, col='ocoef', mthds=mthds, baselines=baselines, palette=palette, figs=figs)

# Write
savefigs(figs, sys.argv[5])
