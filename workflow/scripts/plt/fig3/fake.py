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
df_cat = pd.read_csv(sys.argv[2])
df_cor = pd.read_csv(sys.argv[3]).sort_values('omic', ascending=False)

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Plot
figs = []
fig, ax = plt.subplots(1, 1, figsize=(3, 3), dpi=150, tight_layout=True)
sns.violinplot(
    data=df_knn,
    x='anchor',
    y='k',
    hue='type',
    density_norm='width',
    ax=ax,
)
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
sns.violinplot(
    data=df_cor,
    x='omic',
    y='stat',
    hue='type',
    density_norm='width',
    ax=ax,
)
ax.set_ylabel('Spearman\'s ρ')
ax.set_xlabel('')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)

# Write
savefigs(figs, sys.argv[4])
