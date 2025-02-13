import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp
from matplotlib.patches import Rectangle
import matplotlib.pyplot as plt
import pandas as pd
import mudata as mu
import scanpy as sc
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def fixed_pip(mthds, res, title):
    steps = ['pre', 'c2g', 'tfb', 'mdl']
    fig, axes = plt.subplots(len(mthds), 1, sharex=True, sharey=True, figsize=(2, 1.5 * len(mthds)), dpi=150, tight_layout=True)
    for i, mth in enumerate(mthds):
        df = res[res['mth'] == mth]
        ax = axes[i]
        if i == 0:
            ax.set_title(title)
        sns.boxplot(data=df, x='step', order=steps, y=title, ax=ax, color=palette[mth], fill=None, fliersize=0)
        sns.stripplot(data=df, x='step', order=steps, y=title, ax=ax, color=palette[mth])
        ax.set_ylabel(mth)
        ax.set_ylim(-0.05, 1)
        ax.axhline(y=0.5, ls='--', color='black')
    return fig


def sim_mat(mat, sts, palette):
    h1 = ma.Heatmap(mat, cmap='Purples', width=3, height=3, name="h1", vmax=1, label="Overlap\nCoefficient", cbar_kws=dict(orientation='horizontal'))
    legend=False
    for cat_col in ['mdl', 'tfb', 'c2g', 'pre']:
        cats = sts[cat_col].to_list()
        cat_colors = mp.Colors(cats, palette={k: v for k, v in palette.items() if k in cats}, label=cat_col, legend_kws=dict(title=''))
        if cat_col == 'pre':
            legend=True
        h1.add_top(cat_colors, pad=0, size=0.15, legend=legend)
    h1.add_dendrogram("left", show=False)
    h1.add_dendrogram("top", show=False)
    h1.add_legends(
        side="right",
        stack_by='col',
        stack_size=2,
    )
    h1.render()
    plt.close()
    hax = h1.get_ax(name='h1')
    border = Rectangle((0, 0), 1, 1, fill=False, ec=".1", lw=2, transform=hax.transAxes)
    hax.add_artist(border)
    fig = h1.figure
    fig.set_dpi(150)
    return fig


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']


mdata = mu.read(sys.argv[1])
df_nc = pd.read_csv(sys.argv[2])
df_qc = pd.read_csv(sys.argv[3])
df = pd.read_csv(sys.argv[4])
sts = pd.read_csv(sys.argv[5])
fvsd = pd.read_csv(sys.argv[6])
df_stab = pd.read_csv(sys.argv[7])

# Remove original runs and baselines
df = df[~(df['name_a'].str.startswith('o_') | df['name_b'].str.startswith('o_'))]
df = df[~(df['name_a'].str.split('.', expand=True)[0].isin(baselines) | df['name_b'].str.split('.', expand=True)[0].isin(baselines))]

figs = []
# Add qc plots
fig, ax = plt.subplots(1, 1, figsize=(4, 4), dpi=150)
sc.pl.umap(
    mdata,
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
    ax=ax
)
ax.set_xlabel('Number of cells')
ax.set_ylabel('')
ax.legend().set_visible(False)
ax.legend(loc='lower left', bbox_to_anchor=(0.25, -0.5), frameon=False)
figs.append(fig)
fig, ax = plt.subplots(1, 1, figsize=(2, 1.7), dpi=150, tight_layout=True)
df_qc = df_qc.groupby(['celltype', 'omic'], as_index=False).mean(numeric_only=True)
cmap = {
    'rna': 'purple',
    'atac': '#ffd700'
}
sns.boxplot(
    data=df_qc,
    x='omic',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    fill=None,
    fliersize=0,
)
sns.stripplot(
    data=df_qc,
    x='omic',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
)
ax.set_xlabel('')
ax.set_ylabel('Total counts (log1p)')
ax.legend().set_visible(False)
figs.append(fig)

fig, ax = plt.subplots(1, 1, figsize=(2, 1.7), dpi=150, tight_layout=True)
sns.boxplot(
    data=df_qc,
    x='omic',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
    fliersize=0,
    fill=None,
)
sns.stripplot(
    data=df_qc,
    x='omic',
    y='log1p_n_genes_by_counts',
    hue='omic',
    ax=ax,
    palette=cmap,
)
ax.set_xlabel('')
ax.set_ylabel('Total features (log1p)')
ax.legend().set_visible(False)
figs.append(fig)

for oc in ['tf_oc', 'edge_oc', 'target_oc']:
    mat = df.dropna().pivot(index='name_a', columns='name_b', values=oc).fillna(0)
    mat = mat + mat.T
    np.fill_diagonal(mat.values, 1)
    t_sts = sts.set_index('name').loc[mat.index].rename(columns={'p2g': 'c2g'})
    t_sts[['pre', 'c2g', 'tfb', 'mdl']] = t_sts.reset_index()['name_a'].str.split('.', n=4, expand=True).values
    figs.append(fixed_pip(mthds, fvsd, title=oc))
    figs.append(sim_mat(mat, t_sts, palette))

fig, axes = plt.subplots(2, 1, figsize=(2, 4), tight_layout=True)
ax = axes[0]
sns.barplot(
    data=df_stab,
    y='mth',
    x='ocoeff',
    hue='mth',
    ax=ax,
    palette=palette
)
ax.set_xticks([0, 0.5, 1])
ax.set_xlabel('Edge Overlap\nCoefficient')
ax.set_ylabel('')
ax = axes[1]
sns.barplot(
    data=df_stab,
    y='mth',
    x='stat',
    hue='mth',
    ax=ax,
    palette=palette
)
ax.set_xlabel('Pearson ρ')
ax.set_xticks([0, 0.5, 1])
ax.set_ylabel('')
figs.append(fig)

# Write
savefigs(figs, sys.argv[8], index_pngs=[0, 5, 7, 9])
