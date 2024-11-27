import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import marsilea as ma
import marsilea.plotter as mp
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


def make_sim_mat(df, col_name):
    pivot_df = df.pivot_table(index='name_a', columns='name_b', values=col_name, fill_value=0)
    all_names = pd.unique(df[['name_a', 'name_b']].values.ravel('K'))
    pivot_df = pivot_df.reindex(index=all_names, columns=all_names, fill_value=0)
    similarity_matrix = pivot_df.copy()
    similarity_matrix = (similarity_matrix + similarity_matrix.T)
    np.fill_diagonal(similarity_matrix.values, 1)
    return similarity_matrix


def plot_heatmap(df, col, title, mthds, baselines, figs):
    order = mthds + baselines
    mat = make_sim_mat(df, col).loc[order, order]
    h = ma.Heatmap(mat, cmap='Purples', width=1.5, height=1.5, name=title, label="Overlap\nCoefficient", vmin=0, vmax=1)
    h.add_bottom(mp.Labels(mat.columns))
    h.add_left(mp.Labels(mat.index))
    h.add_top(mp.Title(title))
    h.add_legends()
    h.render()
    plt.close()
    figs.append(h.figure)


def barstats(df, col, title, figs):
    fig, axes = plt.subplots(2, 1, figsize=(1.7, 3), dpi=150, sharex=True, gridspec_kw={'height_ratios': [len(mthds), len(baselines)]})
    ax = axes[0]
    sns.barplot(
        data=df[df['name'].isin(mthds)],
        x=col,
        y='name',
        hue='name',
        orient='h',
        palette=palette,
        ax=ax,
        order=mthds
    )
    ax.set_ylabel('Methods')
    ax = axes[1]
    sns.barplot(
        data=df[df['name'].isin(baselines)],
        x=col,
        y='name',
        hue='name',
        orient='h',
        palette=palette,
        ax=ax,
        order=baselines
    )
    ax.set_ylabel('Baselines')
    ax.set_xlabel(title)
    fig.subplots_adjust(wspace=0.05, hspace=0)
    figs.append(fig)


sims = pd.read_csv(sys.argv[1])
stats = pd.read_csv(sys.argv[2])


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Format names
sims['name_a'] = [n.split('.')[0] for n in sims['name_a']]
sims['name_b'] = [n.split('.')[0] for n in sims['name_b']]
stats['name'] = [n.split('.')[0] for n in stats['name']]

# Filter for o_methods and baselines
msk_a = sims['name_a'].str.startswith('o_') | sims['name_a'].isin(baselines)
msk_b = sims['name_b'].str.startswith('o_') | sims['name_b'].isin(baselines)
sims = sims.loc[msk_a & msk_b, :]
msk = stats['name'].str.startswith('o_') | stats['name'].isin(baselines)
stats = stats.loc[msk, :]

# Remove o_
sims['name_a'] = sims['name_a'].str.replace('o_', '')
sims['name_b'] = sims['name_b'].str.replace('o_', '')
stats['name'] = stats['name'].str.replace('o_', '')

# Plot
figs = []

plot_heatmap(sims, col='tf_oc', title='TFs', mthds=mthds, baselines=baselines, figs=figs)
plot_heatmap(sims, col='edge_oc', title='Edges', mthds=mthds, baselines=baselines, figs=figs)
plot_heatmap(sims, col='target_oc', title='Genes', mthds=mthds, baselines=baselines, figs=figs)

barstats(stats, col='n_tfs', title='Number TFs', figs=figs)
barstats(stats, col='n_edges', title='Number Edges', figs=figs)
barstats(stats, col='n_targets', title='Number Genes', figs=figs)
barstats(stats, col='odegree', title='Regulon size', figs=figs)
barstats(stats, col='betweenc', title='B. centrality', figs=figs)
barstats(stats, col='eigv', title='E. centrality', figs=figs)

# Write
savefigs(figs, sys.argv[3])
