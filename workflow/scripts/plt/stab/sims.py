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
import argparse


def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


def make_sim_mat(df, col_name, prefix):
    pivot_df = df.pivot_table(index=f'{prefix}_a', columns=f'{prefix}_b', values=col_name, fill_value=0)
    all_names = pd.unique(df[[f'{prefix}_a', f'{prefix}_b']].values.ravel('K'))
    pivot_df = pivot_df.reindex(index=all_names, columns=all_names, fill_value=0)
    similarity_matrix = pivot_df.copy()
    similarity_matrix = (similarity_matrix + similarity_matrix.T)
    np.fill_diagonal(similarity_matrix.values, 1)
    return similarity_matrix


def plot_heatmap(df, col, title, prefix, mthds, baselines, figs, method_names):
    order = mthds + baselines
    mat = make_sim_mat(df, col, prefix).loc[order, order]
    h = ma.Heatmap(mat, cmap='Purples', width=3, height=3, name=title, label="Overlap\nCoefficient", vmin=0, vmax=1)
    h.add_bottom(mp.Labels(prettify(mat.columns, method_names)))
    h.add_left(mp.Labels(prettify(mat.index, method_names)))
    h.add_top(mp.Title(title))
    h.add_legends()
    h.render()
    plt.close()
    figs.append(h.figure)


def barstats(df, col, title, figs, method_names):
    # Create a copy with pretty names
    df_plot = df.copy()
    df_plot['name'] = prettify(df_plot['name'].tolist(), method_names)
    pretty_mthds = prettify(mthds, method_names)
    pretty_baselines = prettify(baselines, method_names)
    # Create palette with pretty names
    pretty_palette = {method_names.get(k, k): v for k, v in palette.items()}

    fig, axes = plt.subplots(2, 1, figsize=(1.7, 3), dpi=150, sharex=True, gridspec_kw={'height_ratios': [len(mthds), len(baselines)]})
    ax = axes[0]
    sns.barplot(
        data=df_plot[df_plot['name'].isin(pretty_mthds)],
        x=col,
        y='name',
        hue='name',
        orient='h',
        palette=pretty_palette,
        ax=ax,
        order=pretty_mthds
    )
    ax.set_ylabel('Methods')
    ax = axes[1]
    sns.barplot(
        data=df_plot[df_plot['name'].isin(pretty_baselines)],
        x=col,
        y='name',
        hue='name',
        orient='h',
        palette=pretty_palette,
        ax=ax,
        order=pretty_baselines
    )
    ax.set_ylabel('Baselines')
    ax.set_xlabel(title)
    fig.subplots_adjust(wspace=0.05, hspace=0)
    figs.append(fig)


parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_sims', required=True)
parser.add_argument('-b','--path_stats', required=True)
parser.add_argument('-c','--path_tss', required=True)
parser.add_argument('-d','--path_dts', required=True)
parser.add_argument('-e','--path_net', required=True)
parser.add_argument('-f','--path_out', required=True)
parser.add_argument('-g','--baselines', required=True, nargs='+')
args = parser.parse_args()

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = args.baselines
mthds = [m for m in mthds if m not in baselines]

# Method names mapping for pretty labels
method_names = config['method_names']
method_names['promoters'] = 'Promoters'  # Special case for TSS plot


def prettify(names, mapping):
    """Convert lowercase names to pretty names."""
    return [mapping.get(n, n) for n in names]

# Read
sims = pd.read_csv(args.path_sims)
stats = pd.read_csv(args.path_stats)
tss = pd.read_csv(args.path_tss)
dst = pd.read_csv(args.path_dts)
net = pd.read_csv(args.path_net)

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

# Format TSS
tss['tss_a'] = tss['tss_a'].str.replace('.gz', '')
tss['tss_b'] = tss['tss_b'].str.replace('.gz', '')
tss = tss.groupby(['tss_a', 'tss_b'], as_index=False)['ocoef'].mean()

# Format distances
dst = dst.assign(kb=lambda x: x['dist'] / 1000)

# Format subnet
net['score'] = 1
net['link'] = net['source'] + ' -> ' + net['target']
mat = net.groupby(['link', 'name'])['score'].sum().unstack(fill_value=0)
r_order = net.groupby(['link'])['score'].sum().sort_values(ascending=False).index
c_order = mthds + baselines
for c in c_order:
    if c not in mat:
        mat[c] = 0
mat = mat.loc[r_order, c_order]

# Plot
figs = []

plot_heatmap(sims, col='tf_oc', title='TFs', prefix='name', mthds=mthds, baselines=baselines, figs=figs, method_names=method_names)
plot_heatmap(sims, col='cre_oc', title='CREs', prefix='name', mthds=mthds, baselines=baselines, figs=figs, method_names=method_names)
plot_heatmap(sims, col='target_oc', title='Genes', prefix='name', mthds=mthds, baselines=baselines, figs=figs, method_names=method_names)
plot_heatmap(sims, col='edge_oc', title='Edges', prefix='name', mthds=mthds, baselines=baselines, figs=figs, method_names=method_names)
plot_heatmap(tss, col='ocoef', title='TSS', prefix='tss', mthds=mthds, baselines=['promoters'], figs=figs, method_names=method_names)

fig, ax = plt.subplots(1, 1, figsize=(2, 3), dpi=150)

# Map mth column to pretty names for boxplot
dst_plot = dst.copy()
dst_plot['mth'] = prettify(dst_plot['mth'].tolist(), method_names)
pretty_palette = {method_names.get(k, k): v for k, v in palette.items()}
pretty_order = prettify(list(palette.keys()), method_names)

sns.boxplot(
    data=dst_plot,
    x='kb',
    y='mth',
    hue='mth',
    palette=pretty_palette,
    ax=ax,
    order=pretty_order,
    fliersize=0,
    fill=False,
)
ax.set_xticks([0, 250, 500, 750])
ax.set_xlim(-50, 800)
ax.axvline(x=250, ls='--', color='gray')
ax.set_xlabel('Distance to TSS (kb)')
ax.set_ylabel('')
figs.append(fig)

barstats(stats, col='n_tfs', title='Number TFs', figs=figs, method_names=method_names)
barstats(stats, col='n_cres', title='Number CREs', figs=figs, method_names=method_names)
barstats(stats, col='n_targets', title='Number Genes', figs=figs, method_names=method_names)
barstats(stats, col='n_edges', title='Number Edges', figs=figs, method_names=method_names)
barstats(stats, col='odegree', title='Regulon size', figs=figs, method_names=method_names)
barstats(stats, col='betweenc', title='B. centrality', figs=figs, method_names=method_names)
barstats(stats, col='eigv', title='E. centrality', figs=figs, method_names=method_names)

h = ma.Heatmap(mat, cmap='Purples', width=mat.shape[1] * 0.15, height=mat.shape[0] * 0.15, label="Overlap\nCoefficient", vmin=0, vmax=1)
h.add_bottom(mp.Labels(prettify(mat.columns, method_names)))
h.add_left(mp.Labels(mat.index))
h.render()
plt.close()
figs.append(h.figure)

# Write
savefigs(figs, args.path_out)
