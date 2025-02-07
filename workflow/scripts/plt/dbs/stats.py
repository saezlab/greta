import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import marsilea as ma
import marsilea.plotter as mp
import matplotlib
import numpy as np
from matplotlib import ticker as mticker
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import savefigs


def heatmap(df, typ, title, order, figs):
    mat = df[df['type'] == typ].pivot(index='db_a', columns='db_b', values='ocoeff')
    mat = mat.combine_first(mat.T).fillna(1)
    mat = mat.loc[order, order]
    h = ma.Heatmap(mat, cmap='Purples', width=1.5, height=1.5, name=title, label="Overlap\nCoefficient", vmin=0, vmax=1)
    h.add_bottom(mp.Labels(mat.columns))
    h.add_left(mp.Labels(mat.index))
    h.add_top(mp.Title(title))
    h.add_legends()
    h.render()
    plt.close()
    figs.append(h.figure)


# Read
df = pd.read_csv(sys.argv[1])
oc = pd.read_csv(sys.argv[2])

# Process
meta = df.drop_duplicates(['metric', 'name']).groupby('metric', as_index=False).size()
meta['metric'] = pd.Categorical(meta['metric'], categories=df['metric'].unique(), ordered=True)
meta = meta.sort_values('metric').reset_index(drop=True)
n = meta.shape[0]
heights = (meta['size'] / meta['size'].max()).values

# Plot
figs = []
fig, axes = plt.subplots(
    n, 1, figsize=(2.5, n * 1), dpi=150,
    gridspec_kw={'height_ratios': heights}, sharex=True
)


def plot_bars(data, x, y, hue, ax):
    theme_dict = {**axes_style("whitegrid"), "grid.linestyle": "-", }
    (
        so.Plot(
            data=data,
            x=x,
            y=y,
            color=hue,
        )
        .add(so.Bar(alpha=1), so.Agg(), so.Stack())
        .theme(theme_dict)
        .on(ax)
        .plot()
    )

for i, metric in enumerate(meta['metric']):
    ax = axes[i]
    sns.barplot(
        y='name',
        x='val',
        data=df[df['metric'] == metric],
        hue='type',
        width=0.95,
        ax=ax,
    )
    ax.set_xlabel('')
    ax.set_ylabel(metric)
    ax.set_xscale('log')
    ax.set_xlim(10, 1e9)
    ax.xaxis.set_major_locator(mticker.LogLocator(numticks=999))
    ax.xaxis.set_minor_locator(mticker.LogLocator(numticks=999, subs="auto"))
    ax.legend(loc='center left', bbox_to_anchor=(1, 0.5), frameon=False)
plt.subplots_adjust(wspace=0, hspace=0.5)
figs.append(fig)

order = ['knocktf', 'hpa', 'tfmdb', 'europmc', 'intact', 'chipatlas', 'remap2022', 'unibind']
heatmap(oc, typ='tf', title='TFs', order=order, figs=figs)

order = ['knocktf', 'hall', 'kegg', 'reac', 'prog', 'eqtlcatalogue']
heatmap(oc, typ='gene', title='Genes', order=order, figs=figs)

order = ['chipatlas', 'remap2022', 'unibind', 'encode', 'gwascatalogue', 'blacklist', 'phastcons', 'promoters', 'zhang21']
heatmap(oc, typ='bp', title='Base pairs', order=order, figs=figs)

# Write
savefigs(figs, sys.argv[3])