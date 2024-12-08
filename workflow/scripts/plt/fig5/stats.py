import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import matplotlib
import numpy as np
from matplotlib import ticker as mticker
import sys


# Read
df = pd.read_csv(sys.argv[1])

# Process
meta = df.drop_duplicates(['metric', 'name']).groupby('metric', as_index=False).size()
meta['metric'] = pd.Categorical(meta['metric'], categories=df['metric'].unique(), ordered=True)
meta = meta.sort_values('metric').reset_index(drop=True)
n = meta.shape[0]
heights = (meta['size'] / meta['size'].max()).values

# Plot
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

# Write
fig.savefig(sys.argv[2], bbox_inches='tight')
