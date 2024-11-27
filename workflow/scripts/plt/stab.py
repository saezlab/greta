import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
from utils import read_config, savefigs


def base_stability(df, col, title, mthds, baselines, palette, figs):
    fig, axes = plt.subplots(2, 1, figsize=(1.7, 3), dpi=150, gridspec_kw={'height_ratios': [len(mthds), len(baselines)]})
    ax = axes[0]
    sns.boxplot(data=df[df['cat'].isin(['full']) & df['mth'].isin(mthds)], y='mth', x=col, hue='mth', ax=ax, fill=False, palette=palette)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Methods')
    ax.set_title(title)
    ax.set_xlim(-0.05, 1.05)
    
    ax = axes[1]
    sns.boxplot(data=df[df['cat'].isin(['full']) & df['mth'].isin(baselines)], y='mth', x=col, hue='mth', ax=ax, fill=False, palette=palette)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('Overlap coefficient')
    ax.set_ylabel('Baselines')
    ax.set_xlim(-0.05, 1.05)
    
    fig.subplots_adjust(wspace=0.05, hspace=0)
    figs.append(fig)


def sampled_stability(df, col, ylabel, palette, figs, plot_diag=False):
    fig, axes = plt.subplots(1, 2, figsize=(4, 2), dpi=150, sharex=False)
    axes = axes.ravel()
    ax = axes[0]
    sns.pointplot(data=df[df['cat'].isin(['fixed_nfeats', 'full'])], x='n', y=col, hue='mth', ax=ax, errorbar=None, palette=palette)
    ax.legend().set_visible(False)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylabel(ylabel)
    ax.set_xlabel('Number of features')
    if plot_diag:
        ax.axline([0, 0], [4, 1], linestyle='--', c='gray')
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([int(x.get_text()) * 5 for x in ax.get_xticklabels()])
    
    ax = axes[1]
    sns.pointplot(data=df[df['cat'].isin(['fixed_ncells', 'full'])], x='n', y=col, hue='mth', ax=ax, errorbar=None, palette=palette)
    ax.legend().set_visible(False)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylabel('')
    ax.set_xlabel('Number of cells')
    if plot_diag:
        ax.axline([0, 0], [4, 1], linestyle='--', c='gray')
    ax.set_yticks([])
    ax.set_yticklabels([])
    
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False)
    fig.subplots_adjust(wspace=0., hspace=0)
    figs.append(fig)


# Read
df = pd.read_csv(sys.argv[1])

# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

# Plot
figs = []
base_stability(df, col='s_ocoeff', title='TFs', mthds=mthds, baselines=baselines, palette=palette, figs=figs)
base_stability(df, col='e_ocoeff', title='Edges', mthds=mthds, baselines=baselines, palette=palette, figs=figs)
base_stability(df, col='t_ocoeff', title='Genes', mthds=mthds, baselines=baselines, palette=palette, figs=figs)

sampled_stability(df, col='s_ocoeff', ylabel='TF Overlap coefficient', palette=palette, plot_diag=True, figs=figs)
sampled_stability(df, col='e_ocoeff', ylabel='Edge Overlap coefficient', palette=palette, plot_diag=True, figs=figs)
sampled_stability(df, col='t_ocoeff', ylabel='Gene Overlap coefficient', palette=palette, plot_diag=True, figs=figs)

sampled_stability(df, col='n_sources', ylabel='Number of TFs', palette=palette, figs=figs)
sampled_stability(df, col='n_edges', ylabel='Number of edges', palette=palette, figs=figs)
sampled_stability(df, col='n_targets', ylabel='Number of Genes', palette=palette, figs=figs)
sampled_stability(df, col='r_size', ylabel='Regulon size', palette=palette, figs=figs)

sampled_stability(df, col='h', ylabel='Time (hours)', palette=palette, figs=figs)
sampled_stability(df, col='gb', ylabel='Memory (GBs)', palette=palette, figs=figs)

# Write
savefigs(figs, sys.argv[2])
