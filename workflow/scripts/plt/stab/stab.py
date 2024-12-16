import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def base_stability(df, col, title, mthds, baselines, palette, figs):
    fig, axes = plt.subplots(2, 1, figsize=(1.7, 3), dpi=150, gridspec_kw={'height_ratios': [len(mthds), len(baselines)]})
    ax = axes[0]
    data = df[df['cat'].isin(['full']) & df['mth'].isin(mthds)].copy()
    data['seeds'] = [''.join(sorted([str(int(sa)), str(int(sb))])) for sa, sb in zip(data['seed'], data['other_seed'])]
    data = data.drop_duplicates(['mth', 'seeds'])
    sns.boxplot(data=data, y='mth', x=col, hue='mth', ax=ax, fill=False, palette=palette, fliersize=0)
    sns.stripplot(data=data, y='mth', x=col, hue='mth', ax=ax, palette=palette)
    ax.tick_params(axis='x', rotation=90)
    ax.set_xlabel('')
    ax.set_ylabel('Methods')
    ax.set_title(title)
    ax.set_xlim(-0.05, 1.05)
    
    ax = axes[1]
    data = df[df['cat'].isin(['full']) & df['mth'].isin(baselines)].copy()
    data['seeds'] = [''.join(sorted([str(int(sa)), str(int(sb))])) for sa, sb in zip(data['seed'], data['other_seed'])]
    data = data.drop_duplicates(['mth', 'seeds'])
    sns.boxplot(data=data, y='mth', x=col, hue='mth', ax=ax, fill=False, palette=palette, fliersize=0)
    sns.stripplot(data=data, y='mth', x=col, hue='mth', ax=ax, palette=palette)
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
    ax.set_xlabel('Number of cells')
    if plot_diag:
        ax.axline([0, 0], [4, 1], linestyle='--', c='gray')
    for line in ax.lines:
        line.set_marker('')
    
    ax = axes[1]
    sns.pointplot(data=df[df['cat'].isin(['fixed_ncells', 'full'])], x='n', y=col, hue='mth', ax=ax, errorbar=None, palette=palette)
    ax.legend().set_visible(False)
    ax.tick_params(axis='x', rotation=90)
    ax.set_ylabel('')
    ax.set_xlabel('Number of features')
    if plot_diag:
        ax.axline([0, 0], [4, 1], linestyle='--', c='gray')
    ax.set_yticks([])
    ax.set_yticklabels([])
    ax.set_xticks(ax.get_xticks())
    ax.set_xticklabels([int(x.get_text()) * 5 for x in ax.get_xticklabels()])
    for line in ax.lines:
        line.set_marker('')
    
    handles, labels = ax.get_legend_handles_labels()
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False)
    fig.subplots_adjust(wspace=0., hspace=0)
    figs.append(fig)


def auc(df, typ, title, palette, figs):
    data = df[df['type'] == typ].pivot(index='mth', columns='cat', values='auc').reset_index()
    fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
    sns.scatterplot(data, y='fixed_nfeats', x='fixed_ncells', hue='mth', palette=palette)
    ax.set_xlim(-0.05, 1.05)
    ax.set_ylim(-0.05, 1.05)
    ax.axvline(x=0.5, ls='--', color='gray')
    ax.axhline(y=0.5, ls='--', color='gray')
    ax.set_xlabel('Stability number features')
    ax.set_ylabel('Stability number cells')
    ax.set_xticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_yticks([0, 0.25, 0.5, 0.75, 1])
    ax.set_title(title)
    handles, labels = ax.get_legend_handles_labels()
    ax.legend().set_visible(False)
    fig.legend(handles, labels, loc='center left', bbox_to_anchor=(0.9, 0.5), frameon=False)
    figs.append(fig)


# Read
df = pd.read_csv(sys.argv[1])
ac = pd.read_csv(sys.argv[2])

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

auc(ac, typ='s_ocoeff', title='TFs', palette=palette, figs=figs)
auc(ac, typ='e_ocoeff', title='Edges', palette=palette, figs=figs)
auc(ac, typ='t_ocoeff', title='Genes', palette=palette, figs=figs)

# Write
savefigs(figs, sys.argv[3])
