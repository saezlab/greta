import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

path_repl_wgt = sys.argv[1]
path_repl_cor = sys.argv[2]
repl_wgt = pd.read_csv(path_repl_wgt)
repl_cor = pd.read_csv(path_repl_cor)

figs = []
for mth in repl_wgt['mth'].unique():
    tmp = repl_wgt[repl_wgt['mth'] == mth]
    if tmp.shape[0] > 1:
        fig, ax = plt.subplots(1, 1, figsize=(2, 2), dpi=150)
        max_n = np.max([tmp['score_x'].abs().max(), tmp['score_y'].abs().max()])
        max_n = max_n + (max_n * 0.05)
        sns.histplot(
            data=tmp,
            x='score_x',
            y='score_y',
            cbar=False,
            cmap='magma',
            stat='proportion',
            vmin=0.,
            vmax=1e-2,
            bins=(50, 50),
            cbar_kws=dict(label='Proportion', shrink=0.5, aspect=5, orientation='horizontal')
        )
        ax.set_xlabel('Run A edge score')
        ax.set_ylabel('Run B edge score')
        ax.set_xlim(-max_n, max_n)
        ax.set_ylim(-max_n, max_n)
        ax.set_title(mth)
        figs.append(fig)
        

fig, ax = plt.subplots(1, 1, figsize=(1.5, 1), dpi=150)
order = mthds + baselines
order = [m for m in order if m in repl_cor['mth'].unique()]
sns.boxplot(data=repl_cor, x='stat', y='mth', hue='mth', fill=None, ax=ax, palette=palette, order=order)
sns.stripplot(data=repl_cor, x='stat', y='mth', hue='mth', ax=ax, palette=palette, order=order)
ax.set_xlabel('Pearson ρ')
ax.set_ylabel('')
ax.set_xticks([0, 0.5, 1])
ax.set_xlim(-0.05, 1.05)
figs.append(fig)

# Write
savefigs(figs, sys.argv[3])
