import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-w','--path_repl_wgt', required=True)
parser.add_argument('-b','--baselines', required=True, nargs='+')
parser.add_argument('-c','--path_repl_cor', required=True)
parser.add_argument('-o','--path_out', required=True)
args = parser.parse_args()

# Read config
config = read_config()
palette_old = config['colors']['nets']
palette = {config['method_names'][k]: palette_old[k] for k in palette_old}
mthds = list(config['methods'].keys())
baselines = args.baselines
mthds = [m for m in mthds if m not in baselines]
mth_dict = config['method_names']

path_repl_wgt = args.path_repl_wgt
path_repl_cor = args.path_repl_cor
repl_wgt = pd.read_csv(path_repl_wgt)
repl_cor = pd.read_csv(path_repl_cor)

figs = []
for mth in repl_wgt['mth'].unique():
    tmp = repl_wgt[repl_wgt['mth'] == mth].copy()
    tmp['mth'] = [mth_dict[m] for m in tmp['mth']]
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
        ax.set_title(mth_dict[mth])
        figs.append(fig)
        

fig, ax = plt.subplots(1, 1, figsize=(1.5, 2), dpi=150)
order = mthds + baselines
order = [mth_dict[m] for m in order if m in repl_cor['mth'].unique()]
repl_cor['mth'] = [mth_dict[m] for m in repl_cor['mth']]

sns.boxplot(data=repl_cor, x='stat', y='mth', hue='mth', fill=None, ax=ax, palette=palette, order=order, fliersize=0)
sns.stripplot(data=repl_cor, x='stat', y='mth', hue='mth', ax=ax, palette=palette, order=order, size=0)
ax.set_xlabel('Pearson ρ')
ax.set_ylabel('')
ax.set_xticks([0, 0.5, 1])
ax.set_xlim(-0.05, 1.05)
figs.append(fig)

# Write
savefigs(figs, args.path_out)
