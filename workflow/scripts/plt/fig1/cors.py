import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import scipy.stats as ss
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


def robustness(dname, mth):
    def scatter(data, ax, name):
        n = np.min([1000, data.shape[0]])
        sns.scatterplot(
            data=data.sample(n=n, replace=False, random_state=42),
            x='score_x',
            y='score_y',
            ax=ax,
            color=palette[name.replace(' scores','')]
        )
        ax.set_xlabel('')
        ax.set_ylabel(name)
    seeds = [0, 1, 2]
    df = []
    for i, seed_a in enumerate(seeds):
        seed_a = str(seed_a)
        path_a = f'dts/{dname}/cases/16384_16384_{seed_a}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'
        grn_a = pd.read_csv(path_a)[['source', 'target', 'score']]
        for seed_b in seeds[i + 1:]:
            path_b = f'dts/{dname}/cases/16384_16384_{seed_b}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'
            grn_b = pd.read_csv(path_b)[['source', 'target', 'score']]
            df.append(pd.merge(grn_a, grn_b, how='inner', on=['source', 'target']).assign(comp=f'{seed_a}_{seed_b}'))
    fig, axes = plt.subplots(1, 3, figsize=(6, 2), sharey=True, sharex=True)
    axes = axes.ravel()

    mth = mth.replace('o_', '')
    cors = []
    for i in range(len(df)):
        tmp = df[i]
        s, p = ss.pearsonr(tmp['score_x'], tmp['score_y'])
        scatter(tmp, ax=axes[i], name=f'{mth} scores')
        cors.append([mth, i, s, p])
    cors = pd.DataFrame(cors, columns=['mth', 'pair', 'stat', 'pval'])
    cors['padj'] = ss.false_discovery_control(cors['pval'])
    fig.subplots_adjust(wspace=0.)
    return cors, fig


# Read config
config = read_config()
palette = config['colors']['nets']
mthds = list(config['methods'].keys())
baselines = config['baselines']

path_df = sys.argv[1]
dname = os.path.basename(path_df).split('.')[0]
df = pd.read_csv(path_df)
mthds = df[df['cat'] == 'full'].groupby('mth', as_index=False)['e_ocoeff'].mean()
mthds = mthds[mthds['e_ocoeff'] < 1.]['mth'].values

cors = []
figs = []
for mth in mthds:
    if mth != 'scenic':
        mth = 'o_' + mth
    cor, fig = robustness(dname, mth=mth)
    cors.append(cor)
    figs.append(fig)
cors = pd.concat(cors)

fig, ax = plt.subplots(1, 1, figsize=(1.5, 1), dpi=150)
sns.boxplot(data=cors, x='stat', y='mth', hue='mth', fill=None, ax=ax, palette=palette)
sns.stripplot(data=cors, x='stat', y='mth', hue='mth', ax=ax, palette=palette)
ax.set_xlabel('Pearson ρ')
ax.set_ylabel('')
ax.set_xticks([0, 0.5, 1])
ax.set_xlim(-0.05, 1.05)
figs.append(fig)

# Write
savefigs(figs, sys.argv[2])
