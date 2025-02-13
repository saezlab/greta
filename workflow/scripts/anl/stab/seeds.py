import scipy.stats as sts
import pandas as pd
import numpy as np
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config


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

# Find inter across seeds
seeds = [0, 1, 2]
dfs = []
for mth in mthds:
    if mth not in baselines:
        mth = 'o_' + mth
    df = []
    for i, seed_a in enumerate(seeds):
        seed_a = str(seed_a)
        path_a = f'dts/{dname}/cases/16384_16384_{seed_a}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'
        grn_a = pd.read_csv(path_a)[['source', 'target', 'score']]
        for seed_b in seeds[i + 1:]:
            path_b = f'dts/{dname}/cases/16384_16384_{seed_b}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'
            grn_b = pd.read_csv(path_b)[['source', 'target', 'score']]
            df.append(pd.merge(grn_a, grn_b, how='inner', on=['source', 'target']).assign(comp=f'{seed_a}_{seed_b}'))
    mth = mth.replace('o_', '')
    df = pd.concat(df)
    if df.shape[0] > 1:
        df.insert(0, 'mth', mth)
    else:
        df.loc[0, :] = [np.nan for c in df.columns]
        df['mth'] = mth
    dfs.append(df)
df = pd.concat(dfs)

# Cors
pairs = ['0_1', '0_2', '1_2']
cors = []
for mth in df['mth'].unique():
    tmp = df[df['mth'] == mth]
    for pair in pairs:
        comp = tmp[tmp['comp'] == pair]
        if comp.shape[0] > 1:
            r, p = sts.pearsonr(comp['score_x'], comp['score_y'])
        else:
            r, p = np.nan, 1
        cors.append([mth, r, p, pair])
cors = pd.DataFrame(cors, columns=['mth', 'stat', 'pval', 'comp'])
cors['padj'] = sts.false_discovery_control(cors['pval'])

# Write
df.to_csv(sys.argv[2], index=False)
cors.to_csv(sys.argv[3], index=False)

