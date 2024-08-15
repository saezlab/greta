import pandas as pd
import numpy as np
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--mthds', required=True, nargs='+')
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

mthds = args['mthds']
out_path = args['out_path']

# Process
dataset = os.path.basename(out_path).replace('.csv', '')

def over_coef(df_a, df_b, source='source', target='target'):
    pairs_a = (df_a[source] + '|' + df_a[target]).values.astype('U')
    pairs_b = (df_b[source] + '|' + df_b[target]).values.astype('U')
    inter = np.intersect1d(pairs_a, pairs_b)
    return inter.size / np.min([pairs_a.size, pairs_b.size])


def get_grn_stats(df):
    n_s = df['source'].unique().size
    n_e = df.shape[0]
    n_t = df['target'].unique().size
    n_r = df.groupby(['source']).count()['target'].mean()
    return n_s, n_e, n_t, n_r


def read_mth_cat(dataset, mth, fix_ncells):
    ns = [10, 11, 12, 13, 14]
    seeds = [0, 1, 2]
    dfs = []
    ref = pd.read_csv('datasets/{dataset}/cases/{n_cell}_{n_gene}_{seed}/runs/{mth}.src.csv'.
                      format(dataset=dataset, n_cell=2 ** ns[-1], n_gene=2 ** ns[-1], seed=seeds[0], mth=mth))
    for n in ns:
        n = 2 ** n
        fixed_n = 2 ** 14
        for seed in seeds:
            if fix_ncells:
                n_cell = fixed_n
                n_gene = n
                cat = 'fixed_ncells'
            else:
                n_cell = n
                n_gene = fixed_n
                cat = 'fixed_ngenes'
            
            df = pd.read_csv('benchmarks/{dataset}.{n_cell}_{n_gene}_{seed}.{mth}.src.txt'.
                             format(dataset=dataset, n_cell=n_cell, n_gene=n_gene, seed=seed, mth=mth), sep='\t')
            df = df[['s', 'max_rss']]
            df['s'] = df['s'] / 60 / 60
            df['mth'] = mth
            df['n'] = n
            df['cat'] = cat
            df = df.rename(columns={'s': 'h', 'max_rss': 'mbs'})
            df = df[['mth', 'cat', 'n', 'h', 'mbs']]
            
            net = pd.read_csv('datasets/{dataset}/cases/{n_cell}_{n_gene}_{seed}/runs/{mth}.src.csv'.
                              format(dataset=dataset, n_cell=n_cell, n_gene=n_gene, seed=seed, mth=mth))
            df['ocoeff'] = over_coef(ref, net)
            df[['n_sources', 'n_edges', 'n_targets', 'r_size']] = get_grn_stats(net)
            dfs.append(df)
    dfs = pd.concat(dfs)
    return dfs


df = []
for i, mth in enumerate(mthds):
    n_cells = read_mth_cat(dataset=dataset, mth=mth, fix_ncells=True)
    n_genes = read_mth_cat(dataset=dataset, mth=mth, fix_ncells=False)
    df.append(pd.concat([n_cells, n_genes]))
df = pd.concat(df)

# Write
df.to_csv(out_path, index=False)
