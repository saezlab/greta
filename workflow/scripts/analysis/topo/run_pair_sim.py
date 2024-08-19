import concurrent.futures
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from functools import partial
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import ocoeff, parallel_ocoeff
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--paths', required=True, nargs='+')
parser.add_argument('-t','--stat_path', required=True)
parser.add_argument('-s','--sim_path', required=True)
args = vars(parser.parse_args())

paths = args['paths']
stat_path = args['stat_path']
sim_path = args['sim_path']


def get_grn_stats(df):
    n_s = df['source'].unique().size
    n_e = df.shape[0]
    n_t = df['target'].unique().size
    n_r = df.groupby(['source']).count()['target'].mean()
    return n_s, n_e, n_t, n_r


print('Reading and computing grns stats...')
names = []
dfs = []
n_srcs = []
n_edgs = []
n_trgs = []
n_regs = []
for path in tqdm(paths):
    names.append(os.path.basename(path).replace('.grn.csv', '').replace('.csv', ''))
    df = pd.read_csv(path).drop_duplicates(['source', 'target'], keep='first')
    dfs.append(df)
    n_s, n_e, n_t, n_r = get_grn_stats(df)
    n_srcs.append(n_s)
    n_edgs.append(n_e)
    n_trgs.append(n_t)
    n_regs.append(n_r)

# Store as df
stats = pd.DataFrame()
stats['name'] = names
stats['n_tfs'] = n_srcs
stats['n_edges'] = n_edgs
stats['n_targets'] = n_trgs
stats['mean_reg_size'] = n_regs

print('Computing pairwise overlap coefficients...')
index_pairs = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names))]
names_a = []
names_b = []
tf_coefs = []
edge_coefs = []
target_coefs = []
parallel_ocoeff_with_dfs = partial(parallel_ocoeff, dfs=dfs)
with concurrent.futures.ProcessPoolExecutor() as executor:
    for res in tqdm(executor.map(parallel_ocoeff_with_dfs, index_pairs), total=len(index_pairs)):
        i, j, tf_coef, edge_coef, target_coef = res
        names_a.append(names[i])
        names_b.append(names[j])
        tf_coefs.append(tf_coef)
        edge_coefs.append(edge_coef)
        target_coefs.append(target_coef)

# Store as df
sims = pd.DataFrame()
sims['name_a'] = names_a
sims['name_b'] = names_b
sims['tf_oc'] = tf_coefs
sims['edge_oc'] = edge_coefs
sims['target_oc'] = target_coefs

# Write
stats.to_csv(stat_path, index=False)
sims.to_csv(sim_path, index=False)
