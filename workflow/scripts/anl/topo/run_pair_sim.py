import concurrent.futures
import pandas as pd
import numpy as np
import os
import glob
from tqdm import tqdm
from functools import partial
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import (
    ocoeff,
    get_grn_name,
    get_grn_stats
)
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-t','--stat_path', required=True)
parser.add_argument('-s','--sim_path', required=True)
args = vars(parser.parse_args())

stat_path = args['stat_path']
sim_path = args['sim_path']

dat, case = os.path.basename(stat_path).split('.')[:2]
paths = glob.glob(os.path.join('dts', dat, 'cases', case, 'runs', '*.grn.csv'))

print('Reading and computing grns stats...')
names = []
dfs = []
stats = []
tfs = []
edges = []
genes = []

for path in tqdm(paths):
    name = get_grn_name(path)
    names.append(name)
    df = pd.read_csv(path).drop_duplicates(['source', 'target'], keep='first')
    stat = get_grn_stats(df)
    stats.append([name] + list(stat))
    tfs.append(set(df['source']))
    edges.append(set(df['source'] + '|' + df['target']))
    genes.append(set(df['target']))
    

# Store as df
cols = ['name', 'n_tfs', 'n_edges', 'n_targets', 'odegree', 'betweenc', 'eigv']
stats = pd.DataFrame(stats, columns=cols)

print('Computing pairwise overlap coefficients...')


def set_ocoef(a, b):
    min_s = min(len(a), len(b))
    if min_s == 0:
        return np.nan
    else:
        inter = len(a & b)
        return inter / min_s


names_a = []
names_b = []
tf_coefs = []
edge_coefs = []
target_coefs = []
for i in tqdm(range(len(names))):
    name_a = names[i]
    tf_a = tfs[i]
    ed_a = edges[i]
    gn_a = genes[i]
    for j in range(i, len(names)):
        name_b = names[j]
        tf_b = tfs[j]
        ed_b = edges[j]
        gn_b = genes[j]
        names_a.append(name_a)
        names_b.append(name_b)
        tf_coefs.append(set_ocoef(tf_a, tf_b))
        edge_coefs.append(set_ocoef(ed_a, ed_b))
        target_coefs.append(set_ocoef(gn_a, gn_b))


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
