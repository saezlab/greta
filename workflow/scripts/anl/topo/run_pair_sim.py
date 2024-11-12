import concurrent.futures
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from functools import partial
import sys
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import (
    ocoeff,
    parallel_ocoeff,
    parallel_ocoeff_chunk,
    get_grn_name,
    get_grn_stats
)
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


print('Reading and computing grns stats...')
names = []
dfs = []
stats = []
for path in tqdm(paths):
    name = get_grn_name(path)
    names.append(name)
    df = pd.read_csv(path).drop_duplicates(['source', 'target'], keep='first')
    dfs.append(df)
    stat = get_grn_stats(df)
    stats.append(['name'] + list(stat))

# Store as df
cols = ['name', 'n_tfs', 'n_edges', 'n_targets', 'mean_reg_size', 'tf_out_degree', 'tf_betweenc', 'tf_eigv']
stats = pd.DataFrame(stats, columns=cols)

print('Computing pairwise overlap coefficients...')
chunk_size = 100  # Adjust based on profiling
index_pairs = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names))]
index_pairs_chunks = [index_pairs[i:i + chunk_size] for i in range(0, len(index_pairs), chunk_size)]
names_a = []
names_b = []
tf_coefs = []
edge_coefs = []
target_coefs = []
total_pairs = len(index_pairs)
processed_pairs = 0
with tqdm(total=total_pairs, desc="Processing", unit="pair") as pbar:
    with concurrent.futures.ProcessPoolExecutor(max_workers=64) as executor:
        for chunk_result in tqdm(executor.map(partial(parallel_ocoeff_chunk, dfs=dfs), index_pairs_chunks)):
            for res in chunk_result:
                i, j, tf_coef, edge_coef, target_coef = res
                names_a.append(names[i])
                names_b.append(names[j])
                tf_coefs.append(tf_coef)
                edge_coefs.append(edge_coef)
                target_coefs.append(target_coef)
            processed_pairs += len(chunk_result)
            pbar.update(len(chunk_result))

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
