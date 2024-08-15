import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--paths', required=True, nargs='+')
parser.add_argument('-m','--mthds', required=True, nargs='+')
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

paths = args['paths']
mthds = args['mthds']
out_path = args['out_path']


def get_grn_stats(df):
    n_s = df['source'].unique().size
    n_e = df.shape[0]
    n_t = df['target'].unique().size
    n_r = df.groupby(['source']).count()['target'].mean()
    return n_s, n_e, n_t, n_r


print('Reading and computing grns stats...')
names = []
n_srcs = []
n_edgs = []
n_trgs = []
n_regs = []
for path in tqdm(paths):
    names.append(os.path.basename(path).replace('.grn.csv', ''))
    df = pd.read_csv(path).drop_duplicates(['source', 'target'], keep='first')
    n_s, n_e, n_t, n_r = get_grn_stats(df)
    n_srcs.append(n_s)
    n_edgs.append(n_e)
    n_trgs.append(n_t)
    n_regs.append(n_r)

# Store as df
df = pd.DataFrame()
df['name'] = names
df['n_tfs'] = n_srcs
df['n_edges'] = n_edgs
df['n_targets'] = n_targets
df['mean_reg_size'] = n_regs

# Write
df.to_csv(out_path, index=False)