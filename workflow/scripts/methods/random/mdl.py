import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--p2g_path', required=True)
parser.add_argument('-t','--tfb_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

p2g_path = args['p2g_path']
tfb_path = args['tfb_path']
out_path = args['out_path']

# Read
p2g = pd.read_csv(p2g_path)
tfb = pd.read_csv(tfb_path)
if (p2g.shape[0] == 0) or (tfb.shape[0] == 0):
    grn = pd.DataFrame(columns=['source', 'target'])
    grn.to_csv(path_out, index=False)
    exit()

# Randomize
name = os.path.basename(tfb_path).replace('.tfb.csv', '')
if name != 'random.random.random':
    rng = np.random.default_rng(seed=42)
    for col in p2g.columns:
        p2g[col] = rng.choice(p2g[col], p2g.shape[0], replace=False)
    for col in tfb.columns:
        tfb[col] = rng.choice(tfb[col], tfb.shape[0], replace=False)

# Join
df = pd.merge(tfb[['tf', 'cre']], p2g[['cre', 'gene']], how='inner', on='cre')[['tf', 'gene']]
df = df.sort_values(['tf', 'gene']).rename(columns={'tf': 'source', 'gene': 'target'}).drop_duplicates(['source', 'target'])

# Write
df.to_csv(out_path, index=False)
