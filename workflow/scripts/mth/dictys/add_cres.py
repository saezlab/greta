import pandas as pd
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p', '--path_p2g', required=True)
parser.add_argument('-t', '--path_tfb', required=True)
parser.add_argument('-m', '--path_mdl', required=True)
parser.add_argument('-o', '--path_out', required=True)
args = vars(parser.parse_args())

path_p2g = args['path_p2g']
path_tfb = args['path_tfb']
path_mdl = args['path_mdl']
path_out = args['path_out']

# Read in chunks to reduce memory usage
chunksize = 100_000  # Adjust based on available memory
dtype_dict = {'source': 'category', 'target': 'category', 'score': 'float32', 'pval': 'float32'}
mdl_chunks = pd.read_csv(path_mdl, dtype=dtype_dict, chunksize=chunksize)
mdl = pd.concat(mdl_chunks, ignore_index=True)

# Skip if empty
if mdl.empty:
    grn = pd.DataFrame(columns=['source', 'cre', 'target', 'score', 'pval'])
    grn.to_csv(path_out, index=False)
    os._exit(0)

# Limit to 100k largest absolute scores
mdl['abs_score'] = mdl['score'].abs()
mdl = mdl.nlargest(100_000, 'abs_score', keep='all').reset_index(drop=True)
mdl = mdl.drop(columns=['abs_score'])
tfs = set(mdl['source'].unique())
gns = set(mdl['target'].unique())

# Read relevant columns with filtering
usecols_tfb = ['tf', 'cre']
usecols_p2g = ['cre', 'gene']

tfb_chunks = pd.read_csv(path_tfb, usecols=usecols_tfb, dtype={'tf': 'category', 'cre': 'category'}, chunksize=chunksize)
tfb = pd.concat((chunk[chunk['tf'].isin(tfs)] for chunk in tfb_chunks), ignore_index=True)

p2g_chunks = pd.read_csv(path_p2g, usecols=usecols_p2g, dtype={'cre': 'category', 'gene': 'category'}, chunksize=chunksize)
p2g = pd.concat((chunk[chunk['gene'].isin(gns)] for chunk in p2g_chunks), ignore_index=True)

# Merge in an optimized manner
grn = tfb.merge(p2g, on='cre', how='inner')
grn = grn.rename(columns={'tf': 'source', 'gene': 'target'})
grn = grn.merge(mdl, on=['source', 'target'], how='inner')
grn = grn.sort_values(['source', 'target', 'cre']).reset_index(drop=True)
grn = grn[['source', 'cre', 'target', 'score', 'pval']]

grn.to_csv(path_out, index=False)
