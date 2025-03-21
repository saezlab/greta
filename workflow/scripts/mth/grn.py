import pandas as pd
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i', '--path_input', required=True)
parser.add_argument('-o', '--path_out', required=True)
args = vars(parser.parse_args())

mdl_path = args['path_input']
path_out = args['path_out']

# Find paths
path = os.path.dirname(mdl_path)
names = os.path.basename(mdl_path)
lst = names.replace('.mdl.csv', '').split('.')

# Read in chunks to reduce memory usage
chunksize = 100_000  # Adjust based on available memory
dtype_dict = {'source': 'category', 'target': 'category', 'score': 'float32', 'pval': 'float32'}
mdl_chunks = pd.read_csv(mdl_path, dtype=dtype_dict, chunksize=chunksize)
mdl = pd.concat(mdl_chunks, ignore_index=True)

# Skip if empty
if mdl.empty:
    grn = pd.DataFrame(columns=['source', 'cre', 'target', 'score', 'pval'])
    grn.to_csv(path_out, index=False)
    os._exit(0)

# Limit to 100k largest absolute scores
mdl = mdl.nlargest(100_000, 'score', keep='all').reset_index(drop=True)
tfs = set(mdl['source'].unique())
gns = set(mdl['target'].unique())

# Skip baselines
baselines = {'collectri', 'dorothea', 'random', 'scenic'}
if lst[0] in baselines or lst[0].startswith('o_'):
    mdl.to_csv(path_out, index=False)
    os._exit(0)

# Read paths
pre_name, p2g_name, tfb_name, mdl_name = lst
p2g_path = os.path.join(path, f'{pre_name}.{p2g_name}.p2g.csv')
tfb_path = os.path.join(path, f'{pre_name}.{p2g_name}.{tfb_name}.tfb.csv')

# Read relevant columns with filtering
usecols_tfb = ['tf', 'cre']
usecols_p2g = ['cre', 'gene']

tfb_chunks = pd.read_csv(tfb_path, usecols=usecols_tfb, dtype={'tf': 'category', 'cre': 'category'}, chunksize=chunksize)
tfb = pd.concat((chunk[chunk['tf'].isin(tfs)] for chunk in tfb_chunks), ignore_index=True)

p2g_chunks = pd.read_csv(p2g_path, usecols=usecols_p2g, dtype={'cre': 'category', 'gene': 'category'}, chunksize=chunksize)
p2g = pd.concat((chunk[chunk['gene'].isin(gns)] for chunk in p2g_chunks), ignore_index=True)

# Merge in an optimized manner
grn = tfb.merge(p2g, on='cre', how='inner')
grn = grn.rename(columns={'tf': 'source', 'gene': 'target'})
grn = grn.merge(mdl, on=['source', 'target'], how='inner')
grn = grn.sort_values(['source', 'target', 'cre']).reset_index(drop=True)
grn = grn[['source', 'cre', 'target', 'score', 'pval']]

grn.to_csv(path_out, index=False)
