import pandas as pd
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

mdl_path = args['path_input']
path_out = args['path_out']

# Find paths
path = os.path.dirname(mdl_path)
names = os.path.basename(mdl_path)
lst = names.replace('.mdl.csv', '').split('.')
mdl = pd.read_csv(mdl_path, dtype={'source': 'str', 'target': 'str', 'score': 'float32', 'pval': 'float32'})
# Skip if empty
if len(mdl) == 0:
    grn = pd.DataFrame(columns=['source', 'cre', 'target', 'score', 'pval'])
    grn.to_csv(path_out, index=False)
    os._exit(0)
tfs = set(mdl['source'])
gns = set(mdl['target'])

# Skip baselines
baselines = ['collectri', 'dorothea', 'random', 'scenic']
if all(x == lst[0] for x in lst):
    if (lst[0] in baselines) or lst[0].startswith('o_'):
        mdl.to_csv(path_out, index=False)
        os._exit(0)

# Read
pre_name, p2g_name, tfb_name, mdl_name = lst
p2g_path = os.path.join(path, '{pre}.{p2g}.p2g.csv'.format(pre=pre_name, p2g=p2g_name))
tfb_path = os.path.join(path, '{pre}.{p2g}.{tfb}.tfb.csv'.format(pre=pre_name, p2g=p2g_name, tfb=tfb_name))

# Open files
tfb = pd.read_csv(tfb_path, usecols=['tf', 'cre'], dtype={'tf': 'str', 'cre': 'str'})
tfb = tfb[tfb['tf'].isin(tfs)]
p2g = pd.read_csv(p2g_path, usecols=['cre', 'gene'], dtype={'cre': 'str', 'gene': 'str'})
p2g = p2g[p2g['gene'].isin(gns)]

# Merge dfs (can contain duplicates at cre level)
grn = pd.merge(tfb, p2g, on='cre', how='inner').rename(columns={'tf': 'source', 'gene': 'target'})
grn = pd.merge(grn, mdl, on=['source', 'target'], how='inner').sort_values(['source', 'target', 'cre']).reset_index(drop=True)
grn = grn[['source', 'cre', 'target', 'score', 'pval']]

grn.to_csv(path_out, index=False)
