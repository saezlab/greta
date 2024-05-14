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
pre_name, p2g_name, tfb_name, mdl_name = names.replace('.mdl.csv', '').split('.')
p2g_path = os.path.join(path, '{pre}.{p2g}.p2g.csv'.format(pre=pre_name, p2g=p2g_name))
tfb_path = os.path.join(path, '{pre}.{p2g}.{tfb}.tfb.csv'.format(pre=pre_name, p2g=p2g_name, tfb=tfb_name))

# Open files
tfb = pd.read_csv(tfb_path)
p2g = pd.read_csv(p2g_path)
mdl = pd.read_csv(mdl_path)

# Merge dfs (can contain duplicates at cre level)
grn = pd.merge(tfb[['tf', 'cre']], p2g[['cre', 'gene']], on='cre').rename(columns={'tf': 'source', 'gene': 'target'})
grn = pd.merge(grn, mdl, on=['source', 'target'], how='inner').sort_values(['source', 'target', 'cre']).reset_index(drop=True)

# Write
grn.to_csv(path_out, index=False)
