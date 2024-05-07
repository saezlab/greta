import matplotlib.pyplot as plt  # Import else compiling error
import numpy as np
import pandas as pd
import muon as mu
import celloracle as co
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-g','--path_p2g', required=True)
parser.add_argument('-t','--path_tfb', required=True)
parser.add_argument('-a','--alpha', required=True)
parser.add_argument('-p','--pthr', required=True)
parser.add_argument('-n','--top_n', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
path_p2g = args['path_p2g']
path_tfb = args['path_tfb']
alpha = float(args['alpha'])
pthr = float(args['pthr'])
top_n = int(args['top_n'])
path_out = args['path_out']

# Process base GRN
p2g = pd.read_csv(path_p2g)
tfb = pd.read_csv(path_tfb)
if (p2g.shape[0] == 0) or (tfb.shape[0] == 0):
    grn = pd.DataFrame(columns=['source', 'target', 'score', 'pval'])
    grn.to_csv(path_out, index=False)
    exit()
tfb['score'] = 1
p2g = p2g[['cre', 'gene']]
base_grn = pd.merge(
    p2g,
    tfb
    .pivot(index='cre', columns='tf')
    .fillna(0)
    .droplevel(0, axis=1)
    .reset_index()
)
base_grn = base_grn.rename(columns={'cre': 'peak_id', 'gene': 'gene_short_name'})
base_grn['peak_id'] = base_grn['peak_id'].str.replace('-', '_')

# Init oracle object
oracle = co.Oracle()
oracle.adata = mu.read(path_mdata)['rna'].copy()
oracle.adata.obsm['X_umap'] = np.zeros((oracle.adata.shape[0], 2))
oracle.adata.layers['imputed_count'] = oracle.adata.X
oracle.adata.obs['cluster'] = 'cluster'
oracle.cluster_column_name = 'cluster'
oracle.embedding_name = 'X_umap'
oracle.pcs = np.zeros((oracle.adata.shape[0], 2))
oracle.knn = True
oracle.k_knn_imputation = True
oracle.import_TF_data(TF_info_matrix=base_grn)

# Model TF ~ G
print('Modeling GRN...')
links = oracle.get_links(
    cluster_name_for_GRN_unit="cluster",
    alpha=alpha,
    n_jobs=32,
)
print('Modeling Done!')
print('Filtering links...')
links.filter_links(
    p=pthr,
    weight="coef_abs",
    threshold_number=top_n
)
print('Filtering done!')

# Extract grn
grn = links.filtered_links['cluster'].dropna()[['source', 'target', 'coef_mean', 'p']]
grn = grn.rename(columns={'coef_mean': 'score', 'p': 'pval'})
grn = grn.sort_values(['source', 'target', 'pval'])

# Write
grn.to_csv(path_out, index=False)

print('Done')
os._exit(0)  # Add this else it gets stuck