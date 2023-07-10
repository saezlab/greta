import pandas as pd
import numpy as np
import celloracle as co
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-l','--path_links', required=True)
parser.add_argument('-b','--path_basegrn', required=True)
parser.add_argument('-p','--thr_edge_pval', required=True)
parser.add_argument('-t','--thr_top_edges', required=True)
parser.add_argument('-g','--path_grn', required=True)
parser.add_argument('-r','--path_tri', required=True)
args = vars(parser.parse_args())

path_links = args['path_links']
path_basegrn = args['path_basegrn']
thr_edge_pval = float(args['thr_edge_pval'])
thr_top_edges = int(args['thr_top_edges'])
path_grn = args['path_grn']
path_tri = args['path_tri']

# Filter links and concatenate
links = co.load_hdf5(file_path=path_links)
links.filter_links(p=thr_edge_pval, weight="coef_abs", threshold_number=thr_top_edges)
grn = []
for celltype in links.links_dict.keys():
    tmp = links.filtered_links[celltype].dropna()
    tmp['celltype'] = celltype
    grn.append(tmp)
grn = pd.concat(grn)[['source', 'target', 'coef_mean', 'p', 'celltype']]
grn = grn.rename(columns={'coef_mean': 'weight', 'p': 'pvals'}).reset_index(drop=True)

# Extract triplets and filter by grn
tri = pd.read_csv(path_basegrn, index_col=0)
tri = tri.melt(id_vars=['peak_id', 'gene_short_name'], var_name='TF', value_name='interaction')
tri = tri[tri['interaction'] != 0]
tri = tri[['TF', 'peak_id', 'gene_short_name']].rename(columns={'TF': 'tf', 'peak_id': 'peak', 'gene_short_name': 'gene'})
tfs = grn['source'].unique().astype('U')
gns = grn['target'].unique().astype('U')
tri = tri.loc[np.isin(tri['tf'], tfs)]
tri = tri.loc[np.isin(tri['gene'], gns)].reset_index(drop=True)

# Save both grns
grn.to_csv(path_grn, index=False)
tri.to_csv(path_tri, index=False)