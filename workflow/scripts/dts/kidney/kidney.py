import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_gex', required=True)
parser.add_argument('-b','--path_peaks', required=True)
parser.add_argument('-c','--path_annot', required=True)
parser.add_argument('-e','--path_geneids', required=True)
parser.add_argument('-f','--path_output', required=True)
args = vars(parser.parse_args())

path_gex = args['path_gex']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_geneids = args['path_geneids']
path_output = args['path_output']

# Read
rna = sc.read_h5ad(path_gex)
atac = sc.read_h5ad(path_peaks)
obs = pd.read_csv(path_annot, index_col=0)
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()
# Filter faulty gene symbols
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()
rna.var['gene_ids'] = ensmbls[msk]
# Basic QC
sc.pp.filter_cells(rna, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)
del rna.obs['n_genes']
# Remove duplicated genes based on num of cells
to_remove = []
for dup in rna.var.index[rna.var.index.duplicated()]:
    tmp = rna.var.loc[dup]
    max_idx = tmp.set_index('gene_ids')['n_cells'].idxmax()
    to_remove.extend(tmp['gene_ids'][tmp['gene_ids'] != max_idx].values)
rna = rna[:, ~rna.var['gene_ids'].isin(to_remove)].copy()
del rna.obs
del rna.var
# Match
atac = atac[rna.obs_names].copy()
# Dtype
rna.X = rna.X.astype(np.float32)
# Create mdata
mdata = mu.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)
# Write
mdata.write(path_output)
