import anndata as ad
import pandas as pd
import scanpy as sc
import numpy as np
import argparse
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_rann', required=True)
parser.add_argument('-b','--path_gid', required=True)
parser.add_argument('-c','--samples', required=True, nargs='+')
parser.add_argument('-d','--path_rna', required=True)
parser.add_argument('-e','--path_ann', required=True)
args = vars(parser.parse_args())
path_rann = args['path_rann']
path_gid = args['path_gid']
samples = args['samples']
path_rna = args['path_rna']
path_ann = args['path_ann']
# Read
obs = pd.read_csv(path_rann, index_col=0)
geneids = pd.read_csv(path_gid).set_index('symbol')['id'].to_dict()
# Filter obs
obs = obs[obs['batch'].isin(samples)]
ncells = obs.groupby('celltype').size()
ncells = ncells[ncells >= 100]
obs = obs[obs['celltype'].isin(ncells.index)]
adata = []
for sample in samples:
    print(sample)
    path_sample = f'dts/hg38/lung/{sample}_matrix.h5'
    rna = sc.read_10x_h5(path_sample)
    rna.obs.index = [f'{sample}_' + i.replace('-1', '') for i in rna.obs.index]
    # Basic filter
    sc.pp.filter_cells(rna, min_genes=100)
    sc.pp.filter_genes(rna, min_cells=3)
    # Filter by annot
    inter = rna.obs_names.intersection(obs.index)
    rna = rna[inter, :].copy()
    # Filter faulty gene symbols
    ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
    msk = ensmbls != ''
    rna = rna[:, msk].copy()
    # Remove duplicated genes based on num of cells
    to_remove = []
    for dup in rna.var.index[rna.var.index.duplicated()]:
        tmp = rna.var.loc[dup]
        max_idx = tmp.set_index('gene_ids')['n_cells'].idxmax()
        to_remove.extend(tmp['gene_ids'][tmp['gene_ids'] != max_idx].values)
    rna = rna[:, ~rna.var['gene_ids'].isin(to_remove)].copy()
    adata.append(rna)
adata = ad.concat(adata, join='outer')
del adata.obs
adata = adata[obs.index, :].copy()
# Write
adata.write(path_rna)
obs.to_csv(path_ann)
