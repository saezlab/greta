import h5py
import anndata as ad
import glob
import os
import scipy.sparse as sps
import pandas as pd
import numpy as np
import scanpy as sc
import sys

path_ann = sys.argv[1]
path_gid = sys.argv[2]

# Read
ann = pd.read_csv(path_ann, index_col=0)
geneids = pd.read_csv(path_gid).set_index('symbol')['id'].to_dict()
path_dir = os.path.dirname(path_ann)
fnames = glob.glob(os.path.join(path_dir, '*_filtered_feature_bc_matrix.h5'))
adata = []
for fname in fnames:
    f = h5py.File(fname, 'r')
    sname = os.path.basename(fname).split('_')[1]
    ftype = f['matrix']['features']['feature_type'][:].astype('U')
    msk_col = ftype == 'Gene Expression'
    barcodes = f['matrix']['barcodes'][:].astype('U')
    barcodes = [sname + '_' + b.split('-')[0] for b in barcodes]
    genes = f['matrix']['features']['name'][:].astype('U')[msk_col]
    X = sps.csr_matrix((f['matrix']['data'][:], f['matrix']['indices'][:], f['matrix']['indptr'][:]))
    X = X[:, msk_col]
    obs = pd.DataFrame(index=barcodes)
    var = pd.DataFrame(index=genes)
    rna = ad.AnnData(X=X, obs=obs, var=var)
    # Filter by obs
    msk = rna.obs_names.isin(ann.index)
    rna = rna[msk, :].copy()
    # Filter faulty gene symbols
    ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
    msk = ensmbls != ''
    rna = rna[:, msk].copy()
    rna.var['gene_ids'] = ensmbls[msk]
    # Basic QC
    sc.pp.filter_cells(rna, min_genes=100)
    sc.pp.filter_genes(rna, min_cells=3)
    # Remove duplicated genes based on num of cells
    to_remove = []
    for dup in rna.var.index[rna.var.index.duplicated()]:
        to_remove.append(dup)
    rna = rna[:, ~rna.var_names.isin(to_remove)].copy()
    del rna.obs
    del rna.var
    adata.append(rna)
adata = ad.concat(adata, join='outer')
adata = adata[ann.index, :].copy()
adata.write(os.path.join(path_dir, 'rna.h5ad'))
