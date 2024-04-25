import scanpy as sc
import pandas as pd
import numpy as np
import os
import snapatac2 as snap
from snapatac2.datasets import _datasets, datasets
from pathlib import Path
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-r','--path_rna', required=True)
parser.add_argument('-g','--path_geneids', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-t','--path_tmp', required=True)
parser.add_argument('-f','--path_frags', required=True)
args = vars(parser.parse_args())

# Get args
path_rna = args['path_rna']
path_geneids = args['path_geneids']
organism = args['organism']
path_tmp = args['path_tmp']
path_frags = args['path_frags']

n_jobs=8

# Change default cache dir
if not os.path.exists(path_tmp):
    os.mkdir(path_tmp)
_datasets = datasets()
_datasets.path = Path(path_tmp)

# Read raw rna
rna = sc.read_10x_h5(path_rna)

# Filter faulty gene symbols
path_geneids = os.path.join(path_geneids, organism + '.csv')
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()
ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
msk = ensmbls != ''
rna = rna[:, msk].copy()

# Basic QC
sc.pp.filter_cells(rna, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)

# Remove duplicated genes based on num of cells
to_remove = []
for dup in rna.var.index[rna.var.index.duplicated()]:
    tmp = rna.var.loc[dup]
    max_idx = tmp.set_index('gene_ids')['n_cells'].idxmax()
    to_remove.extend(tmp['gene_ids'][tmp['gene_ids'] != max_idx].values)
rna = rna[:, ~rna.var['gene_ids'].isin(to_remove)].copy()

# General QC
rna.var["mt"] = rna.var_names.str.startswith("MT-")
sc.pp.calculate_qc_metrics(
    rna, qc_vars=["mt"], inplace=True, log1p=False
)
rna = rna[(rna.obs['n_genes_by_counts'] < 5000) & (rna.obs['total_counts'] < 20000) & (rna.obs['pct_counts_mt'] < 25), :].copy()

# Remove doublets
sc.external.pp.scrublet(rna)
rna = rna[rna.obs['doublet_score'] < 0.1].copy()

# Normalize
rna.layers["counts"] = rna.X.copy()
sc.pp.normalize_total(rna, target_sum=10000)
sc.pp.log1p(rna)

# Dim reduction
sc.pp.highly_variable_genes(rna, n_top_genes=2048)
sc.tl.pca(rna)
sc.pp.neighbors(rna)
sc.tl.leiden(rna)
sc.tl.umap(rna)

# Delete potential erythroblasts and empty genes
rna = rna[rna.obs['leiden'] != '13'].copy()
sc.pp.filter_genes(rna, min_cells=3)

# Annotate
annot = {
    '0': 'mono_CD14',
    '1': 'tnaive_CD8',
    '2': 'tnaive_CD4',
    '3': 'tcm_CD4',
    '4': 'tem_CD8',
    '5': 'mono_CD14',
    '6': 'bnaive',
    '7': 'nk',
    '8': 'mono_CD16',
    '9': 'bmemory',
    '10': 'cdc',
    '11': 'treg',
    '12': 'pdc',
}

rna.obs['celltype'] = [annot[l] for l in rna.obs['leiden']]

# Read fragment file
print('Reading fragments')
atac = snap.pp.import_data(
    path_frags,
    chrom_sizes=snap.genome.hg38,
    file=None,
    tempdir=path_tmp,
    sorted_by_barcode=False,
)

# Remove by RNA intersect
atac = atac[atac.obs.index.intersection(rna.obs_names), :].copy()

# Basic QC
snap.metrics.tsse(atac, snap.genome.hg38)
snap.pp.filter_cells(atac, min_counts=5000, min_tsse=10, max_counts=100000)

# Generation tile matrix
snap.pp.add_tile_matrix(atac, n_jobs=n_jobs)
snap.pp.select_features(atac, n_features=250000)

# Filter by doublet score
snap.pp.scrublet(
    atac,
    n_jobs=n_jobs
)
snap.pp.filter_doublets(atac)

# Call peaks
atac.obs['celltype'] = rna.obs['celltype'].copy()
snap.tl.macs3(atac, groupby='celltype', n_jobs=n_jobs, tempdir=path_tmp)
peaks = snap.tl.merge_peaks(atac.uns['macs3'], snap.genome.hg38)
atac = snap.pp.make_peak_matrix(atac, use_rep=peaks['Peaks'])
sc.pp.filter_cells(atac, min_genes=100)
sc.pp.filter_genes(atac, min_cells=3)

# Compute dim reduction
snap.tl.spectral(atac, features=None)
snap.tl.umap(atac, random_state=0)

# Remove by final intersect
inter = rna.obs_names.intersection(atac.obs_names)
rna = rna[inter, :].copy()
atac = atac[inter, :].copy()

# Final QC
sc.pp.filter_genes(atac, min_cells=3)
sc.pp.filter_cells(atac, min_genes=100)
sc.pp.filter_genes(rna, min_cells=3)
sc.pp.filter_cells(rna, min_genes=100)

# Generate shared emmbeding
atac.obs['X_joint'] = snap.tl.multi_spectral([rna, atac], features=None)[1]
snap.tl.umap(atac, use_rep='X_joint')

# Clean objects
rna.obs = rna.obs[['celltype']].copy()
del rna.var
del rna.uns
del rna.varm
del rna.obsp

# Renormalize
rna.X = rna.layers['counts'].copy()
sc.pp.normalize_total(rna, target_sum=10000)
sc.pp.log1p(rna)
del rna.uns
atac.layers['counts'].copy()
sc.pp.normalize_total(atac, target_sum=10000)
sc.pp.log1p(atac)
del atac.uns

atac.write('datasets/pbmc10k/atac.h5ad')
rna.write('datasets/pbmc10k/rna.h5ad')
