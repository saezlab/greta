import scanpy as sc
import pandas as pd
import numpy as np
import anndata as ad
import mudata as md
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_matrix_d1m', required=True)
parser.add_argument('-b','--path_barcodes_d1m', required=True)
parser.add_argument('-c','--path_matrix_d2m', required=True)
parser.add_argument('-d','--path_barcodes_d2m', required=True)
parser.add_argument('-e','--path_gsym', required=True)
parser.add_argument('-f','--path_peaks', required=True)
parser.add_argument('-g','--path_annot', required=True)
parser.add_argument('-i','--path_barmap', required=True)
parser.add_argument('-j','--path_geneids', required=True)
parser.add_argument('-k','--organism', required=True)
parser.add_argument('-l','--path_output', required=True)
args = vars(parser.parse_args())

path_matrix_d1m = args['path_matrix_d1m']
path_barcodes_d1m = args['path_barcodes_d1m']
path_matrix_d2m = args['path_matrix_d2m']
path_barcodes_d2m = args['path_barcodes_d2m']
path_gsym = args['path_gsym']
path_peaks = args['path_peaks']
path_annot = args['path_annot']
path_barmap = args['path_barmap']
path_geneids = args['path_geneids']
organism = args['organism']
path_output = args['path_output']

# Read annots and barmap
obs = pd.read_csv(path_annot, index_col=0)
bar_map = pd.read_csv(path_barmap, sep='\t').set_index('RNA_bc')['ATAC_bc'].to_dict()

# Read gene ids
path_geneids = os.path.join(path_geneids, organism + '.csv')
geneids = pd.read_csv(path_geneids).set_index('symbol')['id'].to_dict()

def read_sample(path_matrix, path_barcodes, path_gsym, bar_map, obs, geneids):
    rna = sc.read_mtx(path_matrix).T
    rna.obs_names = pd.read_csv(path_barcodes, sep='\t', header=None)[0].values
    var = pd.read_csv(path_gsym, sep='\t', header=None).set_index(1).drop(columns=[2]).rename(columns={0: 'gene_ids'})
    var.index.name = None
    rna.var = var

    # Extract sample id
    sample_id = os.path.basename(path_matrix).split('.')[0]

    # Update barcodes
    rna.obs_names = np.array([bar_map[b.split('-1')[0]] for b in rna.obs_names])
    msk = rna.obs_names.isin(obs[obs['batch'] == sample_id].index)
    rna = rna[msk, :].copy()
    # Add celltype annotation
    rna.obs['celltype'] = obs[obs['batch'] == sample_id]['celltype']
    rna.obs['batch'] = sample_id
    rna.obs_names = [sample_id + '_' + o for o in rna.obs_names]
    
    # Filter faulty gene symbols
    ensmbls = np.array([geneids[g] if g in geneids else '' for g in rna.var_names])
    msk = ensmbls != ''
    rna = rna[:, msk].copy()
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
    return rna


# Read samples
rna_d1m = read_sample(path_matrix_d1m, path_barcodes_d1m, path_gsym, bar_map, obs, geneids)
rna_d2m = read_sample(path_matrix_d2m, path_barcodes_d2m, path_gsym, bar_map, obs, geneids)

# Merge
rna = [rna_d1m, rna_d2m]
rna = ad.concat(rna, join='outer')

# Clean
obs = rna.obs.copy()
del rna.uns
del rna.obs

# Read atac data
atac = ad.read_h5ad(path_peaks)
atac = atac[rna.obs_names].copy()

# Create mdata
mdata = md.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)

# Write
mdata.write(path_output)
