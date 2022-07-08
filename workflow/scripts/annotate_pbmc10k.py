import scanpy as sc
import muon as mu
from muon import atac as ac
import argparse
import os


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--input', required=True)
parser.add_argument('-p','--plot', required=True)
parser.add_argument('-o','--output', required=True)
args = vars(parser.parse_args())

data_dir = args['input']
plot = args['plot']
out_dir = args['output']

# Read and create mdata object
mdata = mu.read_10x_h5(os.path.join(data_dir, "filtered_feature_bc_matrix.h5"))
mdata.var_names_make_unique()

# Fix muon #65 bug
mdata.mod['atac'].uns['files'] = dict(mdata.mod['atac'].uns['files'])
mdata.mod['atac'].uns['atac'] = dict(mdata.mod['atac'].uns['atac'])

# Extract RNA
rna = mdata.mod['rna']

## QC
rna.var['mt'] = rna.var_names.str.startswith('MT-')
sc.pp.calculate_qc_metrics(rna, qc_vars=['mt'], percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(rna, 'n_cells_by_counts', lambda x: x >= 3)
mu.pp.filter_obs(rna, 'n_genes_by_counts', lambda x: (x >= 200) & (x < 5000))
mu.pp.filter_obs(rna, 'total_counts', lambda x: x < 15000)
mu.pp.filter_obs(rna, 'pct_counts_mt', lambda x: x < 20)

## Normalize
sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)

## HVG
sc.pp.highly_variable_genes(rna, min_mean=0.02, max_mean=4, min_disp=0.5)

## PCA
rna.raw = rna
sc.pp.scale(rna, max_value=10)
sc.tl.pca(rna, svd_solver='arpack')

## UMAP
sc.pp.neighbors(rna, n_neighbors=10, n_pcs=20)
sc.tl.leiden(rna, resolution=.5)
sc.tl.umap(rna, spread=1., min_dist=.5, random_state=11)

## Remove noisy clusters
mu.pp.filter_obs(rna, "leiden", lambda x: ~x.isin(["9", "15", "12", "16"]))

## Annotate
new_cluster_names = {
    "0": "CD4+ memory T",
    "1": "intermediate mono",
    "2": "CD4+ naïve T",
    "3": "CD8+ naïve T",
    "4": "CD14 mono",
    "5": "CD8+ activated T",
    "6": "memory B",
    "7": "NK",
    "8": "CD16 mono",
    "10": "naïve B",
    "11": "mDC",
    "13": "MAIT",
    "14": "pDC"
}
rna.obs['celltype'] = rna.obs.leiden.astype("str").values
rna.obs.celltype = rna.obs.celltype.astype("category")
rna.obs.celltype = rna.obs.celltype.cat.rename_categories(new_cluster_names)
rna.obs.celltype.cat.reorder_categories([
    'CD4+ naïve T', 'CD4+ memory T', 'MAIT',
    'CD8+ naïve T', 'CD8+ activated T', 'NK',
    'naïve B', 'memory B',
    'CD14 mono', 'intermediate mono', 'CD16 mono',
    'mDC', 'pDC'], inplace=True)

# Extract ATAC
atac = mdata.mod['atac']

## QC
sc.pp.calculate_qc_metrics(atac, percent_top=None, log1p=False, inplace=True)
mu.pp.filter_var(atac, 'n_cells_by_counts', lambda x: x >= 10)
mu.pp.filter_obs(atac, 'n_genes_by_counts', lambda x: (x >= 2000) & (x <= 15000))
mu.pp.filter_obs(atac, 'total_counts', lambda x: (x >= 4000) & (x <= 40000))

## Normalize
atac.layers["counts"] = atac.X
ac.pp.tfidf(atac, scale_factor=1e4)
sc.pp.normalize_per_cell(atac, counts_per_cell_after=1e4)
sc.pp.log1p(atac)

## HVP
sc.pp.highly_variable_genes(atac, min_mean=0.05, max_mean=1.5, min_disp=.5)

## LSI
atac.raw = atac
ac.tl.lsi(atac)
atac.obsm['X_lsi'] = atac.obsm['X_lsi'][:,1:]
atac.varm["LSI"] = atac.varm["LSI"][:,1:]
atac.uns["lsi"]["stdev"] = atac.uns["lsi"]["stdev"][1:]

## UMAP
sc.pp.neighbors(atac, use_rep="X_lsi", n_neighbors=10, n_pcs=30)
sc.tl.leiden(atac, resolution=.5)
sc.tl.umap(atac, spread=1.5, min_dist=.5, random_state=20)

## Remove noisy clusters
mu.pp.filter_obs(atac, "leiden", lambda x: ~x.isin(["13", "14", "15"]))

## Annotate
new_cluster_names = {
    "0": "intermediate mono",
    "1": "CD4+ memory T",
    "2": "CD8+ naïve T",
    "3": "CD4+ naïve T",
    "4": "CD14 mono",
    "5": "CD8+ activated T",
    "6": "NK",
    "7": "memory B",
    "8": "CD16 mono",
    "9": "naïve B",
    "10": "mDC",
    "11": "MAIT",
    "12": "pDC",
}

atac.obs['celltype'] = atac.obs.leiden.astype("str").values
atac.obs.celltype = atac.obs.celltype.astype("category").cat.rename_categories(new_cluster_names)
atac.obs.celltype.cat.reorder_categories([
    'CD4+ naïve T', 'CD4+ memory T', 'MAIT',
    'CD8+ naïve T', 'CD8+ activated T', 'NK',
    'naïve B', 'memory B',
    'CD14 mono', 'intermediate mono', 'CD16 mono',
    'mDC', 'pDC'], inplace=True)

# Integration

## Intersect
mu.pp.intersect_obs(mdata)
mdata.update()

## MOFA
mu.tl.mofa(mdata, verbose=True)

## UMAP
sc.pp.neighbors(mdata, use_rep="X_mofa")
sc.tl.leiden(mdata, key_added='leiden_joint', resolution=0.5)
sc.tl.umap(mdata, min_dist=.2, spread=1., random_state=10)

## Annotate
new_cluster_names = {
    "0": "CD8+ naïve T",
    "1": "CD4+ naïve T",
    "2": "CD14 mono",
    "3": "CD4+ memory T",
    "4": "intermediate mono",
    "5": "CD4+ memory T",
    "6": "CD14 mono",
    "7": "CD8+ activated T",
    "8": "NK", 
    "9": "memory B",
    "10": "CD16 mono",
    "11": "naïve B",
    "12": "CD8+ activated T",
    "13": "mDC",
    "14": "MAIT",
    "15": "pDC"
}
mdata.obs['celltype'] = mdata.obs.leiden_joint.astype("str")
mdata.obs.celltype = mdata.obs.celltype.map(new_cluster_names).astype("category")
mdata.obs.celltype.cat.reorder_categories([
    'CD4+ naïve T', 'CD4+ memory T', 'MAIT',
    'CD8+ naïve T', 'CD8+ activated T', 'NK',
    'naïve B', 'memory B',
    'CD14 mono', 'intermediate mono', 'CD16 mono',
    'mDC', 'pDC'], inplace=True)

## Figure
fig = mu.pl.umap(mdata, color="celltype", legend_loc="on data", frameon=False, return_fig=True, show=False)
fig.set_facecolor('white')
fig.savefig(plot)

## Write
mdata.write(out_dir)