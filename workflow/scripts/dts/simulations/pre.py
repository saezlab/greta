import pandas as pd
import numpy as np
import anndata as ad
import mudata as mu
import scanpy as sc
import scipy.sparse as sps
import os
import glob
import sys

path_dataset = sys.argv[1]

path_gex = os.path.join(path_dataset, 'expression.csv')
path_acc = os.path.join(path_dataset, 'peaks.csv')
path_p2g = os.path.join(path_dataset, 'region_to_gene.csv')
path_tfb = os.path.join(path_dataset, 'region_to_tf.csv')

rna = ad.AnnData(pd.read_csv(path_gex, index_col=0).T)
rna.var.index.name = None
atac = ad.AnnData(pd.read_csv(path_acc, index_col=0).T)
atac.X = np.ceil(atac.X)
atac.var.index.name = None

def transform_cre(r):
    r = int(r.replace('region', ''))
    return f'chrZ-{r}-{r + 500}'

atac.var_names = [transform_cre(r) for r in atac.var_names]

rna.X = sps.csr_matrix(rna.X)
atac.X = sps.csr_matrix(atac.X)

sc.pp.filter_cells(rna, min_genes=200)
sc.pp.filter_genes(rna, min_cells=3)
sc.pp.filter_cells(atac, min_genes=200)
sc.pp.filter_genes(atac, min_cells=3)
inter = rna.obs_names.intersection(atac.obs_names)
rna = rna[inter, :].copy()
atac = atac[inter, :].copy()

sc.pp.normalize_total(rna, target_sum=1e4)
sc.pp.log1p(rna)
sc.pp.normalize_total(atac, target_sum=1e4)
sc.pp.log1p(atac)

del rna.obs
del rna.var
del rna.uns
del atac.obs
del atac.var
del atac.uns

obs = rna.obs.copy()
obs['batch'] = 'smpl'
obs['celltype'] = 'ctype'

mdata = mu.MuData(
    {'rna': rna, 'atac': atac,},
    obs=obs
)

p2g = pd.read_csv(path_p2g)
p2g = pd.melt(p2g, id_vars=['region'], var_name='gene', value_name='score')
p2g = p2g[p2g['score'] != 0]
p2g['region'] = [transform_cre(r) for r in p2g['region']]
p2g = p2g.rename(columns={'region': 'cre'}).assign(pval=0.01).reset_index(drop=True)

tfb = pd.read_csv(path_tfb)
tfb = pd.melt(tfb, id_vars=['region'], var_name='tf', value_name='score')
tfb = tfb[tfb['score'] != 0]
tfb['region'] = [transform_cre(r) for r in tfb['region']]
tfb = tfb.rename(columns={'region': 'cre'}).assign(pval=0.01).reset_index(drop=True)

# Delete tmp files:
for path_f  in glob.glob(os.path.join(path_dataset, '*.csv')):
    os.remove(path_f)

# Create outputs
mdata.write(os.path.join(path_dataset, 'mdata.h5mu'))
p2g.to_csv(os.path.join(path_dataset, 'p2g.csv'), index=False)
tfb.to_csv(os.path.join(path_dataset, 'tfb.csv'), index=False)
