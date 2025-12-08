import h5py
import scipy.sparse as sps
import pandas as pd
import numpy as np
import anndata as ad
import scanpy as sc
import scglue
from itertools import chain
import os
import sys

path_gex = sys.argv[1]
path_acc = sys.argv[2]
path_gid = sys.argv[3]
path_ann = sys.argv[4]
path_out = sys.argv[5]

def read_h5ad(path_data):
    with h5py.File(path_data, 'r') as f:
        obs_names = f['matrix']['barcodes'][:].astype(str)
        var_names = f['matrix']['features']['name'][:].astype(str)
        obs = pd.DataFrame(index=obs_names)
        var = pd.DataFrame(index=var_names)
        X = sps.csr_matrix((f['matrix']['data'][:], f['matrix']['indices'][:], f['matrix']['indptr'][:]))
    return ad.AnnData(X=X, obs=obs, var=var)

rna = sc.read_10x_h5(path_gex)
atac = read_h5ad(path_acc)

# Process rna
geneids = pd.read_csv(path_gid).set_index('symbol')['id'].to_dict()
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
del rna.var
# Filter rna by annotation
obs = pd.read_csv(path_ann, index_col=0).set_index('X')
obs.index.name = None
rna = rna[obs.index, :].copy()
rna.layers["counts"] = rna.X.copy()
sc.pp.highly_variable_genes(rna, n_top_genes=2000, flavor="seurat_v3")
sc.pp.normalize_total(rna)
sc.pp.log1p(rna)
sc.pp.scale(rna)
sc.tl.pca(rna, n_comps=100, svd_solver="auto")

sc.pp.filter_cells(atac, min_genes=100)
sc.pp.filter_genes(atac, min_cells=3)
scglue.data.lsi(atac, n_components=100, n_iter=15)

import urllib.request
url = "ftp://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_38/gencode.v38.primary_assembly.annotation.gtf.gz"
path_gtf = os.path.join(os.path.dirname(path_gex), 'gencode.gtf.gz')
urllib.request.urlretrieve(url, path_gtf)
scglue.data.get_gene_annotation(
    rna, gtf=path_gtf,
    gtf_by="gene_name"
)
msk_rna = rna.var.loc[:, ["chrom", "chromStart", "chromEnd"]].dropna().index
rna = rna[:, msk_rna].copy()
os.remove(path_gtf)

split = atac.var_names.str.split(r"[:-]")
atac.var["chrom"] = split.map(lambda x: x[0])
atac.var["chromStart"] = split.map(lambda x: x[1]).astype(int)
atac.var["chromEnd"] = split.map(lambda x: x[2]).astype(int)

guidance = scglue.genomics.rna_anchored_guidance_graph(rna, atac)
scglue.graph.check_graph(guidance, [rna, atac])

scglue.models.configure_dataset(
    rna, "NB", use_highly_variable=True,
    use_layer="counts", use_rep="X_pca"
)
scglue.models.configure_dataset(
    atac, "NB", use_highly_variable=True,
    use_rep="X_lsi"
)
guidance_hvf = guidance.subgraph(chain(
    rna.var.query("highly_variable").index,
    atac.var.query("highly_variable").index
)).copy()

glue = scglue.models.fit_SCGLUE(
    {"rna": rna, "atac": atac}, guidance_hvf,
    fit_kws={"directory": os.path.dirname(path_gex)}
)

X_rna = glue.encode_data("rna", rna)
X_rna = pd.DataFrame(X_rna, index=['rna_' + o.split('-')[0] for o in rna.obs_names])
X_atac = glue.encode_data("atac", atac)
X_atac = pd.DataFrame(X_atac, index=['atac_' + o.split('-')[0] for o in atac.obs_names])
X = pd.concat([X_rna, X_atac])

X.to_csv(path_out)
