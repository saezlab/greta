import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
import sys


mdata = mu.read(sys.argv[1])


def get_qc_omic(mdata, omic):
    adata = mdata.mod[omic]
    adata.X = adata.layers['counts']
    obs, _ = sc.pp.calculate_qc_metrics(
        adata, percent_top=None, log1p=True
    )
    qc = obs.assign(omic=omic)
    qc = pd.merge(qc.reset_index(names='barcode'), mdata.obs.reset_index(names='barcode')[['barcode', 'celltype']], on=['barcode'], how='inner')
    return qc


def extract_n_cells(mdata):
    return mdata.obs.groupby('celltype', as_index=False).size().sort_values('celltype')


# Compute qc
omics = ['rna', 'atac']
n_ctps = extract_n_cells(mdata)
qc = []
for omic in omics:
    qc.append(get_qc_omic(mdata, omic))
qc = pd.concat(qc)

# Write
qc.to_csv(sys.argv[2], index=False)
n_ctps.to_csv(sys.argv[3], index=False)
