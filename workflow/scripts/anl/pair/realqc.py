import mudata as mu
import pandas as pd
import numpy as np
import scanpy as sc
import sys


pmdata = mu.read(sys.argv[1])
nmdata = mu.read(sys.argv[2])


def get_qc_omic(mdata, omic, tpe):
    adata = mdata.mod[omic]
    adata.X = adata.layers['counts']
    obs, _ = sc.pp.calculate_qc_metrics(
        adata, percent_top=None, log1p=True
    )
    qc = obs.assign(omic=omic, type=t)
    qc = pd.merge(qc.reset_index(names='barcode'), mdata.obs.reset_index(names='barcode')[['barcode', 'celltype']], on=['barcode'], how='inner')
    return qc


def extract_n_cells(mdata, tpe):
    return mdata.obs.groupby('celltype', as_index=False).size().sort_values('celltype').assign(type=tpe)


# Compute qc
types = ['paired', 'upaired']
omics = ['rna', 'atac']
qc = []
n_ctps = []
for mdata, t in zip([pmdata, nmdata], types):
    n_ctps.append(extract_n_cells(mdata, t))
    for omic in omics:
        qc.append(get_qc_omic(mdata, omic, t))
qc = pd.concat(qc)
n_ctps = pd.concat(n_ctps)

# Write
qc.to_csv(sys.argv[3], index=False)
n_ctps.to_csv(sys.argv[4], index=False)
