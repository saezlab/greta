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
    return obs.assign(omic=omic, type=t)


# Compute qc
types = ['paired', 'upaired']
omics = ['rna', 'atac']
qc = []
for mdata, t in zip([pmdata, nmdata], types):
    for omic in omics:
        qc.append(get_qc_omic(mdata, omic, t))
qc = pd.concat(qc)
qc = qc.reset_index(names='barcode')

# Write
qc.to_csv(sys.argv[3], index=False)
