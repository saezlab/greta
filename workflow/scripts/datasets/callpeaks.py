import os
import snapatac2 as snap
import anndata as ad
from snapatac2.datasets import _datasets, datasets
from pathlib import Path
import pandas as pd
import numpy as np
import scipy
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f','--path_frags', nargs='+', required=True)
parser.add_argument('-a','--path_annot', required=True)
parser.add_argument('-t','--path_tmp', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_frags = args['path_frags']
path_annot = args['path_annot']
path_tmp = args['path_tmp']
path_output = args['path_output']

print(path_frags)

if __name__ == '__main__':
    n_jobs = 16
    print('N of cores:', n_jobs)
    
    # Change default cache dir
    if not os.path.exists(path_tmp):
        os.mkdir(path_tmp)
    _datasets = datasets()
    _datasets.path = Path(path_tmp)
    
    # Find sample_ids
    sample_ids = [os.path.basename(p).split('.')[0] for p in path_frags]
    tmp_files = [os.path.join(path_tmp, p + '.frags.h5ad') for p in sample_ids]
    
    # Read and create h5ad fragment files
    _ = snap.pp.import_data(
        path_frags,
        chrom_sizes=snap.genome.hg38,
        file=tmp_files,
        tempdir=path_tmp,
        sorted_by_barcode=False,
    )
    del _

    # Filter by annotation
    annot = pd.read_csv(path_annot, index_col=0)
    annot.index = [sample_id + '_' + barcode.split('-1')[0] for barcode, sample_id in zip(annot.index, annot['batch'])]
    uns = None
    type_frags = None
    lst_obs = []
    lst_frags = []
    for sample_id, tmp_file in zip(sample_ids, tmp_files):
        tmp = ad.read_h5ad(tmp_file, backed='r')
        tmp.obs.index = [sample_id + '_' + barcode.split('-1')[0] for barcode in tmp.obs.index]
        # Filter by annotation
        obs = pd.merge(tmp.obs, annot, left_index=True, right_index=True)
        # Add uns
        uns = tmp.uns['reference_sequences']
        # Add type frags
        type_frags = list(tmp.obsm.keys())[0]
        # Add frags
        lst_frags.append(tmp[obs.index, :].obsm[type_frags])
        # Add obs
        lst_obs.append(obs)
    atac = ad.AnnData(obs=pd.concat(lst_obs))
    atac.uns = {'reference_sequences': uns}
    atac.obsm[type_frags] = scipy.sparse.vstack(lst_frags)
    
    # Call and merge peaks
    snap.tl.macs3(atac, groupby='celltype', replicate='batch', n_jobs=n_jobs, tempdir=path_tmp)
    peaks = snap.tl.merge_peaks(atac.uns['macs3'], snap.genome.hg38)
    atac = snap.pp.make_peak_matrix(atac, use_rep=peaks['Peaks'])
    
    # Clean
    del atac.obs
    del atac.var

    # Update format peaks
    new_var_names = []
    for p in atac.var_names:
        p = p.replace(':', '-')
        seq, start, end = p.split('-')
        end = int(end) - 1
        p = '{0}-{1}-{2}'.format(seq, start, end)
        new_var_names.append(p)
    atac.var_names = new_var_names
    
    # Write
    atac.write(path_output)
