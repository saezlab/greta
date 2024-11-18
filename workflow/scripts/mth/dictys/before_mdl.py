import pandas as pd
import numpy as np
import mudata as mu
import sys
import os
import dictys


# Read and process gex
rna_path = os.path.join(sys.argv[1], 'mod', 'rna')
rna = mu.read(rna_path)
rna.X = rna.layers['counts']
rna = rna.to_df().T
rna.to_csv(sys.argv[2], sep='\t', compression='gzip')
name_pre = sys.argv[1].split('/runs/')[1].split('.')[0]
if 'dictys' not in name_pre:
    dictys.preproc.qc_reads(sys.argv[2], sys.argv[2], 50, 10, 0, 200, 100, 0)
rna = pd.read_csv(sys.argv[2], header=0, index_col=0, sep='\t')

# Read and process peaks
use_peaks = bool(sys.argv[3])
if use_peaks:
    peaks = pd.read_csv(sys.argv[4])['cre'].unique()
else:
    atac_path = os.path.join(sys.argv[1], 'mod', 'atac')
    peaks = mu.read(atac_path).var_names
peaks = np.array([p.replace('-', ':') for p in peaks])
peaks = pd.DataFrame(np.zeros((peaks.size, 1)), index=peaks, columns=['placeholder'])
peaks.to_csv(sys.argv[5], sep='\t', compression='gzip')

# Read tfb
tfb = pd.read_csv(sys.argv[6])
tfb['cre'] = tfb['cre'].str.replace('-', ':')
tfb = tfb[tfb['tf'].isin(rna.index) & tfb['cre'].isin(peaks.index)]
output_tfb = tfb.rename(columns={'cre': 'loc', 'tf': 'TF'})[['TF', 'loc', 'score']]
output_tfb.to_csv(sys.argv[7], sep='\t', index=False)
