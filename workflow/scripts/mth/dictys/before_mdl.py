import pandas as pd
import numpy as np
import mudata as mu
import sys
import os
import dictys




# Read gex
rna_path = os.path.join(sys.argv[2], 'mod', 'rna')
rna = mu.read(rna_path)

# Modify the 'rna' object and save the result as a compressed file
rna.X = rna.layers['counts']
rna = rna.to_df().T
rna.to_csv(sys.argv[4], sep='\t', compression='gzip')

# Extract the prefix from the path and check if 'dictys' is not in the name
name_pre = sys.argv[2].split('/runs/')[1].split('.')[0]
if 'dictys' not in name_pre:
    dictys.preproc.qc_reads(sys.argv[4], sys.argv[4], 50, 10, 0, 200, 100, 0)
rna = pd.read_csv(sys.argv[4], header=0, index_col=0, sep='\t')

# Read tfb
tfb = pd.read_csv(sys.argv[1])
tfb['cre'] = tfb['cre'].str.replace('-', ':')

# Filter tfb for g in rna
tfb = tfb[tfb['tf'].isin(rna.index)]

# Store peaks as tsv
atac_path = os.path.join(sys.argv[2], 'mod', 'atac')
peaks = mu.read(atac_path).var_names
peaks = np.array([p.replace('-', ':') for p in peaks])
peaks = pd.DataFrame(np.zeros((peaks.size, 1)), index=peaks, columns=['placeholder'])
peaks.to_csv(sys.argv[3], sep='\t', compression='gzip')

# Rename columns and save the modified DataFrame
output_tfb = tfb.rename(columns={'cre': 'loc', 'tf': 'TF'})[['TF', 'loc', 'score']]
output_tfb.to_csv(sys.argv[5], sep='\t', index=False)
