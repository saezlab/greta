import pandas as pd
import pyranges as pr
import numpy as np
import sys

path_db = sys.argv[1]  #'dbs/hg38/c2g/eqtlcatalogue/eqtlcatalogue.bed.gz'
path_tss = sys.argv[2]  #'dbs/hg38/cre/promoters/promoters.bed.gz'
path_out = sys.argv[3]

db = pr.read_bed(path_db).df
tss = pr.read_bed(path_tss).df

# Get TSS midpoint per gene
tss['TSS_mid'] = (tss['Start'] + tss['End']) // 2
tss_mid = tss.groupby('Name')['TSS_mid'].first().reset_index()

# Merge and compute distance
res = db.merge(tss_mid, on='Name', how='left')
res['Distance'] = ((res['Start'] + res['End']) // 2 - res['TSS_mid']).abs()

# Filter
res = res[res['Distance'] >= 5_000]
res = res[['Chromosome', 'Start', 'End', 'Name', 'Score']]

# Write
res.to_csv(path_out, index=False, sep='\t', header=None, compression="gzip")
