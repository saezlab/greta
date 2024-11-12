import matplotlib.pyplot as plt  # Import else compiling error
import pandas as pd
import numpy as np
from celloracle import motif_analysis as ma
import celloracle as co
import mudata as mu
import os
import re
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--path_data', required=True)
parser.add_argument('-a','--all_peaks', required=True)
parser.add_argument('-c','--connections', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-t','--thr', required=True)
parser.add_argument('-p','--path_out', required=True)
args = vars(parser.parse_args())

path_data = args['path_data']
path_all_peaks = args['all_peaks']
path_connections = args['connections']
organism = args['organism']
thr_coaccess = float(args['thr'])
path_out = args['path_out']

# Process organism
organism = re.search(r'^dbs/([^/]+)/.*$', organism).group(1)

# Load scATAC-seq peak list
peaks = pd.read_csv(path_all_peaks, index_col=0).x.values.astype('U')
peaks = np.char.replace(peaks, '-', '_')

# Load Cicero coaccessibility scores
cicero_connections = pd.read_csv(path_connections, index_col=0)
cicero_connections['Peak1'] = np.char.replace(cicero_connections['Peak1'].values.astype('U'), '-', '_')
cicero_connections['Peak2'] = np.char.replace(cicero_connections['Peak2'].values.astype('U'), '-', '_')

# Extract tss information
tss_annotated = ma.get_tss_info(
    peak_str_list=peaks,
    ref_genome=organism
)

# Integrate
integrated = ma.integrate_tss_peak_with_cicero(
    tss_peak=tss_annotated,
    cicero_connections=cicero_connections
)

# Process
integrated = integrated[integrated['coaccess'] >= thr_coaccess]
integrated['peak_id'] = integrated['peak_id'].str.replace('_', '-')
integrated = integrated.rename(columns={'peak_id': 'cre', 'gene_short_name': 'gene', 'coaccess': 'score'})
integrated = integrated.sort_values(['cre', 'score'], ascending=[True, False])

# Remove unexpressed genes
genes = mu.read(os.path.join(path_data, 'rna')).var.index.values.astype('U')
integrated = integrated[integrated['gene'].isin(genes)]

# Write
integrated.to_csv(path_out, index=False)
