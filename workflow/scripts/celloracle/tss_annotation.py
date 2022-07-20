import pandas as pd
import numpy as np
from celloracle import motif_analysis as ma
import celloracle as co
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--all_peaks', required=True)
parser.add_argument('-c','--connections', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-t','--thr', required=True)
parser.add_argument('-p','--prc_peaks', required=True)
args = vars(parser.parse_args())

path_all_peaks = args['all_peaks']
path_connections = args['connections']
organism = args['organism']
thr_coaccess = args['thr']
path_prc_peaks = args['prc_peaks']

# # Load scATAC-seq peak list
peaks = pd.read_csv(path_all_peaks, index_col=0).x.values
peaks = peaks.x.values

# Load Cicero coaccessibility scores
cicero_connections = pd.read_csv("cicero_connections.csv", index_col=0)
cicero_connections.head()

# Determine genome
if organism == 'human':
    genome = 'hg38'
elif organism == 'mouse':
    genome = 'mm10'

# Extract tss information
tss_annotated = ma.get_tss_info(peak_str_list=peaks, ref_genome=genome)

# Integrate
integrated = ma.integrate_tss_peak_with_cicero(tss_peak=tss_annotated,
                                               cicero_connections=cicero_connections)

# Filter by coaccessability threshold
peak = integrated[integrated.coaccess >= thr_coaccess]
peak = peak[["peak_id", "gene_short_name"]].reset_index(drop=True)

# Save
peak.to_csv(path_prc_peaks)

