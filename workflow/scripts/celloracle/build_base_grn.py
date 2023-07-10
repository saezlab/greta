import pandas as pd
import numpy as np
import celloracle as co
from celloracle import motif_analysis as ma
from celloracle.utility import save_as_pickled_object
import h5py
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_tfi', required=True)
parser.add_argument('-t','--thr_score', required=True)
parser.add_argument('-g','--path_grn', required=True)
args = vars(parser.parse_args())

path_tfi = args['path_tfi']
thr_score = float(args['thr_score'])
path_grn = args['path_grn']

# Read tfi
tfi = ma.load_TFinfo(path_tfi)

# Reset filtering
tfi.reset_filtering()

# Do filtering
tfi.filter_motifs_by_score(threshold=thr_score)

# Format post-filtering results.
tfi.make_TFinfo_dataframe_and_dictionary(verbose=True)

# Get final base GRN
df = tfi.to_dataframe()

# Save
df.to_csv(path_grn)
