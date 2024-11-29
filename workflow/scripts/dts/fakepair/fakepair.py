import pandas as pd
import mudata as mu
import argparse
import sys
import os


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-b','--path_barmap', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
path_barmap = args['path_barmap']
path_output = args['path_output']

# Read
mdata = mu.read('dts/pitupair/annotated.h5mu')
barmap = pd.read_csv('dts/fakepitupair/barmap.csv')

# Format RNA barmap
barmap.loc[:, 'RNA'] = ['smpl_' + b.replace('-1', '') for b in barmap['RNA']]

# Make sure intersection of all
inter = set(barmap['RNA']) & set(barmap['ATAC']) & set(mdata.obs_names)
msk = barmap['ATAC'].isin(inter) & barmap['RNA'].isin(inter)
barmap = barmap.loc[msk, :].reset_index(drop=True)
mdata = mdata[list(inter), :].copy()

# Create new fake object
fmdata = mdata[barmap['ATAC'], :].copy()

# Populate with predicted RNA
fmdata.mod['rna'].X = mdata.mod['rna'][barmap['RNA'].values, :].X

# Update metadata
obs = barmap.set_index('ATAC')
obs.index.name = None
fmdata.obs = obs

# Write
fmdata.write(path_output)
