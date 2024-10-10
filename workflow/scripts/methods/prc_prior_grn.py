import pyranges as pr
import pandas as pd
import numpy as np
import os
import mudata as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--grn_path', required=True)
parser.add_argument('-d','--data_path', required=True)
parser.add_argument('-p','--proms_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
data_path = args['data_path']
proms_path = args['proms_path']
out_path = args['out_path']

# Read
grn = pd.read_csv(grn_path)
genes = mu.read(os.path.join(data_path, 'mod', 'rna')).var_names.astype('U')
peaks = mu.read(os.path.join(data_path, 'mod', 'atac')).var_names.astype('U')
proms = pr.read_bed(proms_path)

# Transform peaks
peaks = pd.DataFrame(peaks, columns=['cre'])
peaks[['Chromosome', 'Start', 'End']] = peaks['cre'].str.split('-', n=2, expand=True)
peaks = pr.PyRanges(peaks[['Chromosome', 'Start', 'End']])

# Filter by genes
grn = grn[grn['source'].astype('U').isin(genes) & grn['target'].astype('U').isin(genes)]
proms = proms[proms.Name.astype('U').isin(genes)]

# Filter by peaks
proms = proms.overlap(peaks)
proms.cre = proms.df['Chromosome'].astype(str) + '-' + proms.df['Start'].astype(str) + '-' + proms.df['End'].astype(str)
proms = proms.df[['cre', 'Name']].rename(columns={'Name': 'target'})

# Merge
grn = pd.merge(grn, proms, how='inner')[['source', 'cre', 'target', 'weight']]
grn = grn.sort_values(['source', 'target', 'cre']).rename(columns={'weight': 'score'})

# Filter regulons with less than 5 targets
n_targets = grn.groupby(['source']).size().reset_index(name='counts')
n_targets = n_targets[n_targets['counts'] > 5]
grn = grn[grn['source'].isin(n_targets['source'])]

# Write
grn.to_csv(out_path, index=False)
