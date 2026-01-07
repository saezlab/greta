import numpy as np
import pandas as pd
import mudata as mu
import pyranges as pr
import scipy.stats as sts
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_data', required=True)
parser.add_argument('-t','--path_tf', required=True)
parser.add_argument('-p','--path_proms')
parser.add_argument('-m','--mode', required=True)
parser.add_argument('-r','--thr_r2', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_data = args['path_data']
path_tf = args['path_tf']
path_proms = args['path_proms']
mode = args['mode']
thr_r2 = float(args['thr_r2'])
path_out = args['path_out']

# Read
mdata = mu.read(path_data)
genes = mdata.mod['rna'].var_names.astype('U')
peaks = mdata.mod['atac'].var_names.astype('U')
tfs = set(pd.read_csv(path_tf, header=None).loc[:, 0].values.astype('U'))
tfs = np.array(list(set(genes) & tfs))

# Compute grn
x, y = mdata.mod['rna'][:, tfs].X, mdata.mod['rna'].X
if mode == 'spearman':
    x = sts.rankdata(x, axis=0)
    y = sts.rankdata(y, axis=0)
grn = np.corrcoef(x=x, y=y, rowvar=False)
grn = pd.DataFrame(grn[:tfs.size, tfs.size:], index=tfs, columns=genes)
grn = grn.reset_index(names='source').melt(id_vars='source', var_name='target', value_name='score')

# Filter by r2 and self regulation
grn = grn[grn['score'].abs() > thr_r2]
grn = grn[grn['source'] != grn['target']]

# Transform peaks
peaks = pd.DataFrame(peaks, columns=['cre'])
peaks[['Chromosome', 'Start', 'End']] = peaks['cre'].str.split('-', n=2, expand=True)
peaks = pr.PyRanges(peaks[['Chromosome', 'Start', 'End']])

# Filter by genes
if path_proms:
    proms = pr.read_bed(path_proms)
    proms = proms[proms.Name.astype('U').isin(genes)]
    
    # Filter by peaks
    proms = proms.nearest(peaks)
    proms = proms.df[proms.df['Distance'] == 0]
    proms['cre'] = proms['Chromosome'].astype(str) + '-' + proms['Start_b'].astype(str) + '-' + proms['End_b'].astype(str)
    proms = proms[['cre', 'Name']].rename(columns={'Name': 'target'}).drop_duplicates()
    
    # Merge
    grn = pd.merge(grn, proms, how='inner')[['source', 'cre', 'target', 'score']]
    grn = grn.sort_values(['source', 'target', 'cre'])

# Filter regulons with less than 5 targets
n_targets = grn.groupby(['source']).size().reset_index(name='counts')
n_targets = n_targets[n_targets['counts'] > 5]
grn = grn[grn['source'].isin(n_targets['source'])]

# Write
grn.to_csv(path_out, index=False)
