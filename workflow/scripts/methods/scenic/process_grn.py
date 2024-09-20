import pyranges as pr
import pandas as pd
import numpy as np
import os
import mudata as mu
import argparse




# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--grn_path', required=True)
parser.add_argument('-i','--data', required=True)
parser.add_argument('-p','--proms_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
data_path = args['data']
proms_path = args['proms_path']
out_path = args['out_path']

# Read
grn = pd.read_csv(grn_path, index_col=False, sep='\t')
genes = mu.read(os.path.join(data_path, 'mod', 'rna')).var_names.astype('U')
proms = pr.read_bed(proms_path)

# Filter by genes
grn = grn[grn['TF'].astype('U').isin(genes) & grn['target'].astype('U').isin(genes)]
proms = proms[proms.Name.astype('U').isin(genes)]


proms['cre'] = proms.df['Chromosome'].astype(str) + '-' + proms.df['Start'].astype(str) + '-' + proms.df['End'].astype(str)
proms = proms.df[['cre', 'Name']].rename(columns={'Name': 'target'})
grn.rename(columns={'TF': 'source'}, inplace=True)




# Merge
grn = pd.merge(grn, proms, how='inner')[['source', 'cre', 'target', 'importance']]
grn = grn.sort_values(['source', 'target', 'cre']).rename(columns={'importance': 'score'})

# Write
grn.to_csv(out_path, index=False)



