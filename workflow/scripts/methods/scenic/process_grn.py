import pyranges as pr
import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--grn_path', required=True)
parser.add_argument('-p','--proms_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
proms_path = args['proms_path']
out_path = args['out_path']

# Read
grn = pd.read_csv(grn_path, index_col=False, sep='\t').rename(columns={'TF': 'source', 'importance': 'score'})
proms = pr.read_bed(proms_path)
proms['cre'] = proms.df['Chromosome'].astype(str) + '-' + proms.df['Start'].astype(str) + '-' + proms.df['End'].astype(str)
proms = proms.df[['cre', 'Name']].rename(columns={'Name': 'target'})

# Merge
grn = pd.merge(grn, proms, how='inner')[['source', 'cre', 'target', 'importance']]
grn = grn.sort_values(['source', 'target', 'cre'])

# Write
grn.to_csv(out_path, index=False)
