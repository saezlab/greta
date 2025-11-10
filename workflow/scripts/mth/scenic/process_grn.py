import pyranges as pr
import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--grn_path', required=True)
parser.add_argument('-p','--proms_path')
parser.add_argument('-r','--reg_path')
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
proms_path = args['proms_path']
out_path = args['out_path']
reg_path = args['reg_path']

# Read
grn = pd.read_csv(grn_path, index_col=False, sep='\t').rename(columns={'TF': 'source', 'importance': 'score'})
#grn = grn[grn["score"] > 0.001]
top_n = int(np.ceil(grn.shape[0] * 0.1))
grn = grn.sort_values('score', ascending=False).head(top_n)

if proms_path:
    # Merge
    proms = pr.read_bed(proms_path).df
    proms['cre'] = proms['Chromosome'].astype(str) + '-' + proms['Start'].astype(str) + '-' + proms['End'].astype(str)
    proms = proms[['cre', 'Name']].rename(columns={'Name': 'target'})
    grn = pd.merge(grn, proms, on='target', how='inner')
    cols = ['source', 'cre', 'target', 'score']
else:
    cols = ['source', 'target', 'score']

# Filter by enriched TFs
if reg_path:
    reg = pd.read_csv(reg_path)
    reg = reg.iloc[2:, [0, 8]]
    reg.columns = ['source', 'target']
    reg['target'] = reg['target'].str.split(',')
    reg_exp = reg.explode('target')
    reg_exp['target'] = reg_exp['target'].str.replace(r"[\[\(\)' ]", "", regex=True)
    grn = pd.merge(grn, reg_exp, on=['source', 'target'], how='inner')

# Format
grn = grn[cols]
grn = grn.sort_values(['source', 'target'])

# Write
grn.to_csv(out_path, index=False)
