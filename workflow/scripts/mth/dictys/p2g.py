import argparse, os, sys
import pandas as pd
import numpy as np
import mudata as md


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--path_data', required=True)
parser.add_argument('-t', '--tmp_path', required=True)
parser.add_argument('-p', '--path_out', required=True)
parser.add_argument('-g', '--gene_annotation', required=True)
parser.add_argument('-e', '--ext', required=True)
args = vars(parser.parse_args())

path_data = args['path_data']
tmp_path = args['tmp_path']
annot = args['gene_annotation']
path_out = args['path_out']
distance = int(args['ext'])


# Write the RNA matrix and ATAC matrix to working directory
rna_filename = os.path.join(tmp_path, "expression.tsv.gz")
atac_filename = os.path.join(tmp_path, "atac_peak.tsv.gz")
dist_filename = os.path.join(tmp_path, "tssdist.tsv.gz")
data = md.read(path_data)
rna_X = pd.DataFrame(np.array(data['rna'].X).T, columns=data['rna'].obs.index, index=data['rna'].var.index)
rna_X.to_csv(rna_filename, sep="\t", compression="gzip")

atac_peak_names = [n.replace('-', ':') for n in data['atac'].var.index]
atac_X = pd.DataFrame(np.zeros((data['atac'].var.index.shape[0], 1)), index=atac_peak_names, columns=['placeholder'])
atac_X.to_csv(atac_filename, sep="\t", compression="gzip")

# Identify all peaks that are within Xbp of annotated TSS
os.system(f'python3 -m dictys chromatin tssdist --cut {distance} {rna_filename} {atac_filename} {annot} {dist_filename}')

# Convert distance to score for p2g
df = pd.read_csv(dist_filename, sep='\t').rename(columns={'region': 'cre', 'target': 'gene', 'dist': 'score'})
df['score'] = -np.abs(df['score'])
df['cre'] = df['cre'].str.replace(':', '-')
df = df.sort_values('score', ascending=False).reset_index(drop=True).reset_index(names='rank')
df['score'] = (1 - (df['rank'] / df['rank'].max()))
df[['cre', 'gene', 'score']].to_csv(path_out, index=False)
