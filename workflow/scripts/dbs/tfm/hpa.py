import pandas as pd
import numpy as np
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-t','--tfs_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
tfs_path = args['tfs_path']
out_path = args['out_path']

# Read
c_cols = ['Tissue expression cluster', 'Cell line expression cluster', 'Single cell expression cluster']
df = pd.read_csv(inp_path, sep='\t').dropna(subset=c_cols)

# Read tfs
tfs = pd.read_csv(tfs_path, sep='\t', header=None).values.ravel().astype('U')

# Filter
msk_evd = (df['Evidence'] == 'Evidence at protein level').values
msk_loc = np.array(['Nucle' in str(s) for s in df['Subcellular location']])
msk_tissue =  ~df['Tissue expression cluster'].str.contains('Non-specific -')
msk_celine =  ~df['Cell line expression cluster'].str.contains('Non-specific -')
msk_cetype =  ~df['Single cell expression cluster'].str.contains('Non-specific -')
msk_tf = df['Gene'].isin(tfs)
msk = msk_tf & msk_evd & msk_loc & msk_tissue & msk_celine & msk_cetype
df = df.loc[msk, :]

# Format names
for col in c_cols:
    df[col] = [s.split(':')[1].split('-')[0].strip() for s in df[col]]
df['ctype'] = [','.join(sorted(set(lst))) for lst in df[c_cols].values]
df = df.rename(columns={'Gene': 'gene'})[['gene', 'ctype']]
df = df.sort_values(['gene', 'ctype'])

# Write
df.to_csv(out_path, sep='\t', index=False, header=None)
