import mudata as mu
import pandas as pd
import numpy as np
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-t','--tf_path', required=True)
parser.add_argument('-p','--p2g_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
tf_path = args['tf_path']
p2g_path = args['p2g_path']
out_path = args['out_path']
n_tfs_per_cre = 2

# Read genes and intersect with tfs
genes = mu.read(inp_path + 'mod/rna').var_names
tfs = pd.read_csv(tf_path, header=None).loc[:, 0].values.astype('U')
tfs = np.intersect1d(genes, tfs)

# Sample random tf-cre interactions
p2g = pd.read_csv(p2g_path)
cres = p2g.cre.unique().astype('U')
rng = np.random.default_rng(seed=42)
df = []

for cre in cres:
    r_tfs = rng.choice(tfs, n_tfs_per_cre)
    for tf in r_tfs:
        df.append([cre, tf, 1])
df = pd.DataFrame(df, columns=['cre', 'tf', 'score']).sort_values(['cre', 'tf'])

# Write
df.to_csv(out_path, index=False)
