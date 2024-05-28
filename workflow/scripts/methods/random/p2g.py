import mudata as mu
import pandas as pd
import numpy as np
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
out_path = args['out_path']

g_perc = 0.05
n_cres = 3

# Read features
genes = mu.read(inp_path + 'mod/rna').var_names
cres = mu.read(inp_path + 'mod/atac').var_names

# Randomly sample genes
rng = np.random.default_rng(seed=42)
n = int(np.round(genes.size * g_perc))
genes = rng.choice(genes, n, replace=False)

# Randomly sample peak-gene connections
df = []
for g in genes:
    r_cres = rng.choice(cres, n_cres, replace=False)
    for cre in r_cres:
        df.append([cre, g, 1])
df = pd.DataFrame(df, columns=['cre', 'gene', 'score'])
df = df.sort_values(['cre', 'gene']).drop_duplicates(['cre', 'gene'])

# Write
df.to_csv(out_path, index=False)
 