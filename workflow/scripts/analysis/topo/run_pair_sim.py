import concurrent.futures
import pandas as pd
import numpy as np
import os
from tqdm import tqdm
from ..utils import ocoeff, parallel_ocoeff
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--paths', required=True, nargs='+')
parser.add_argument('-m','--mthds', required=True, nargs='+')
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

paths = args['paths']
mthds = args['mthds']
out_path = args['out_path']

print('Reading grns...')
names = []
dfs = []
for path in tqdm(paths):
    df = pd.read_csv(path).drop_duplicates(['source', 'target'], keep='first')
    dfs.append(df)
    names.append(os.path.basename(path).replace('.grn.csv', ''))

print('Computing pairwise overlap coefficients...')
index_pairs = [(i, j) for i in range(len(names)) for j in range(i + 1, len(names))]
names_a = []
names_b = []
tf_coefs = []
edge_coefs = []
target_coefs = []
with concurrent.futures.ProcessPoolExecutor() as executor:
    for res in tqdm(executor.map(parallel_ocoeff, index_pairs), total=len(index_pairs)):
        i, j, tf_coef, edge_coef, target_coef = res
        names_a.append(names[i])
        names_b.append(names[j])
        tf_coefs.append(tf_coef)
        edge_coefs.append(edge_coef)
        target_coefs.append(target_coef)

# Store as df
df = pd.DataFrame()
df['name_a'] = names_a
df['name_b'] = names_b
df['tf_oc'] = tf_coefs
df['edge_oc'] = edge_coefs
df['target_oc'] = target_coefs

# Write
df.to_csv(out_path, index=False)
