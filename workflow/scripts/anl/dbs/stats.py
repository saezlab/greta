import pandas as pd
import re
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--paths_prt', required=True, nargs='+')
parser.add_argument('-g','--paths_gst', required=True, nargs='+')
parser.add_argument('-m','--paths_tfm', required=True, nargs='+')
parser.add_argument('-t','--paths_tfp', required=True, nargs='+')
parser.add_argument('-b','--paths_tfb', required=True, nargs='+')
parser.add_argument('-c','--paths_cre', required=True, nargs='+')
parser.add_argument('-e','--paths_c2g', required=True, nargs='+')
parser.add_argument('-o','--path_out', required=True)
args = parser.parse_args()

df = []

# prt
for path in args.paths_prt:
    name = re.search(r'prt/([^/]+)/', path).group(1)
    meta = pd.read_csv(path)
    n_tfs = meta['TF'].unique().size
    n_exp = meta.shape[0]
    df.append(['prt', 'tfs', name, n_tfs])
    df.append(['prt', 'exp', name, n_exp])

df.append(['sss', 'TFs', 'dataset', 1000]) # Around 1k TFs (Lambert)
df.append(['sss', 'celltypes', 'dataset', 15]) # Around 15 cell types per dataset
df.append(['omc', 'genes', 'dataset', 16384])
df.append(['omc', 'cres', 'dataset', 65536])

# gst
for path in args.paths_gst:
    name = re.search(r'gst/([^/]+)\.csv$', path).group(1)
    gst = pd.read_csv(path)
    n_gst = gst['source'].unique().size
    n_gns = gst['target'].unique().size
    df.append(['gst', 'gst', name, n_gst])
    df.append(['gst', 'gns', name, n_gns])

def get_cats(col_cat):
    cats = set()
    for c in col_cat:
        cats.update(set(c.split(',')))
    return cats

# tfm
for path in args.paths_tfm:
    name = re.search(r'tfm/([^/]+)/', path).group(1)
    tfm = pd.read_csv(path, sep='\t', header=None)
    n_tfs = tfm[0].unique().size
    n_cats = len(get_cats(tfm[1]))
    df.append(['tfm', 'tfs', name, n_tfs])
    df.append(['tfm', 'cat', name, n_cats])


# tfp
for path in args.paths_tfp:
    name = re.search(r'tfp/([^/]+)/', path).group(1)
    tfp = pd.read_csv(path, sep='\t', header=None)
    n_tfs = len(set(tfp[0]) | set(tfp[1]))
    n_prs = len(tfp)
    df.append(['tfp', 'tfs', name, n_tfs])
    df.append(['tfp', 'prs', name, n_prs])

# tfb
for path in args.paths_tfb:
    name = re.search(r'tfb/([^/]+)/', path).group(1)
    tfb = pd.read_csv(path, sep='\t', header=None)
    n_tfs = tfb[3].unique().size
    n_cats = len(get_cats(tfb[4]))
    df.append(['tfb', 'tfs', name, n_tfs])
    df.append(['tfb', 'cat', name, n_cats])

# cre
for path in args.paths_cre:
    name = re.search(r'cre/([^/]+)/', path).group(1)
    cre = pd.read_csv(path, sep='\t', header=None)
    n_cre = cre.shape[0]
    nbp = (1 + cre[2] - cre[1]).sum()
    df.append(['cre', 'ncr', name, n_cre])
    df.append(['cre', 'nbp', name, nbp])

# c2g
for path in args.paths_c2g:
    name = re.search(r'c2g/([^/]+)/', path).group(1)
    c2g = pd.read_csv(path, sep='\t', header=None)
    n_gns = c2g[3].unique().size
    n_cats = len(get_cats(tfb[4]))
    df.append(['c2g', 'gns', name, n_gns])
    df.append(['c2g', 'cat', name, n_cats])

# Merge
df = pd.DataFrame(df, columns=['metric', 'type', 'name', 'val'])

# Write
df.to_csv(args.path_out, index=False)
