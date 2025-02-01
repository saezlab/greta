import pyranges as pr
import h5py
import pandas as pd
import mudata as mu
import sys

def get_pr(index):
    df = index.str.replace(':', '-') 
    df = df.str.split('-').tolist()
    df = pd.DataFrame(df, columns=['Chromosome', 'Start', 'End'])
    return pr.PyRanges(df)


def get_vars(df):
    chrm = df.df['Chromosome'].astype(str)
    strt = df.df['Start'].astype(str)
    end = df.df['End'].astype(str)
    return pd.Index(chrm + ':' + strt + '-' + end)


path_pre = sys.argv[1]
path_p2g = sys.argv[2]
path_motifs = sys.argv[3]
path_out = sys.argv[4]

# Read
motifs = mu.read(path_motifs)
p2g = pd.read_csv(path_p2g)
if p2g.shape[0] == 0:
    tfb = pd.DataFrame(columns=['cre', 'tf', 'score'])
    tfb.to_csv(path_out, index=False)
    exit()

# Subset by tf genes
with h5py.File(path_pre, 'r') as f:
    genes = f['mod']['rna']['var']['_index'][:].astype('U')
tf_msk = motifs.var_names.isin(genes)
motifs = motifs[:, tf_msk]

# Subset by regions
mtf_pr = get_pr(motifs.obs_names)
p2g_pr = get_pr(pd.Index(p2g['cre'].unique()))
mtf_cres = get_vars(mtf_pr.overlap(p2g_pr))
motifs = motifs[mtf_cres, :].copy()
motifs.X = motifs.X.tocoo()

# Build df
tfb = pd.DataFrame()
tfb['cre'] = motifs.obs_names[motifs.X.row]
tfb['tf'] = motifs.var_names[motifs.X.col]
tfb['score'] = 5.
tfb['cre'] = tfb['cre'].str.replace(':', '-')

# Transform cres to pre's format
tfb_pr = get_pr(pd.Index(tfb['cre'].unique()))
over = tfb_pr.join(p2g_pr, how='left').df
over['cre'] = over['Chromosome'].astype(str) + '-' + over['Start'].astype(str) + '-' + over['End'].astype(str)
over['p2g_cre'] = over['Chromosome'].astype(str) + '-' + over['Start_b'].astype(str) + '-' + over['End_b'].astype(str)
over = over[['cre', 'p2g_cre']]
tfb = pd.merge(tfb, over, on='cre')
tfb = tfb[['p2g_cre', 'tf', 'score']].rename(columns={'p2g_cre': 'cre'})

# Write
tfb.to_csv(path_out, index=False)
