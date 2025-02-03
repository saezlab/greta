import scipy.sparse as scs
import anndata as ad
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

# Find shared regions
mtf_pr = get_pr(motifs.obs_names)
p2g_pr = get_pr(pd.Index(p2g['cre'].unique()))
inter = mtf_pr.join(p2g_pr)
inter_motifs = get_vars(inter[['Chromosome', 'Start', 'End']])
inter_p2g = get_vars(pr.PyRanges(inter.df[['Chromosome', 'Start_b', 'End_b']].rename(columns={'Start_b': 'Start', 'End_b': 'End'})))

# Create matching motif anndata
new_motifs = ad.AnnData(
    var=pd.DataFrame(index=motifs.var_names),
    obs=pd.DataFrame(index=inter_p2g),
    X=scs.csr_matrix((inter_p2g.size, motifs.var_names.size))
)
new_motifs[inter_p2g, :].X = motifs[inter_motifs, :].X

# Build df
new_motifs.X = new_motifs.X.tocoo()
tfb = pd.DataFrame()
tfb['cre'] = new_motifs.obs_names[new_motifs.X.row]
tfb['tf'] = new_motifs.var_names[new_motifs.X.col]
tfb['score'] = 5.
tfb['cre'] = tfb['cre'].str.replace(':', '-')

# Write
tfb.to_csv(path_out, index=False)
