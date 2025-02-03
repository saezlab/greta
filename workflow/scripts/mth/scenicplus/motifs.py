from scenicplus.triplet_score import get_max_rank_of_motif_for_each_TF
from pycistarget.motif_enrichment_cistarget import cisTargetDatabase
from pycistarget.utils import load_motif_annotations
import scipy.sparse as ss
import numpy as np
import pandas as pd
import pyranges as pr
import anndata as ad
import sys


def get_pr(index):
    df = index.str.replace(':', '-') 
    df = df.str.split('-').tolist()
    df = pd.DataFrame(df, columns=['Chromosome', 'Start', 'End'])
    return pr.PyRanges(df)


def get_motifs_for_TF(tf_names, annotation_to_use, motif_to_tf):
    motif_to_tf = motif_to_tf.fillna("")[annotation_to_use].agg(", ".join, axis = 1).apply(lambda x: [x for x in x.split(", ") if len(x) > 0])
    motif_to_tf = motif_to_tf.loc[[len(x) > 0 for x in motif_to_tf]]
    tf_to_motif = motif_to_tf.explode().reset_index().drop_duplicates().groupby(0)["MotifID"].apply(lambda x: ','.join(list(x)))
    tf_names = pd.Index(tf_names)
    tf_names = tf_names.intersection(tf_to_motif.index)
    return tf_to_motif.loc[tf_names].to_dict()


path_tfb = sys.argv[1]
path_db = sys.argv[2]
path_out = sys.argv[3]

# Read
tfb = pd.read_csv(path_tfb)
tfb['cre'] = tfb['cre'].str.replace('-', ':', 1)
var_names = tfb['tf'].unique()
obs_names = tfb['cre'].unique()

# Create anndata
motifs = ad.AnnData(
    obs=pd.DataFrame(index=obs_names),
    var=pd.DataFrame(index=var_names),
    X=ss.lil_matrix((obs_names.size, var_names.size), dtype=bool),
)
for cre, tfs in tfb.groupby('cre')['tf'].apply(lambda x: np.array(x)).items():
    motifs[cre, tfs].X = True
motifs.X = ss.csr_matrix(motifs.X)

# Find motif annots
motif_to_tf = load_motif_annotations(
    specie = "homo_sapiens",
    version = "v10nr_clust",
    motif_similarity_fdr = 0.001,
    orthologous_identity_threshold = 0.0)

# Remove ann if they are not in db or cres not in db
ctx_db = cisTargetDatabase(
    fname=path_db,
    region_sets=get_pr(motifs.obs_names)
)
inter = motif_to_tf.index.intersection(ctx_db.db_rankings.index)
motif_to_tf = motif_to_tf.loc[inter]
motif_to_tf.index.name = 'MotifID'

# Find motif anns per tf gene name
tf_to_motif = get_motifs_for_TF(
    tf_names = motifs.var_names,
    annotation_to_use = ["Direct_annot", "Orthology_annot"],
    motif_to_tf = motif_to_tf
)

# Subset and add annots
m_msk = motifs.var_names.isin(tf_to_motif)
motifs = motifs[:, m_msk].copy()
motifs.var.loc[:, 'motifs'] = [tf_to_motif[v] for v in motifs.var_names]

# Remove regions not found in db
df = get_max_rank_of_motif_for_each_TF(motifs, path_db)
inter = motifs.obs_names.intersection(df.index)
motifs = motifs[inter, :].copy()

# Write
motifs.write(path_out)
