import mudata as mu
import numpy as np
import pyranges as pr
import pandas as pd
import decoupler as dc
import scipy.stats as st
from tqdm import tqdm
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--pair_path', required=True)
parser.add_argument('-b','--npair_path', required=True)
parser.add_argument('-c','--cor_path', required=True)
parser.add_argument('-d','--sts_path', required=True)
args = parser.parse_args()


# Read
pair = mu.read(args.pair_path)
npair = mu.read(args.npair_path)

def extract_omic(mdata, omic):
    omic_adata = mdata.mod[omic].copy()
    omic_adata.obs = mdata.obs
    return omic_adata

def compute_corrs(pair, npair, omic):
    omic_pair = extract_omic(pair, omic)
    omic_npair = extract_omic(npair, omic)
    mean_pair = dc.get_pseudobulk(
        adata=omic_pair,
        sample_col='celltype',
        groups_col=None,
        mode='mean',
        min_cells=0,
        min_counts=0
    )
    mean_npair = dc.get_pseudobulk(
        adata=omic_npair,
        sample_col='celltype',
        groups_col=None,
        mode='mean',
        min_cells=0,
        min_counts=0
    )
    inter_ctypes = np.intersect1d(mean_pair.obs.index, mean_npair.obs.index)
    mean_pair = mean_pair[inter_ctypes, :].copy()
    mean_npair = mean_npair[inter_ctypes, :].copy()
    if omic == 'atac':
        pair_feat = pr.PyRanges(pd.DataFrame([p.split('-') for p in mean_pair.var_names], columns=['Chromosome', 'Start', 'End']))
        npair_feat = pr.PyRanges(pd.DataFrame([p.split('-') for p in mean_npair.var_names], columns=['Chromosome', 'Start', 'End']))
        overs = pair_feat.join(npair_feat).df
        pair_feat_str = (overs['Chromosome'].astype(str) + '-' + overs['Start'].astype(str) + '-' + overs['End'].astype(str)).values
        npair_feat_str = (overs['Chromosome'].astype(str) + '-' + overs['Start_b'].astype(str) + '-' + overs['End_b'].astype(str)).values
        mean_pair, mean_npair = mean_pair[:, pair_feat_str].copy(), mean_npair[:, npair_feat_str].copy()
        inter_size = pair_feat.overlap(npair_feat).df.shape[0]
        min_size =  np.min([pair_feat.df.shape[0], npair_feat.df.shape[0]])
    else:
        inter = np.intersect1d(mean_pair.var_names, mean_npair.var_names)
        mean_pair, mean_npair = mean_pair[:, inter].copy(), mean_npair[:, inter].copy()
        inter_size = inter.size
        min_size = np.min([mean_pair.var_names.size, mean_npair.var_names.size])
    ocoeff = inter_size / min_size
    df_cor = []
    for j in tqdm(range(mean_pair.shape[1])):
        var_name = mean_pair.var_names[j]
        x = mean_pair.X[:, j].ravel()
        y = mean_npair.X[:, j].ravel()
        stat, pval = st.spearmanr(x, y)
        df_cor.append(['var', omic, var_name, stat, pval])
    for i in range(mean_pair.shape[0]):
        obs_name = mean_pair.obs_names[i]
        x = mean_pair.X[i, :].ravel()
        y = mean_npair.X[i, :].ravel()
        stat, pval = st.spearmanr(x, y)
        df_cor.append(['obs', omic, obs_name, stat, pval])
    df_cor = pd.DataFrame(df_cor, columns=['type', 'omic', 'name', 'stat', 'pval'])
    
    return df_cor, pd.DataFrame([[omic, inter_size, min_size, ocoeff]], columns=['omic', 'inter', 'min_size', 'ocoeff'])

df_cor = []
df_sts = []
for omic in ['rna', 'atac']:
    cor, o_stats = compute_corrs(pair, npair, omic)
    df_cor.append(cor)
    df_sts.append(o_stats)
df_cor = pd.concat(df_cor)
df_sts = pd.concat(df_sts)

def compute_fdr_type(df_cor, typ):
    msk = df_cor['type'] == typ
    if 'padj' not in df_cor.columns:
        df_cor['padj'] = 1.
    df_cor.loc[msk, 'padj'] = st.false_discovery_control(df_cor.loc[msk, 'pval'])
compute_fdr_type(df_cor, typ='obs')
compute_fdr_type(df_cor, typ='var')

# Write
df_cor.to_csv(args.cor_path, index=False)
df_sts.to_csv(args.sts_path, index=False)
