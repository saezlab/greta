import mudata as mu
import numpy as np
import pandas as pd
import pyranges as pr
from tqdm import tqdm
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
parser.add_argument('-t','--tf_path', required=True)
parser.add_argument('-c','--cg_path', required=True)
parser.add_argument('-g','--g_perc', required=True)
parser.add_argument('-n','--scale', required=True)
parser.add_argument('-r','--tf_g_ratio', required=True)
parser.add_argument('-w','--w_size', required=True)
parser.add_argument('-s','--seed', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']
tf_path = args['tf_path']
cg_path = args['cg_path']
g_perc = float(args['g_perc'])
scale = float(args['scale'])
tf_g_ratio = float(args['tf_g_ratio'])
w_size = int(args['w_size'])
seed = int(args['seed'])
out_path = args['out_path']


def run_pre(mdata, seed):
    # Split
    atac = mdata.mod['atac'].copy()
    rna = mdata.mod['rna'].copy()
    # Shuffle genes
    rng = np.random.default_rng(seed=seed)
    genes = rng.choice(rna.var_names, rna.var_names.size, replace=False)
    rna.var_names = genes
    
    # Shuffle cres
    cres = rng.choice(atac.var_names, atac.var_names.size, replace=False)
    atac.var_names = cres
    
    # Shuffle obs
    obs = mdata.obs.copy()
    for col in obs.columns:
        obs[col] = rng.choice(obs[col], obs.shape[0], replace=False)
    mdata.obs = obs
    
    # Update
    mdata.mod['rna'] = rna
    mdata.mod['atac'] = atac

    return mdata


def get_window(gannot, target, w_size):
    target_gr = gannot[gannot.df['Name'] == target]
    tss = target_gr.df['Start'].values[0] + 1000  # Assumes window of 1000+/-
    tss_window = pr.from_dict({
        "Chromosome": target_gr.Chromosome,
        "Start": tss - w_size,
        "End": tss + w_size,
    })
    return tss_window


def get_cres_pr(mdata):
    cres = [c.split('-') for c in mdata.mod['atac'].var_names.values.astype('U')]
    cres = pr.PyRanges(pd.DataFrame(cres, columns=['Chromosome', 'Start', 'End']))
    return cres


def get_overlap_cres(gene, cres, w_size):
    wnd = get_window(gannot, target=gene, w_size=w_size)
    o_cres = cres.overlap(wnd)
    if len(o_cres) > 0:
        o_cres = (
            o_cres.df
            .assign(name=lambda x: x['Chromosome'].astype(str) + '-' + x['Start'].astype(str) + '-' + x['End'].astype(str))
            ['name'].values
        )
        return o_cres


def run_p2g(mdata, gannot, g_perc, scale, w_size, seed):
    # Read features
    genes = mdata.mod['rna'].var_names.values.astype('U')
    cres = get_cres_pr(mdata)
    g_in_ann = list(set(gannot.df['Name'].values))
    genes = genes[np.isin(genes, g_in_ann)]
    print(genes.size)
    
    # Randomly sample genes
    rng = np.random.default_rng(seed=seed)
    n = int(np.round(genes.size * g_perc))
    genes = rng.choice(genes, n, replace=False)
    n_cres = np.ceil(rng.exponential(scale=scale, size=genes.size))
    
    # Randomly sample peak-gene connections
    df = []
    for i in tqdm(range(genes.size)):
        n_cre = int(n_cres[i])
        g = genes[i]
        o_cres = get_overlap_cres(g, cres, w_size)
        if o_cres is not None:
            r_cres = rng.choice(o_cres, np.min([n_cre, o_cres.size]), replace=False)
            for cre in r_cres:
                df.append([cre, g, 1])
    df = pd.DataFrame(df, columns=['cre', 'gene', 'score'])
    df = df.sort_values(['cre', 'gene']).drop_duplicates(['cre', 'gene'])
    return df


def run_tfb(mdata, p2g, tfs, scale, seed):
    # Read genes and intersect with tfs
    genes = mdata.mod['rna'].var_names
    tfs = np.intersect1d(genes, tfs)

    # Sample random tf-cre interactions
    if p2g.shape[0] == 0:
        tfb = pd.DataFrame(columns=['cre', 'tf', 'score'])
        return tfb
    cres = p2g.cre.unique().astype('U')
    rng = np.random.default_rng(seed=seed)
    n_tfs_per_cres = np.ceil(rng.exponential(scale=scale, size=cres.size))
    df = []
    for i in range(cres.size):
        n_tfs_per_cre = int(n_tfs_per_cres[i])
        cre = cres[i]
        r_tfs = rng.choice(tfs, n_tfs_per_cre)
        for tf in r_tfs:
            df.append([cre, tf, 1])
    tfb = pd.DataFrame(df, columns=['cre', 'tf', 'score']).sort_values(['cre', 'tf']).drop_duplicates(['cre', 'tf'])
    return tfb


def run_mdl(p2g, tfb, tf_g_ratio, seed):
    if (p2g.shape[0] == 0) or (tfb.shape[0] == 0):
        df = pd.DataFrame(columns=['source', 'target'])
        return df
    # Join
    df = pd.merge(tfb[['tf', 'cre']], p2g[['cre', 'gene']], how='inner', on='cre')[['tf', 'cre', 'gene']]
    df = df.sort_values(['tf', 'cre', 'gene']).rename(columns={'tf': 'source', 'gene': 'target'})
    df['score'] = 1.

    tfs = df['source'].unique().astype('U')
    genes = df['target'].unique().astype('U')
    n = int(np.round(genes.size * tf_g_ratio))

    rng = np.random.default_rng(seed=seed)
    tfs = rng.choice(tfs, n, replace=False)
    df = df.loc[df['source'].astype('U').isin(tfs), :]

    return df.reset_index(drop=True)


# Run random
mdata = mu.read(inp_path)
mdata = run_pre(mdata, seed)
gannot = pr.read_bed(cg_path)
p2g = run_p2g(mdata, gannot, g_perc, scale, w_size, seed)
tfs = pd.read_csv(tf_path, header=None).loc[:, 0].values.astype('U')
tfb = run_tfb(mdata, p2g, tfs, scale, seed)
grn = run_mdl(p2g, tfb, tf_g_ratio, seed)

# Write
grn.to_csv(out_path, index=False)
