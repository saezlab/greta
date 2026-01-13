import pandas as pd
import pyranges as pr
from tqdm import tqdm
import os
import glob
import argparse


# Parse args
parser = argparse.ArgumentParser()
parser.add_argument('-g', '--path_cmp', required=True)
parser.add_argument('-b','--baselines', required=True, nargs='+')
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()
path_cmp = args.path_cmp
baselines = args.baselines
path_out = args.path_out

# Set variables
org, dname, case = os.path.basename(path_cmp).split('.')[:3]
path_glob = os.path.join('dts', org, dname, 'cases', case, 'runs', '*.grn.csv')
print(path_glob)
path_grns = glob.glob(path_glob)
def compute_dist_tss(path, mth, is_baseline):
    grn = pd.read_csv(path)
    grn = grn.drop_duplicates(['cre', 'target'])
    grn[['Chromosome', 'Start', 'End']] = grn['cre'].str.split('-', expand=True)
    grn = pr.PyRanges(grn[['Chromosome', 'Start', 'End', 'target']].rename(columns={'target': 'Name'}))
    if is_baseline:
        path_tss = 'dbs/hg38/cre/promoters/promoters.bed.gz'
    else:
        path_tss = f'dbs/hg38/gen/tss/{mth}.bed.gz'
    tss = pd.read_csv(path_tss, sep='\t', header=None)
    tss.columns = ['Chromosome', 'Start', 'End', 'Name']
    tss = pr.PyRanges(tss)
    genes = grn.df['Name'].unique().astype('U')
    dists = []
    for g in genes:
        g_grn = grn[grn.Name == g]
        g_tss = tss[tss.Name == g]
        if g_grn.df.shape[0] > 0 and g_tss.df.shape[0] > 0:
            dists.append(g_grn.nearest(g_tss, overlap=True).df[['Chromosome', 'Start', 'End', 'Distance']].assign(gene=g))
    dists = pd.concat(dists).rename(columns={'Distance': 'dist'})
    dists['mth'] = mth
    dists['cre'] = dists['Chromosome'].astype(str) + '-' + dists['Start'].astype(str) + '-' + dists['End'].astype(str)
    dists = dists[['mth', 'cre', 'gene', 'dist']]
    return dists

# Compute dists
dists = []
path_grns = [p for p in path_grns if os.path.basename(p).startswith('o_')]
print(path_grns)
for path_grn in tqdm(path_grns):
    mth = os.path.basename(path_grn).split('.')[0]  # Assume all stp equal
    mth = mth.replace('o_', '')
    is_baseline = mth in baselines
    print(mth, is_baseline)
    dists.append(compute_dist_tss(path_grn, mth, is_baseline))
dists = pd.concat(dists)

# Write
dists.to_csv(path_out, index=False)
