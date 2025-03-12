import pandas as pd
import numpy as np
import pyranges as pr
import mudata as mu
from tqdm import tqdm
import re
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import load_cats, f_beta_score
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--grn_path', required=True)
parser.add_argument('-b','--resource_path', required=True)
parser.add_argument('-d','--grp', default='None')
parser.add_argument('-f','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
resource_path = args['resource_path']
grp = args['grp']
out_path = args['out_path']

if grp == 'None':
    grp = None

grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')
dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(data_path))))
case = os.path.basename(os.path.dirname(data_path))
resource_name = os.path.basename(resource_path).replace('.bed', '')

def read_grn(grn_path, grp=None):
    grn = pd.read_csv(grn_path)
    if 'cre' in grn.columns:
        if grn.shape[0] > 0:
            grn[['Chromosome', 'Start', 'End']] = grn['cre'].str.split('-', n=2, expand=True)
            if grp is not None:
                grn = grn[['Chromosome', 'Start', 'End', grp]].rename(columns={grp: 'Name'})
            else:
                grn = grn[['Chromosome', 'Start', 'End']]
        else:
            grn = None
    else:
        grn = None
    return pr.PyRanges(grn)

grn = read_grn(grn_path=grn_path, grp=grp)

if grn.df.shape[0] > 0:
    # Read resource and filter by cats
    db = pr.read_bed(resource_path)
    cats = load_cats(dataset, case)
    if resource_name in cats:
        cats = [re.escape(c) for c in cats[resource_name]]
        print('Filtering for {0} cats'.format(len(cats)))
        db = db[db.df['Score'].str.contains('|'.join(cats))]
    
    # Filter genomic resource by measured CREs and or genes
    genes = mu.read(os.path.join(data_path, 'mod', 'rna')).var_names.astype('U')
    peaks = mu.read(os.path.join(data_path, 'mod', 'atac')).var_names.astype('U')
    peaks = pd.DataFrame(peaks, columns=['cre'])
    peaks[['Chromosome', 'Start', 'End']] = peaks['cre'].str.split('-', n=2, expand=True)
    peaks = pr.PyRanges(peaks[['Chromosome', 'Start', 'End']])
    db = db.overlap(peaks)
    
    if grp is not None:
        # Remove features that are in GRN but not in db
        db = db[db.df.Name.astype('U').isin(genes)]
        grn_feats = grn.df['Name'].unique()
        db_feats = db.df['Name'].unique()
        features = np.setdiff1d(grn_feats, db_feats)
        grn = grn[~grn.df['Name'].isin(features)]
    
        tps = 0
        fps = 0
        fns = 0
        if grn.df.shape[0] > 0:
            for feature in db_feats:
                f_grn = grn[grn.df['Name'] == feature]
                f_db = db[db.df['Name'] == feature]
                tps += f_grn.overlap(f_db).df.shape[0]
                fps += f_grn.overlap(f_db, invert=True).df.shape[0]
                fns += f_db.overlap(f_grn, invert=True).df.shape[0]
    elif resource_name != 'blacklist':
        tps = grn.overlap(db).df.shape[0]
        fps = grn.overlap(db, invert=True).df.shape[0]
        fns = db.overlap(grn, invert=True).df.shape[0]
    else:
        tps = grn.overlap(db, invert=True).df.shape[0]
        fps = grn.overlap(db).df.shape[0]
        fns = peaks.overlap(db, invert=True).df.shape[0]
    if tps > 0:
        prc = tps / (tps + fps)
        rcl = tps / (tps + fns)
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
