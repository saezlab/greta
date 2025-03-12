import pandas as pd
import numpy as np
import mudata as mu
from tqdm import tqdm
import sys
import os
import re
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import load_cats, f_beta_score
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--grn_path', required=True)
parser.add_argument('-b','--resource_path', required=True)
parser.add_argument('-f','--out_path', required=True)
args = vars(parser.parse_args())

grn_path = args['grn_path']
resource_path = args['resource_path']
out_path = args['out_path']


grn_name = os.path.basename(grn_path).replace('.grn.csv', '')
data_path = os.path.join(os.path.dirname(os.path.dirname(grn_path)), 'mdata.h5mu')
dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(data_path))))
case = os.path.basename(os.path.dirname(data_path))
resource_name = os.path.basename(resource_path).replace('.csv', '')

# Read grn
grn = pd.read_csv(grn_path)

if grn.shape[0] > 0:
    # Read resource and filter by cats
    db = pd.read_csv(resource_path, header=None, sep='\t')
    db.columns = ['gene', 'ctype']
    cats = load_cats(dataset, case)
    if resource_name in cats:
        cats = [re.escape(c) for c in cats[resource_name]]
        print('Filtering for {0} cats'.format(len(cats)))
        db = db[db['ctype'].str.contains('|'.join(cats))]
    
    # Filter resource by measured genes
    genes = mu.read(os.path.join(data_path, 'mod', 'rna')).var_names.astype('U')
    db = db[db['gene'].astype('U').isin(genes)]
    
    # Compute evaluation
    y_pred = grn['source'].unique().astype('U')
    y = db['gene'].unique().astype('U')
    tp = np.intersect1d(y_pred, y).size
    if tp > 0.:
        fp = np.setdiff1d(y_pred, y).size
        fn = np.setdiff1d(y, y_pred).size
        prc = tp / (tp + fp)
        rcl = tp / (tp + fn)
        f01 = f_beta_score(prc, rcl)
    else:
        prc, rcl, f01 = 0., 0., 0.,
    df = pd.DataFrame([[grn_name, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
else:
    df = pd.DataFrame([[grn_name, np.nan, np.nan, np.nan]], columns=['name', 'prc', 'rcl', 'f01'])

# Write
df.to_csv(out_path, index=False)
