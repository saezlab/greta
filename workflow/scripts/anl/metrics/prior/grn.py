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
grn_name_simple = grn_name.split('.')[0].replace('o_', '')
if grn.shape[0] > 0 and grn_name_simple not in ['collectri', 'dorothea']:
    # Read resource and filter by cats
    db = pd.read_csv(resource_path)
    
    # Filter resource by measured genes
    genes = mu.read(os.path.join(data_path, 'mod', 'rna')).var_names.astype('U')
    db = db[db['source'].astype('U').isin(genes) & db['target'].astype('U').isin(genes)]
    
    # Compute evaluation
    set_db = set(db['source'] + '|' + db['target'])
    set_grn = set(grn['source'] + '|' + grn['target'])
    tp = len(set_db & set_grn)
    if tp > 0.:
        fp = len(set_grn - set_db)
        fn = len(set_db - set_grn)
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
