import pandas as pd
import numpy as np
import yaml
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(os.path.curdir), 'workflow', 'scripts', 'analysis')))
from utils import make_combs
from tqdm import tqdm
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i','--sims_path', required=True)
args = vars(parser.parse_args())

sims_path = args['sims_path']

# Read
dataset = os.path.basename(sims_path).split('.')[0].replace('pair', '')
case = os.path.basename(sims_path).split('.')[1]

def compute_stats(ref_path, prd_path):
    on = ['source', 'target']
    ref = pd.read_csv(ref_path).drop_duplicates(on)
    prd = pd.read_csv(prd_path).drop_duplicates(on)

    if (ref.shape[0] > 0) and (prd.shape[0] > 0):

        ref = set([s + '|' + t for s, t in zip(ref['source'], ref['target'])])
        prd = set([s + '|' + t for s, t in zip(prd['source'], prd['target'])])
    
        tp = len(ref & prd)
        if tp > 0:
            fp = len(prd - ref)
            fn = len(ref - prd)
            
            prc = tp / (tp + fp)
            rcl = tp / (tp + fn)
            
            beta = 0.1
            f01 = ((1 - beta**2) * prc * rcl) / ((prc * beta**2) + rcl)
        else:
            prc, rcl, f01 = 0., 0., 0.
    else:
        prc, rcl, f01 = np.nan, np.nan, np.nan
    return prc, rcl, f01


with open('config/config.yaml') as stream:
    config = yaml.safe_load(stream)
mthds = list(config['methods'].keys())


combs_scores = make_combs(
    path='analysis/metrics/pred/pair/{d}/{d}.{c}/'.format(d=dataset, c=case),
    mthds=mthds,
    name='scores',
)

combs_pair = make_combs(
    path='datasets/{d}pair/cases/{c}/runs/'.format(d=dataset, c=case),
    mthds=mthds,
    name='grn',
)

combs_npair = make_combs(
    path='datasets/{d}npair/cases/{c}/runs/'.format(d=dataset, c=case),
    mthds=mthds,
    name='grn',
)


with tqdm(total=len(combs_pair)) as pbar:
    for p_path, n_path, s_path in zip(combs_pair, combs_npair, combs_scores):
        grn_name_p = os.path.basename(p_path).replace('.grn.csv', '')
        grn_name_n = os.path.basename(n_path).replace('.grn.csv', '')
        assert grn_name_p == grn_name_n
        prc, rcl, f01 = compute_stats(p_path, n_path)
        df = pd.DataFrame([[grn_name_p, prc, rcl, f01]], columns=['name', 'prc', 'rcl', 'f01'])
        df.to_csv(s_path, index=False)
        pbar.update(1)
