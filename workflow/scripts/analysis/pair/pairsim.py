import pandas as pd
import numpy as np
import yaml
import sys
import os
from tqdm import tqdm
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import (
    ocoeff,
)
import glob
import argparse


parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_inp', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_inp = args['path_inp']
path_out = args['path_out']

# Find paths
dname, case = os.path.basename(path_inp).split('.')[:2]
dname = dname.replace('pair', '')
path_pair = sorted(glob.glob(f'datasets/{dname}pair/cases/{case}/runs/*.grn.csv'))
path_npair = sorted(glob.glob(f'datasets/{dname}npair/cases/{case}/runs/*.grn.csv'))

# Compute ocoef
df = []
for i in tqdm(range(path_pair)):
    p_path, n_path = path_pair[i], path_npair[i]
    assert os.path.basename(p_path) == os.path.basename(n_path)
    p_grn, n_grn = pd.read_csv(p_path), pd.read_csv(n_path)
    val = ocoeff(p_grn, n_grn, on=['source', 'target'])
    df.append([p_m_name, val])
df = pd.DataFrame(df, columns=['mth', 'ocoef'])
