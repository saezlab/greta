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
parser.add_argument('-a','--path_a', required=True)
parser.add_argument('-b','--path_b', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_a = args['path_a']
path_b = args['path_b']
path_out = args['path_out']

# Find paths
dname_a, case_a = os.path.basename(path_a).split('.')[:2]
dname_b, case_b = os.path.basename(path_b).split('.')[:2]
dname_a = dname_a.replace('pair', '')
dname_b = dname_b.replace('pair', '')
path_pair = sorted(glob.glob(f'dts/{dname_a}pair/cases/{case_b}/runs/*.grn.csv'))
path_npair = sorted(glob.glob(f'dts/{dname_b}pair/cases/{case_b}/runs/*.grn.csv'))

# Compute ocoef
df = []
for i in tqdm(range(len(path_pair))):
    p_path, n_path = path_pair[i], path_npair[i]
    assert os.path.basename(p_path) == os.path.basename(n_path)
    p_grn, n_grn = pd.read_csv(p_path), pd.read_csv(n_path)
    val = ocoeff(p_grn, n_grn, on=['source', 'target'])
    df.append([os.path.basename(p_path).replace('.grn.csv', ''), val])
df = pd.DataFrame(df, columns=['mth', 'ocoef'])

# Write
df.to_csv(path_out, index=False)
