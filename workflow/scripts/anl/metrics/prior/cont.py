#!/usr/bin/env python
# coding: utf-8

# In[ ]:


from itertools import combinations
from tqdm import tqdm
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import pandas as pd
from joblib import Parallel, delayed
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tf_path', required=True)
parser.add_argument('-m','--grn_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

tf_path = args['tf_path']
grn_path = args['grn_path']
out_path = args['out_path']

# Read data
tfs = pd.read_csv(tf_path, header=None)
grn = pd.read_csv(grn_path)
grn = grn[grn["source"].isin(tfs[0])]

def process_pair(tf1, tf2, method):
    targets_tf1 = set(method[method['source'] == tf1]['target'])
    targets_tf2 = set(method[method['source'] == tf2]['target'])
    a = len(targets_tf1 & targets_tf2)  # Intersection of targets 
    # Skip pairs with no intersection
    if a == 0:
        return None
    # Calculate remaining contingency table components
    b = len(targets_tf1 - targets_tf2)  # Targets unique to tf1
    c = len(targets_tf2 - targets_tf1)  # Targets unique to tf2
    d = len(set(method['target']) - (targets_tf1 | targets_tf2))  # All other target
    # Create contingency table
    contingency_table = [[a, b], [c, d]]
    # Perform Fisher's exact test
    _, p_value = fisher_exact(contingency_table, alternative='greater')
    return {
        'TF1': tf1,
        'TF2': tf2,
        'a (Intersection)': a,
        'b (TF1 only)': b,
        'c (TF2 only)': c,
        'd (Others)': d,
        'p-value': p_value
    }

def calculate_contingency(method):
    unique_tfs = method['source'].unique()
    pairs = list(combinations(unique_tfs, 2))
    contingency_results = Parallel(n_jobs=-1)(
        delayed(process_pair)(tf1, tf2, method) for tf1, tf2 in tqdm(pairs, desc="Processing pairs", unit="pair")
    )
    contingency_results = [result for result in contingency_results if result is not None]
    contingency_df = pd.DataFrame(contingency_results)
    if not contingency_df.empty:
        p_values = contingency_df['p-value']
        _, fdr_corrected_p_values, _, _ = multipletests(p_values, method='fdr_bh')
        contingency_df['p-value (FDR)'] = fdr_corrected_p_values
        contingency_df = contingency_df[contingency_df['p-value (FDR)'] < 2.2e-16]
    return contingency_df

filt_cont = calculate_contingency(grn)

filt_cont.to_csv(out_path)

