#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from itertools import combinations
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
from tqdm import tqdm  
from joblib import Parallel, delayed
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--db_path', required=True)
parser.add_argument('-c','--cont_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

db_path = args['db_path']
cont_path = args['cont_path']
out_path = args['out_path']

db = pd.read_csv(db_path)
filt_cont = pd.read_csv(cont_path)

def calculate_f01(results_df, database, beta):
    if results_df.empty:
        # Return a DataFrame with 0s if the contingency matrix is empty
        return pd.DataFrame([{
            "TP": 0,
            "FP": 0,
            "FN": len(database),  # All database pairs are considered false negatives
            "Precision": 0,
            "Recall": 0,
            "F0.1 Score": 0
        }])
        
    significant_pairs = set(map(tuple, results_df[['TF1', 'TF2']].apply(sorted, axis=1)))
    database_pairs = set(map(tuple, database[['TF1', 'TF2']].apply(sorted, axis=1)))

    # Compute metrics
    TP = len(significant_pairs & database_pairs)
    FP = len(significant_pairs - database_pairs)
    FN = len(database_pairs - significant_pairs)

    precision = TP / (TP + FP) if (TP + FP) > 0 else 0
    recall = TP / (TP + FN) if (TP + FN) > 0 else 0
    f01_score = (1 + beta**2) * precision * recall / ((beta**2 * precision) + recall) if (precision + recall) > 0 else 0

    return pd.DataFrame([{
        "TP": TP,
        "FP": FP,
        "FN": FN,
        "Precision": precision,
        "Recall": recall,
        "F0.1 Score": f01_score
    }])

f01_results = calculate_f01(filt_cont, db, beta=0.1)
f01_results.to_csv(out_path)


