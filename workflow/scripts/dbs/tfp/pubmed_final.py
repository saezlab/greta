#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
from scipy.stats import fisher_exact
from statsmodels.stats.multitest import multipletests
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-p','--pairs_path', required=True)
parser.add_argument('-t','--tf_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

pairs_path = args['pairs_path']
tf_path = args['tf_path']
out_path = args['out_path']

# Read data
tf_pairs = pd.read_csv(pairs_path, sep="\t")
tf_counts = pd.read_csv(tf_path, sep="\t")

# Filter tf pairs that are 0s and preprocess
tf_pairs['Abstracts'] = pd.to_numeric(tf_pairs['Abstracts'])
filt_tf_pairs = tf_pairs[tf_pairs['Abstracts'] > 0].copy()
filt_tf_pairs[['TF1', 'TF2']] = filt_tf_pairs['TF_Pair'].str.split(',', expand=True)
filt_tf_pairs.drop(columns=['TF_Pair'], inplace=True)

def calculate_fisher_and_filter(tf_counts, filt_tf_pairs, fdr_threshold=0.05, p_value_cap=2.2e-16):
    # Prepare data
    tf_counts_dict = dict(zip(tf_counts['TF'], tf_counts['Abstracts']))
    total_abstracts = tf_counts['Abstracts'].sum()

    # Calculate Fisher test for each pair
    results = []
    for _, row in filt_tf_pairs.iterrows():
        tf1, tf2, pair_count = row['TF1'], row['TF2'], row['Abstracts']
        tf1_count, tf2_count = tf_counts_dict.get(tf1, 0), tf_counts_dict.get(tf2, 0)
        a, b, c, d = pair_count, tf1_count - pair_count, tf2_count - pair_count, total_abstracts - (pair_count + tf1_count + tf2_count - pair_count)
        if a + b + c > total_abstracts or d < 0:
            continue
        try:
            odds_ratio, p_value = fisher_exact([[a, b], [c, d]], alternative='greater')
        except ValueError:
            continue
        results.append({'TF1': tf1, 'TF2': tf2, 'P-Value': p_value, 'Odds Ratio': odds_ratio})
    
    results_df = pd.DataFrame(results)
    results_df['FDR_Adjusted_P'] = multipletests(results_df['P-Value'], method='fdr_bh')[1]
    final_pairs = results_df[results_df['FDR_Adjusted_P'] < p_value_cap][['TF1', 'TF2']]
    return final_pairs

final_results = calculate_fisher_and_filter(tf_counts, filt_tf_pairs)

final_results.to_csv(out_path)



