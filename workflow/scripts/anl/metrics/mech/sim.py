#!/usr/bin/env python
# coding: utf-8

# In[4]:


# Load libraries
import os
import argparse
import pandas as pd
import numpy as np
import pyboolnet
from pyboolnet import file_exchange, trap_spaces
from statsmodels.stats.multitest import multipletests
from scipy.stats import fisher_exact
from sklearn.preprocessing import StandardScaler

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-g','--grn_path', required=True)
parser.add_argument('-d','--data_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())
grn_path = args['grn_path']
data_path = args['data_path']
out_path = args['out_path']

# Load data
data = pd.read_csv(data_path, index_col = 0)
grns = pd.read_csv(grn_path, index_col = 0)

# Load TF marker data    
# This is what can be changed


# Get the top 15 regulons x cell type based on mat
scaler = StandardScaler()
scaled_mat = pd.DataFrame(scaler.fit_transform(data), index=data.index, columns=data.columns)
ranked_mat = scaled_mat.rank(axis=1, ascending=False, method="min")

# Extract top 15 rankings and convert to dictionary
top_15_rankings_df = ranked_mat.apply(lambda row: row.nsmallest(15).index, axis=1).apply(pd.Series)
top_15_rankings_df.columns = [f"Rank_{i+1}" for i in range(top_15_rankings_df.shape[1])]
gold_standard_dict = top_15_rankings_df.apply(lambda row: row.dropna().tolist(), axis=1).to_dict()
all_tfs_universe = set(tf for tfs in gold_standard_dict.values() for tf in tfs)

# Calculate steady states
def analyze_df(grn, mat):
    # Generate Boolean rules
    rules = {}
    for target in grn['target'].unique():
        target_df = grn[grn['target'] == target]
        activators = target_df[target_df['score'] == 1]['source'].tolist()
        inhibitors = target_df[target_df['score'] == -1]['source'].tolist()
        activator_rule = " | ".join(activators) if activators else ""   # We are using OR not AND.
        inhibitor_rule = " | ".join([f"!{inh}" for inh in inhibitors]) if inhibitors else ""
        if activator_rule and inhibitor_rule:
            rule = f"{activator_rule} & {inhibitor_rule}"
        elif activator_rule:
            rule = activator_rule
        elif inhibitor_rule:
            rule = inhibitor_rule
        rules[target] = rule

    bnet = "\n".join([f"{target}, {rule}" for target, rule in rules.items()])
    primes = file_exchange.bnet2primes(bnet)

    # Compute steady states
    steady_states = trap_spaces.compute_steady_states(primes, max_output=1000000)  #Change max_output?
    df_steady = pd.DataFrame(steady_states)

    return df_steady
df_steady = analyze_df(grns, data)

def compute_adjusted_pvalues(df_steady, gold_standard_dict, all_tfs_universe):
    all_results = []

    for row_idx, row in df_steady.iterrows():
        selected_tfs = set(row[row == 1].index.tolist())
        row_results = []

        for cell_type, gold_standard_tfs in gold_standard_dict.items():
            gold_standard_set = set(gold_standard_tfs)

            selected_in_gold = len(selected_tfs & gold_standard_set)
            selected_not_in_gold = len(selected_tfs) - selected_in_gold
            non_selected_in_gold = len(gold_standard_set) - selected_in_gold
            non_selected_not_in_gold = len(all_tfs_universe - (selected_tfs.union(gold_standard_set)))

            contingency_table = [
                [selected_in_gold, selected_not_in_gold],
                [non_selected_in_gold, non_selected_not_in_gold]
            ]

            odds_ratio, p_value = fisher_exact(contingency_table, alternative="greater")

            row_results.append({
                "Row": row_idx,
                "Cell Type": cell_type,
                "Odds Ratio": odds_ratio,
                "P-value": p_value
            })

        row_df = pd.DataFrame(row_results)
        row_df["Adjusted P-value"] = multipletests(row_df["P-value"], method="fdr_bh")[1]
        all_results.append(row_df)

    final_results_df = pd.concat(all_results, ignore_index=True)
    return final_results_df.pivot(index="Row", columns="Cell Type", values="Adjusted P-value")

df_final = compute_adjusted_pvalues(df_steady, gold_standard_dict, all_tfs_universe)

df_final.to_csv(out_path)



# In[ ]:




