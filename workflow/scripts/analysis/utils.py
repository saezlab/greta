import pandas as pd
import numpy as np


def ocoeff(df_a, df_b, on=['source', 'target']):
    """Compute overlap coefficient between two dfs"""
    a_size, b_size = df_a.shape[0], df_b.shape[0]
    if (a_size > 0) and (b_size > 0):
        inter = pd.merge(df_a, df_b, on=on, how='inner')
        i_size = inter.shape[0]
        coeff = i_size / np.min([a_size, b_size])
    else:
        coeff = 0.
    return coeff


def parallel_ocoeff(index_pair):
    i, j = index_pair
    tf_oc = ocoeff(dfs[i], dfs[j], on=['source'])
    edge_oc = ocoeff(dfs[i], dfs[j], on=['source', 'target'])
    target_oc = ocoeff(dfs[i], dfs[j], on=['target'])
    return i, j, tf_oc, edge_oc, target_oc
