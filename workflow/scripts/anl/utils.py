import pandas as pd
import numpy as np
import os


def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


def get_grn_name(grn_path):
    name = os.path.basename(grn_path).replace('.grn.csv', '').replace('.csv', '')
    return name


def get_grn_stats(grn):
    import igraph as ig
    if len(grn) == 0:
        return np.nan, np.nan, np.nan, np.nan, np.nan, np.nan
    n_s = grn['source'].unique().size
    n_e = grn.shape[0]
    n_t = grn['target'].unique().size

    g = ig.Graph.TupleList(list(zip(grn['source'], grn['target'])), directed=True)
    tf_bet = np.mean(g.betweenness())
    tf_odg = grn.groupby(['source']).size().mean()
    if not g.is_acyclic():
        tf_eig = np.mean(g.eigenvector_centrality())
    else:
        tf_eig = 0.
    
    return n_s, n_e, n_t, tf_odg, tf_bet, tf_eig


def ocoeff(df_a, df_b, on=['source', 'target']):
    """Compute overlap coefficient between two dfs"""
    tmp_a, tmp_b = df_a.drop_duplicates(on), df_b.drop_duplicates(on)
    a_size, b_size = tmp_a.shape[0], tmp_b.shape[0]
    if (a_size > 0) and (b_size > 0):
        inter = pd.merge(tmp_a, tmp_b, on=on, how='inner')
        i_size = inter.shape[0]
        coeff = i_size / np.min([a_size, b_size])
    else:
        coeff = 0.
    return coeff
