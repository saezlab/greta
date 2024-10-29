import pandas as pd
import numpy as np
import os


def get_grn_name(grn_path):
    name = os.path.basename(grn_path).replace('.grn.csv', '').replace('.csv', '')
    pre, p2g, tfb, mdl = name.split('.')
    if (pre == p2g) & (p2g == tfb) & (tfb == mdl):
        name = pre
    return name


def get_grn_stats(grn):
    import igraph as ig
    n_s = grn['source'].unique().size
    n_e = grn.shape[0]
    n_t = grn['target'].unique().size
    n_r = grn.groupby(['source']).count()['target'].mean()
    
    tfs = set(grn['source']) & set(grn['target'])
    msk = grn['source'].isin(tfs) & grn['target'].isin(tfs)
    tf_grn = grn.loc[msk, :]
    tf_g = ig.Graph.TupleList(list(zip(tf_grn['source'], tf_grn['target'])), directed=True)
    tf_bet = np.mean(tf_g.betweenness())
    tf_odg = np.mean(tf_g.outdegree())
    if not tf_g.is_acyclic():
        tf_eig = np.mean(tf_g.eigenvector_centrality())
    else:
        tf_eig = 0.
    
    return n_s, n_e, n_t, n_r, tf_odg, tf_bet, tf_eig


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


def parallel_ocoeff(index_pair, dfs):
    i, j = index_pair
    tf_oc = ocoeff(dfs[i], dfs[j], on=['source'])
    edge_oc = ocoeff(dfs[i], dfs[j], on=['source', 'target'])
    target_oc = ocoeff(dfs[i], dfs[j], on=['target'])
    return i, j, tf_oc, edge_oc, target_oc


def parallel_ocoeff_chunk(index_pairs_chunk, dfs):
    return [parallel_ocoeff(pair, dfs) for pair in index_pairs_chunk]


def make_combs(path, mthds, name):
    from itertools import product
    s = '{0}.{1}.{2}.{3}.' + name + '.csv'
    combinations = product(mthds, repeat=4)
    strings = []
    for combo in combinations:
        formatted_string = path + s.format(*combo)
        strings.append(formatted_string)

    # Add src
    s = '{m}.{m}.{m}.{m}.' + name + '.csv'
    for m in mthds:
        formatted_string = path + s.format(m='o_' + m)
        strings.append(formatted_string)

    # Add indv nets
    s = path + 'random.random.random.random.' + name + '.csv'
    strings.append(s)
    s = path + 'scenic.scenic.scenic.scenic.' + name + '.csv'
    strings.append(s)
    s = path + 'collectri.collectri.collectri.collectri.' + name + '.csv'
    strings.append(s)
    s = path + 'dorothea.dorothea.dorothea.dorothea.' + name + '.csv'
    strings.append(s)

    return strings
