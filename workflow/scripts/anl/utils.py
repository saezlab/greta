import pandas as pd
import numpy as np
import pyranges as pr
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
    n_c = grn['cre'].unique().size
    n_t = grn['target'].unique().size
    tgrn = grn.drop_duplicates(['source', 'target'])
    n_e = tgrn.shape[0]

    g = ig.Graph.TupleList(list(zip(tgrn['source'], tgrn['target'])), directed=True)
    tf_bet = np.mean(g.betweenness())
    tf_odg = tgrn.groupby(['source']).size().mean()
    if not g.is_acyclic():
        tf_eig = np.mean(g.eigenvector_centrality())
    else:
        tf_eig = 0.
    
    return n_s, n_c, n_t, n_e, tf_odg, tf_bet, tf_eig



def _map_regions(regions_a: np.ndarray, regions_b: np.ndarray) -> pd.DataFrame:
    """Map overlapping genomic regions between two sets of peaks."""

    def to_df(regions):
        data = [r.split("-") for r in regions]
        return pd.DataFrame(data, columns=["Chromosome", "Start", "End"]).assign(
            Start=lambda x: x.Start.astype(int),
            End=lambda x: x.End.astype(int),
            region=regions,
        )

    pr_a = pr.PyRanges(to_df(regions_a))
    pr_b = pr.PyRanges(to_df(regions_b))
    joined = pr_b.join(pr_a, suffix="_a")
    if joined.empty:
        return pd.DataFrame(columns=["region_a", "region_b"])
    return joined.df[["region", "region_a"]].rename(columns={"region": "region_b"})


def ocoeff(df_a: pd.DataFrame, df_b: pd.DataFrame, on: list, use_overlap: bool = False) -> float:
    tmp_a, tmp_b = df_a.drop_duplicates(on), df_b.drop_duplicates(on)
    a_size, b_size = tmp_a.shape[0], tmp_b.shape[0]
    if (a_size > 0) and (b_size > 0):
        if use_overlap and len(on) == 1 and on[0] == "cre":
            regions_a = tmp_a["cre"].unique()
            regions_b = tmp_b["cre"].unique()
            mapping = _map_regions(regions_a, regions_b)
            if mapping.empty:
                i_size = 0
            else:
                # Count overlapping regions from the smaller set
                if a_size <= b_size:
                    i_size = mapping["region_a"].nunique()
                else:
                    i_size = mapping["region_b"].nunique()
        else:
            inter = pd.merge(tmp_a, tmp_b, on=on, how="inner")
            i_size = inter.shape[0]
        coeff = i_size / np.min([a_size, b_size])
    else:
        coeff = 0.0
    return coeff
