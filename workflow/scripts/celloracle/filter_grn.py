import pandas as pd
import numpy as np
import celloracle as co
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-l','--path_links', required=True)
parser.add_argument('-p','--thr_edge_pval', required=True)
parser.add_argument('-t','--thr_top_edges', required=True)
parser.add_argument('-o','--path_grns', required=True)
args = vars(parser.parse_args())

path_links = args['path_links']
thr_edge_pval = float(args['thr_edge_pval'])
thr_top_edges = int(args['thr_top_edges'])
path_grns = args['path_grns']

# Load links
links = co.load_hdf5(file_path=path_links)

# Filter
links.filter_links(p=thr_edge_pval, weight="coef_abs", threshold_number=thr_top_edges)

# Write GRN per cluster
for celltype in links.links_dict.keys():
    links.links_dict[cluster].to_csv(os.path.join(path_grns, "{0}.csv".format(celltype))

