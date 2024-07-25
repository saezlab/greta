import os
import argparse
os.environ[ 'NUMBA_CACHE_DIR' ] = '/tmp/'

import pandas as pd
import muon as mu

import circe as ci

from distributed import LocalCluster, Client
from arboreto.algo import grnboost2

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--path_mudata', required=True)
parser.add_argument('-r', '--path_grnboost2', required=True)
parser.add_argument('-a', '--path_atacnet', required=True)
parser.add_argument('-t', '--distance_threshold', required=True, type=int)
parser.add_argument('-c', '--n_cores', required=True, type=int)
parser.add_argument('-n', '--tf_names', required=True)
args = vars(parser.parse_args())


def run_grnboost2(expression_data, tf_names, n_cores=1):

    local_cluster = LocalCluster(n_workers=n_cores, threads_per_worker=1)
    custom_client = Client(local_cluster)

    rna_network = grnboost2(
        expression_data=expression_data,
        tf_names=tf_names,
        client_or_address=custom_client)
    custom_client.close()

    return rna_network


path_mudata = args['path_mudata']
path_grnboost2 = args['path_grnboost2']
path_atacnet = args['path_atacnet']
distance_threshold = args['distance_threshold']
n_cores = args['n_cores']
tf_names = args['tf_names']

# Load the data
mudata = mu.read_h5mu(path_mudata)

if __name__ == '__main__':
    print('Computing networks...')
    tf_names = pd.read_csv(tf_names, header=None)[0].tolist()

    rna_network = run_grnboost2(
        expression_data=mudata["rna"].to_df(),
        tf_names=tf_names,
        n_cores=n_cores)
    rna_network.to_csv(path_grnboost2)

    # Create the atacnet network
    atac = ci.add_region_infos(mudata["atac"], sep=(':', '-'))
    mudata = mu.MuData(
        {"rna": mudata["rna"],
         "atac": atac})
    ci.compute_atac_network(
        mudata["atac"],
        window_size=distance_threshold)

    atac_network = ci.extract_atac_links(mudata["atac"])
    atac_network.to_csv(path_atacnet)



