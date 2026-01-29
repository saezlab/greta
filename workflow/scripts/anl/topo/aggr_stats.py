import pandas as pd
import scipy.stats as sts
import os
import argparse


def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--paths_stats', required=True, nargs='+')
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

paths_stats = args['paths_stats']
path_out = args['path_out']

config = read_config()
mth_dict = config['method_names']
dts_dict = config['dts_names']

df = []
for path_stats in paths_stats:
    dts_name = os.path.basename(path_stats).split('.')[1]
    if dts_name in dts_dict:
        dts_name = dts_dict[dts_name]
        stat = pd.read_csv(path_stats)
        stat = stat[stat['name'].str.startswith('o_')]
        stat['name'] = [mth_dict[n.split('.')[0].replace('o_', '')] for n in stat['name']]
        stat = stat.drop(columns=['betweenc', 'eigv'])
        stat['dts'] = dts_name
        df.append(stat)
df = pd.concat(df)

sc_dict = {
    'n_tfs': '# TFs',
    'n_cres': '# CREs',
    'n_targets': '# Genes',
    'n_edges': '# Edges',
    'odegree': 'Regulon size'
}
df = df.rename(columns=sc_dict)

# Write
df.to_csv(path_out, index=False)
