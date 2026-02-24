import pandas as pd
import sys

path_inp = sys.argv[1]  # 'anl/stab/pitupair.ovc.csv'
path_dsd = sys.argv[2]
path_out = sys.argv[3]

uses_gpu = ['dictys', 'scdori', 'scgpt']

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


config = read_config()
mth_dict = config['method_names']

df = pd.read_csv(path_inp)
df = df[df['cat'] == 'full']
df = df.groupby(['mth'], as_index=False)[['e_ocoeff', 'h', 'gb']].mean()
df['use_gpu'] = df['mth'].isin(uses_gpu)
df['mth'] = [mth_dict[n] for n in df['mth']]
df = df.rename(columns={'e_ocoeff': 'stability'})
dsd = pd.read_csv(path_dsd).groupby('name')['oc_edge'].mean()
df['stability'] = [dsd[m] for m in df['mth']]
df.to_csv(path_out, index=False)
