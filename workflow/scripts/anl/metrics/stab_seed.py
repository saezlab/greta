import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
from glob import glob
import sys

path_out = sys.argv[1]

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config

config = read_config()
mth_dict = config['method_names']
class_dict = config['class_names']
task_names = config['task_names']
db_names = config['dbs_names']

paths = glob('anl/metrics/*/*/*/hg38.pitupair.16384_16384_*/*.scores.csv')
df = []
for path in paths:
    tmp = pd.read_csv(path)[['prc', 'rcl', 'f01']]
    path_lst = path.split('/')
    tmp['class'] = class_dict[path_lst[2]]
    tmp['task'] = task_names[path_lst[3]]
    tmp['db'] = db_names[path_lst[4]]
    _, _, tmp['case'] = path_lst[5].split('.')
    tmp['name'] = mth_dict[path_lst[6].split('.')[0].replace('o_', '')]
    row = ['class', 'task', 'db', 'case', 'name', 'prc', 'rcl', 'f01']
    tmp = tmp[row]
    df.append(tmp)
df = pd.concat(df)
df = df.fillna(0)
df.to_csv(path_out, index=False)
