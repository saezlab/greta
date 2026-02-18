import pandas as pd
import glob
import sys
import yaml

path_out = sys.argv[1]


def read_config(path_config='config/config.yaml'):
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config


config = read_config()
class_dict = config['class_names']
db_dict = config['dbs_names']
task_dict = config['task_names']
dts_dict = config['dts_names']
mth_dict = config['method_names']

lst_paths = glob.glob('anl/metrics/*/*/*/*.*.all.scores.csv')
df = []
for path in lst_paths:
    spath = path.split('/')
    class_metric = spath[2]
    task_metric = spath[3]
    name_db = spath[4]
    name_org = spath[-1].split('.')[0]
    name_dts = spath[-1].split('.')[1]
    tmp = pd.read_csv(path)
    msk = tmp['name'].str.startswith('o_')
    tmp = tmp.loc[msk, :].copy()
    tmp['class'] = class_dict[class_metric]
    tmp['task'] = task_dict[task_metric]
    tmp['db'] = db_dict[name_db]
    tmp['dts'] = dts_dict[name_dts]
    tmp['org'] = name_org
    df.append(tmp)
cols = ['class', 'task', 'db', 'org', 'dts', 'name', 'prc', 'rcl', 'f01']
df = pd.concat(df, axis=0).loc[:, cols].reset_index(drop=True)
df['name'] = [mth_dict[n.split('.')[0].replace('o_', '')] for n in df['name']]

# Write
df.to_csv(path_out, index=False)
