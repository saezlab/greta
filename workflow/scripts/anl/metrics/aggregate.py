import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True, nargs='+')
parser.add_argument('-a','--add_info', required=False, action="store_true")
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

df_paths = args['path_input']
add_info = args['add_info']
path_out = args['path_out']

df = []
for df_path in df_paths:
    tmp = pd.read_csv(df_path)
    if add_info:
        dts = os.path.basename(os.path.dirname(df_path))
        db = os.path.basename(os.path.dirname(os.path.dirname(df_path)))
        task = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(df_path))))
        metric = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(os.path.dirname(df_path)))))
        tmp[['metric', 'task', 'db', 'dts']] = [metric, task, db, dts]
        tmp = tmp[['metric', 'task', 'db', 'dts', 'name', 'prc', 'rcl', 'f01']]
    df.append(tmp)
df = pd.concat(df)

# Write
df.to_csv(path_out, index=False)