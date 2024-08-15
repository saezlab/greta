import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True, nargs='+')
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

df_paths = args['path_input']
path_out = args['path_out']

df = []
for df_path in df_paths:
    tmp = pd.read_csv(df_path)
    df.append(tmp)
df = pd.concat(df)

# Write
df.to_csv(path_out, index=False)
