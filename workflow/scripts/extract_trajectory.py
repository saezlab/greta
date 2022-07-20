import pandas as pd
import numpy as np
import muon as mu
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-c','--celltypes', required=True)
parser.add_argument('-o','--path_output', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
celltypes = args['celltypes']
path_output = args['path_output']

# Read object
mdata = mu.read(inp)

# Process celltypes
celltypes = celltypes.split(';')

# Filter
mdata = mdata[np.isin(mdata.obs['celltype'], celltypes)]

# Save
mdata.write(path_output)

