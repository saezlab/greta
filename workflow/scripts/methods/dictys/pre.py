import argparse, os, sys
import pandas as pd
import numpy as np
import mudata as md


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_input', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_input = args['path_input']
path_out = args['path_out']

# Read
mdata = mu.read(path_input)

# Write
mdata.write(path_out)