import pandas as pd
import numpy as np
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--inpath_annot', required=True)
parser.add_argument('-b', '--outpath_annot', required=True)
parser.add_argument('-c', '--path_samples', required=True)
args = vars(parser.parse_args())

inpath_annot = args['inpath_annot']
outpath_annot = args['outpath_annot']
path_samples = args['path_samples']


# Get samples 

samples = [x for x in os.listdir(path_samples) if x.endswith(".h5")]
samples = [x.split('_')[-1] for x in samples]
samples = [x.split('.')[0] for x in samples]

annot = pd.read_csv(inpath_annot)

annot[annot['batch'].isin(samples)].to_csv(outpath_annot)
