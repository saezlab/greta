import pandas as pd
import numpy as np
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--path_rannot', required=True)
parser.add_argument('-b','--samples', required=True, nargs='+')
parser.add_argument('-c', '--path_annot', required=True)
args = vars(parser.parse_args())

path_rannot = args['path_rannot']
samples = args['samples']
path_annot = args['path_annot']

annot = pd.read_csv(path_rannot)
annot = annot[annot['batch'].isin(samples)]
annot['barcode'] = annot['batch'] + '_' + annot['barcode']
annot = annot.set_index('barcode', drop=True)
annot.to_csv(path_annot, header=True)
