import pandas as pd
import numpy as np
import zipfile
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_annot', required=True)
args = vars(parser.parse_args())

path_annot = args['path_annot']

archive = zipfile.ZipFile(path_annot, 'r')
obs = pd.read_csv(archive.open('multiome/snATAC/cells.tsv'), sep='\t')
obs = obs[['barcode', 'sample', 'cluster']].set_index('barcode').rename(columns={'sample': 'batch', 'cluster': 'celltype'})
obs.index.name = None
obs = obs[obs['batch'] != 'D2']
annot = {
    1: 'Fibroblast',
    2: 'Fibroblast-like',
    3: 'Fibroblast-like',
    4: 'Fibroblast-like',
    5: 'Fibroblast-like',
    6: 'Keratinocyte-like',
    7: 'hOSK',
    8: 'xOSK',
    9: 'Intermediate',
    10: 'Partially-reprogrammed',
    11: 'Intermediate',
    12: 'Intermediate',
    13: 'Pre-iPSC',
    14: 'Pre-iPSC',
    15: 'iPSC',
}
obs['celltype'] = [annot[c] for c in obs['celltype']]
obs.to_csv(path_annot)
