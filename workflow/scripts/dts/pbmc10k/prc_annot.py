import os
import snapatac2 as snap
from snapatac2.datasets import _datasets, datasets
from pathlib import Path
import pandas as pd
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_annot', required=True)
args = vars(parser.parse_args())

path_annot = args['path_annot']

# Change default cache dir
_datasets = datasets()
_datasets.path = Path('/tmp/')

# Download
rna = snap.read(snap.datasets.pbmc10k_multiome(modality='RNA', type='h5ad'), backed=None)

# Extract annot
rna.obs['batch'] = 'smpl'
annot = rna.obs.rename(columns={'cell_type': 'celltype'})[['batch', 'celltype']]
annot.index.name = None
annot.index = ['smpl_' + i.replace('-1', '') for i in annot.index]

# Write
annot.to_csv(path_annot)
