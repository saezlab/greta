import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-h','--path_h5ad', required=True)
parser.add_argument('-a','--path_annot', required=True)
args = vars(parser.parse_args())

path_h5ad = args['path_h5ad']
path_annot = args['path_annot']

# Load data
adata = sc.read_h5ad(path_h5ad)

obs = adata.obs
obs = obs[['combinedID', 'scANVI_predictions']].rename(columns={'combinedID': 'batch', 'scANVI_predictions': 'celltype'})
obs.index.name = None
obs = obs[obs['batch'] != 'na']
obs['batch'] = obs['batch'].cat.remove_categories(['na'])
obs.index = obs.index.to_series().apply(lambda x: x.rsplit('_', 1)[-1])

obs.to_csv(path_annot)