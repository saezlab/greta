import scanpy as sc
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_h5ad', required=True)
parser.add_argument('-a','--path_atac', required=True)
parser.add_argument('-o','--path_annot', required=True)
args = vars(parser.parse_args())

path_h5ad = args['path_h5ad']
path_atac = args['path_atac']
path_annot = args['path_annot']

# load atac data
adata = ad.read_h5ad(path_atac)

# filter for left ventricle
atac_obs = adata.obs[adata.obs['region'] == 'LV']

# rename indexes and neural cell type for ATAC
atac_obs = atac_obs[['combinedID', 'cell_type']].rename(columns={'combinedID': 'batch_atac', 'cell_type': 'celltype_atac'})
atac_obs.index.name = None
atac_obs.index = atac_obs.index.to_series().apply(lambda x: x.rsplit('_', 1)[-1])
atac_obs['celltype_atac'] = [c if c != 'Neural cell' else 'Neuronal cell' for c in atac_obs['celltype_atac']]
atac_obs = atac_obs.reset_index().rename(columns={'index': 'index'})

# load rna data
adataRNA = ad.read_h5ad(path_h5ad)

# filter for left ventricle
rna_obs = adataRNA.obs[adataRNA.obs['region'] == 'LV']

# rename indexes for RNA
rna_obs = rna_obs[['combinedID', 'scANVI_predictions']].rename(columns={'combinedID': 'batch_rna', 'scANVI_predictions': 'celltype_rna'})
rna_obs.index.name = None
rna_obs = rna_obs[rna_obs['batch_rna'] != 'na']
rna_obs['batch_rna'] = rna_obs['batch_rna'].cat.remove_categories(['na'])
rna_obs.index = rna_obs.index.to_series().apply(lambda x: x.rsplit('_', 1)[-1])
rna_obs = rna_obs.reset_index().rename(columns={'index': 'index'})

# merge rna and atac annotation
obs = pd.merge(rna_obs, atac_obs, how='inner') 
obs = obs[(obs['celltype_rna'] == obs['celltype_atac']) & (obs['batch_rna'] == obs['batch_atac'])] 
obs = obs[['index', 'batch_rna', 'celltype_rna']].rename(columns={'batch_rna': 'batch', 'celltype_rna': 'celltype'}).set_index('index')  # Rename
obs.index.name = None 

obs.to_csv(path_annot)