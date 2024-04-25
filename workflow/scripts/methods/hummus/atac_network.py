import numpy as np
import pandas as pd
import anndata as ad
import scanpy as sc
import atacnet as an
import mudata as mu

# Parameters
distance_threshold = snakemake.params['distance_threshold']

# Load data
#atac = ad.read_csv(snakemake.input["atac"], delimiter='\t')
#print(atac.var_names)
#
#
## Add region infos
#an.add_region_infos(atac, sep=('_', '_'))
#print(atac)
#print(atac.X)
#print(atac.X.sum(0).min(), atac.X.sum(0).min())
#
#atac = an.metacells.compute_metacells(atac,
 #                              k = snakemake.params["number_cells_per_clusters"])



mudata = mu.read_h5mu(snakemake.input["mudata"])
atac = mudata["atac"].copy()
del mudata
print(atac)
# Add region infos
an.add_region_infos(atac, sep=('-', '-'))


# Compute network
an.compute_atac_network(
    atac,
    window_size=distance_threshold,
    unit_distance=1000,
    distance_constraint=distance_threshold/2,
    n_samples=100,
    n_samples_maxtry=500,
    max_alpha_iteration=100,
)

# Save results
peak_layer = an.extract_atac_links(
    atac,
    columns=['peak1', 'peak2', 'score'])

peak_layer.to_csv(snakemake.output[0], sep='\t', index=False)

