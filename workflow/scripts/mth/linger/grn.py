import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

out_dir = args['out_dir']
path_out = args['path_out'] 

# Load three cell population GRNs
alpha = pd.read_csv(f'{out_dir}/cell_population_TF_RE_binding.txt', sep='\t', index_col=0)
beta  = pd.read_csv(f'{out_dir}/cell_population_cis_regulatory.txt', sep='\t', header=None)
gamma = pd.read_csv(f'{out_dir}/cell_population_trans_regulatory.txt', sep='\t', index_col=0)
beta.columns = ['RE', 'TG', 'beta']

alpha.index.name = 'RE'
alpha_long = alpha.reset_index().melt(id_vars='RE', var_name='TF', value_name='alpha')

gamma.index.name = 'TG'
gamma_long = gamma.reset_index().melt(id_vars='TG', var_name='TF', value_name='gamma')

# Filter to keep top 10% (reduce the combinatorial explosion of later steps)
alpha_long = alpha_long[alpha_long['alpha'] > alpha_long['alpha'].quantile(0.90)]
beta       = beta[beta['beta'] > beta['beta'].quantile(0.90)]
gamma_long = gamma_long[gamma_long['gamma'] > gamma_long['gamma'].quantile(0.90)]

# Build triplets
triplets = pd.merge(alpha_long, beta, on='RE')
triplets = pd.merge(triplets, gamma_long, on=['TF', 'TG'])

# Compute score
triplets['score'] = triplets['alpha'] * triplets['beta'] * triplets['gamma']
triplets = triplets[['TF', 'RE', 'TG', 'score']].sort_values('score', ascending=False)

# Reformat 
triplets.columns = ['source', 'cre', 'target', 'score']
triplets['cre'] = triplets['cre'].str.replace(':', '-')

triplets.to_csv(path_out, index=False)