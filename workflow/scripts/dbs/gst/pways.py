import decoupler as dc
import pandas as pd
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--path_inp', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_reac = args['path_reac']
path_hall = args['path_hall']
path_kegg = args['path_kegg']
path_tfs = args['path_tfs']
path_prg = args['path_prg']
path_out_hall = args['path_out_hall']
path_out_kegg = args['path_out_kegg']
path_out_prg = args['path_out_prg']
path_out_reac = args['path_out_reac']

# Process hallmark
hall = dc.read_gmt(path_hall)
hall['source'] = hall['source'].str.replace('HALLMARK_', '')

# Process kegg
kegg = dc.read_gmt(path_kegg)
kegg['source'] = kegg['source'].str.replace('KEGG_', '')

# Process progeny
prg = pd.read_csv(path_prg)
prg = prg.rename(columns={'gene': 'target', 'pathway': 'source', 'p.value': 'pval'})
prg = prg[['source', 'target', 'weight', 'pval']]
prg = prg[prg['pval'] < 0.05]
prg = prg.sort_values(['source', 'pval'])
prg = prg.rename(columns={'source': 'pathway', 'target': 'gene'})

# Process reactome
reac = dc.read_gmt(path_reac)
reac['source'] = reac['source'].str.replace('REACTOME_', '')

# Write
kegg.to_csv(path_out_kegg, index=False)
prg.to_csv(path_out_prg, index=False)
hall.to_csv(path_out_hall, index=False)
reac.to_csv(path_out_reac, index=False)
