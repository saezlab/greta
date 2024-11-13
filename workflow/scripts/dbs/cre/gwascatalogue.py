import pandas as pd
import numpy as np
from tqdm import tqdm
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--inp_path', required=True)
args = vars(parser.parse_args())

inp_path = args['inp_path']

# Read tsv input
df = pd.read_csv(inp_path, sep='\t', dtype={9: 'str', 12: 'str', 23: 'str', 26: 'str'})

# Remove nans
df = df[~df['CHR_POS'].isna()]
df = df[~df['SNP_ID_CURRENT'].isna()]
df = df[~df['MAPPED_TRAIT_URI'].isna()]

# Drop one special case with multiple snp ids
df = df[~df['SNP_ID_CURRENT'].astype(str).str.contains(';')]

# Split urls and obtain key, take care of multiple terms separated by commas
df['MAPPED_TRAIT_URI'] = [', '.join([x.split('/')[-1] for x in url.split(',')]) for url in df['MAPPED_TRAIT_URI']]

# Exctracts the risk allele and sets anything else from ATGC to unknown
str_alleles = []
bases = np.array(['A', 'T', 'G', 'C'])
for snp in tqdm(df['STRONGEST SNP-RISK ALLELE']):
    snp = snp.split('-')[-1].upper()
    has_bases = np.all(np.isin([l for l in snp], bases))
    if has_bases and snp != '':
        str_alleles.append(snp)
    else:
        str_alleles.append('?')
df['STRONGEST SNP-RISK ALLELE'] = str_alleles
df['CHR_POS_2'] = df['CHR_POS'].copy()

# Subset by important cols
cols = ['CHR_ID', 'CHR_POS', 'CHR_POS_2', 'STRONGEST SNP-RISK ALLELE',
        'P-VALUE', 'MAPPED_TRAIT', 'MAPPED_TRAIT_URI', 'PUBMEDID']
df = df[cols]

# Transform to correct data types
df['CHR_ID'] = 'chr' + df['CHR_ID'].astype(str)
df['CHR_POS'] = df['CHR_POS'].astype(int)
df['CHR_POS_2'] = df['CHR_POS_2'].astype(int)
df['P-VALUE'] = df['P-VALUE'].astype(float)
df['MAPPED_TRAIT'] = df['MAPPED_TRAIT'].astype(str)
df['MAPPED_TRAIT_URI'] = df['MAPPED_TRAIT_URI'].astype(str)
df['PUBMEDID'] = df['PUBMEDID'].astype(str)

# Summarize when multiple p-values are given
df = df.groupby(list(df.columns[df.columns != 'P-VALUE'])).mean(numeric_only=True).reset_index()

# Rename and sort
df = df.rename(columns={
    'CHR_ID': 'chr_id',
    'CHR_POS': 'chr_start',
    'CHR_POS_2': 'chr_end',
    'STRONGEST SNP-RISK ALLELE': 'eff_allele',
    'MAPPED_TRAIT': 'trait_name',
    'MAPPED_TRAIT_URI': 'trait_uri',
    'PUBMEDID': 'pubmedid',
    'P-VALUE': 'pval'
})

# Save
df = df[['chr_id', 'chr_start', 'chr_end', 'eff_allele', 'trait_name']]
df['trait_name'] = df['trait_name'].str.strip()
df.to_csv(inp_path, index=False, header=None, sep='\t')
