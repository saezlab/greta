#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
from unipressed import IdMappingClient
import time
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-t','--tfs_path', required=True)
parser.add_argument('-i','--db_path', required=True)
parser.add_argument('-o','--out_path', required=True)
args = vars(parser.parse_args())

tfs_path = args['tfs_path']
db_path = args['db_path']
path_out = args['out_path']


# Load data
# We load tf list and intact database
tfs = pd.read_csv(tfs_path, header=None)
db = pd.read_csv(db_path, sep="\t")

# Extract TF list
tfs_list = tfs[0].to_list()

# We map gene names to uniprotkb
request = IdMappingClient.submit(
    source="GeneCards", dest="UniProtKB", ids=tfs_list
)
time.sleep(1)

# Retrieve the results and convert to a DataFrame
results = list(request.each_result())
df = pd.DataFrame(results)

# We extract UniProt IDs from the database
db['Interactor A UniProt'] = db['#ID(s) interactor A'].str.extract(r'uniprotkb:(\w+)')
db['Interactor B UniProt'] = db['ID(s) interactor B'].str.extract(r'uniprotkb:(\w+)')

# We filter db for TF presents.
db_uni = db[
    db['Interactor A UniProt'].isin(df['to']) &
    db['Interactor B UniProt'].isin(df['to'])
]

# Merge with df to get 'TF A' and 'TF B' names
merged = (
    db_uni
    .merge(df[['from', 'to']], left_on='Interactor A UniProt', right_on='to', how='left')
    .rename(columns={'from': 'TF1'})
    .drop(columns='to')
    .merge(df[['from', 'to']], left_on='Interactor B UniProt', right_on='to', how='left')
    .rename(columns={'from': 'TF2'})
    .drop(columns='to')
)

# Filter out self-interactions and duplicates (where A-B is the same as B-A)
final_filtered = (
    merged[merged['TF1'] != merged['TF2']]
    .assign(sorted_pair=lambda x: x.apply(lambda row: tuple(sorted([row['TF1'], row['TF2']])), axis=1))
    .drop_duplicates(subset='sorted_pair')
    .drop(columns='sorted_pair')
)

# Filter by confidence score and select only 'Gene A' and 'Gene B'
final = (
    final_filtered
    .assign(Confidence_Score=lambda x: x['Confidence value(s)'].str.extract(r'intact-miscore:(\d+\.\d+)').astype(float))
    [lambda x: x['Confidence_Score'] > 0.9]   # Change confidence score.
    [['TF1', 'TF2']]
)


final.to_csv(path_out)

