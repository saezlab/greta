#!/usr/bin/env python
# coding: utf-8

# In[12]:


import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--path_out', required=True)
parser.add_argument('-i', '--path_input', required=True)
args = parser.parse_args()
out_path = args.path_out
input_path = args.path_input

# Read file
bed_df = pd.read_csv(input_path, sep='\t', header=None)

# Process columns
bed_df.columns = ['Chromosome', 'Start', 'End', 'Name', 'score', 'strand']
bed_df = bed_df[['Chromosome', 'Start', 'End', 'Name']]

# Save file
bed_df.to_csv(out_path, sep="\t")
