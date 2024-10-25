#!/usr/bin/env python
# coding: utf-8

# In[ ]:


import pandas as pd
import requests
import argparse

# Initiate args
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()
out_path = args.path_out 


# Download bed file
url = "https://github.com/morris-lab/CellOracle/blob/e5ae78e93272da7d772378e60ae6cd4602f24be6/celloracle/motif_analysis/tss_ref_data/hg38_tss_info.bed?raw=true"
response = requests.get(url)

# Convert bed to csv file
from io import StringIO

bed_data = StringIO(response.text)  # Convert the response content into a file-like object
col_names = ['Chromosome', 'Start', 'End', 'Name', 'Pos', 'Strand']  # Specify the correct column names
bed_df = pd.read_csv(bed_data, sep='\t', names=col_names)
bed_df = bed_df[["Chromosome", "Start", "End", "Name"]]

# Save file
bed_df.to_csv(out_path, index=False)

