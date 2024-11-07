import pandas as pd
import requests
from io import StringIO
import argparse

# Initiate args
parser = argparse.ArgumentParser()
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()
out_path = args.path_out 

# Download bed file
url = "https://github.com/morris-lab/CellOracle/blob/e5ae78e93272da7d772378e60ae6cd4602f24be6/celloracle/motif_analysis/tss_ref_data/hg38_tss_info.bed?raw=true"
response = requests.get(url)
bed = pd.read_csv(StringIO(response.text), sep='\t', header=None)[[0, 1, 2, 3]].dropna().sort_values([0, 1, 2])

# Save file
bed.to_csv(out_path, sep="\t", index=False, header=False)
