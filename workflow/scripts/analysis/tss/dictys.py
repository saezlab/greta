import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-o', '--path_out', required=True)
parser.add_argument('-i', '--path_input', required=True)
args = parser.parse_args()
out_path = args.path_out
input_path = args.path_input

# Read file
bed = pd.read_csv(input_path, sep='\t', header=None)

# Process columns
bed.columns = ['Chromosome', 'Start', 'End', 'Name', 'score', 'strand']
bed = bed[['Chromosome', 'Start', 'End', 'Name']]
bed['Start'] = bed['Start'] - 1
bed['End'] = bed['End'] - 1

# Save file
bed.to_csv(out_path, sep="\t", index=False, header=None)
