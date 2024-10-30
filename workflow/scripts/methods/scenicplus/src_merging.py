import pandas as pd
import argparse

# Parse arguments
parser = argparse.ArgumentParser()
parser.add_argument('-e', '--grn_extended',  type=str, help='Path to extended GRN tsv file')
parser.add_argument('-d', '--grn_direct', type=str, help='Path to direct GRN tsv file')
parser.add_argument('-o', '--output', type=str, help='Path to output merged GRN csv file')
args = vars(parser.parse_args())


extended = pd.read_csv(args["grn_extended"], sep="\t")
direct = pd.read_csv(args["grn_direct"], sep="\t")
extended = extended[["TF", "Gene", "rho_TF2G"]]
direct = direct[["TF", "Gene", "rho_TF2G"]]

# Merge the two GRNs
merged = pd.concat([extended, direct])

# Average the values for the same TF-Gene pairs
merged = merged.groupby(["TF", "Gene"]).mean().reset_index()

# Save the merged GRN
merged = merged.rename({"TF": "source", "Gene":"target", "rho_TF2G":"score"}, axis=1)
merged.to_csv(args["output"], sep=",", index=False)

