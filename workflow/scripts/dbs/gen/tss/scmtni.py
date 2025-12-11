import pandas as pd
import sys

path_inp = sys.argv[1]
path_out = sys.argv[2]

bed = pd.read_csv(path_inp, sep='\t', header=None)
bed = bed[[0, 1, 2, 6]].dropna().drop_duplicates()

bed.to_csv(path_out, sep="\t", index=False, header=None, compression="gzip")
