from gtfparse import read_gtf
import pandas as pd
import sys

path_inp = sys.argv[1]
path_out = sys.argv[2]

bed = read_gtf(path_inp)
bed = pd.DataFrame(bed, columns=bed.columns)
bed = bed[bed["gene_type"] == "protein_coding"]
bed = bed[['seqname', 'start', 'end', 'gene_name']].dropna()
bed = bed.drop_duplicates()

bed.to_csv(path_out, sep="\t", index=False, header=None, compression="gzip")
