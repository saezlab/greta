import argparse
import gzip
import re
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--path_genes', required=True)
parser.add_argument('-i', '--path_gtf', required=True)
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()

gene_set = set(pd.read_csv(args.path_genes, header=None)[0].tolist())
print(f"LINGER gene universe: {len(gene_set)} genes")

records = []
with gzip.open(args.path_gtf, 'rt') as f:
    for line in f:
        fields = line.rstrip('\n').split('\t')
        if fields[2] != 'transcript':
            continue
        gene_name = re.search(r'gene_name "([^"]+)"', fields[8]).group(1)
        records.append([fields[0], int(fields[3]), gene_name])

# mirror `get_TSS_ensembl` in LINGER_tr_fast.py 
df = pd.DataFrame(records, columns=['chr', 'start', 'symbol'])
df = df.groupby(['chr', 'symbol'], as_index=False)['start'].min()
df = df[df['symbol'] != '']       
df['start'] -= 1                   

std_chroms = [f'chr{i}' for i in range(1, 23)] + ['chrX']
df = df[df['symbol'].isin(gene_set) & df['chr'].isin(std_chroms)]
print(f"Matched {df['symbol'].nunique()}/{len(gene_set)} genes")

bed = pd.DataFrame({
    'chr':   df['chr'],
    'start': df['start'],
    'end':   df['start'],   
    'name':  df['symbol']
})
bed = bed.sort_values(['chr', 'start'])
bed.to_csv(args.path_out, sep='\t', index=False, header=False, compression='gzip')
print(f"Written {len(bed)} records to {args.path_out}")