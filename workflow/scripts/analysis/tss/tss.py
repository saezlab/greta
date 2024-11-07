import pandas as pd
from tqdm import tqdm
import pyranges as pr
import os
import argparse


# Parse args
parser = argparse.ArgumentParser()
parser.add_argument('-t', '--path_tss', required=True, nargs='+')
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()
path_tss = args.path_tss
out_path = args.path_out


# Read
names = []
pr_tss = []
for path in path_tss:
    name = os.path.basename(path).replace('.bed', '')
    tss = pd.read_csv(path, sep='\t', header=None)
    tss.columns = ['Chromosome', 'Start', 'End', 'Name']
    tss = tss.sort_values(['Chromosome', 'Start', 'End', 'Name'])
    tss = pr.PyRanges(tss)
    names.append(name)
    pr_tss.append(tss)


# Find shared genes
genes = set()
genes = genes.union(pr_tss[0].Name)
for tss in pr_tss[1:]:
    genes = genes.intersection(tss.Name)

# Find genomic overlap coef
def overlap_coef_per_gene(gene, tss_a, tss_b):
    ftss_a = tss_a[tss_a.Name == gene].merge()
    ftss_b = tss_b[tss_b.Name == gene].merge()
    if ftss_a.empty or ftss_b.empty:
        raise ValueError('Gene has to be in tss')
    overlap = ftss_a.intersect(ftss_b)
    if overlap.empty:
        return 0.
    else:
        l = overlap.length
        if l == 0:
            return 1
        else:
            return l / min(ftss_a.length, ftss_b.length)


df = []
for i in range(len(names)):
    tss_a = pr_tss[i]
    tss_a = tss_a[tss_a.Name.isin(genes)]
    name_a = names[i]
    for j in range(i + 1, len(names)):
        tss_b = pr_tss[j]
        tss_b = tss_b[tss_b.Name.isin(genes)]
        name_b = names[j]
        print(f"Comparison: {name_a} vs {name_b}")
        for gene in tqdm(list(genes)):
            val = overlap_coef_per_gene(gene, tss_a, tss_b)
            df.append([name_a, name_b, gene, val])
df = pd.DataFrame(df, columns=['tss_a', 'tss_b', 'gene', 'ocoef'])

# Write
df.to_csv(out_path, index=False)
