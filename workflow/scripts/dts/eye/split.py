import gzip
import pandas as pd
import glob
import os
import argparse
import sys
from tqdm import tqdm


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_ann', required=True)
parser.add_argument('-s','--samples', nargs='+', required=True)
args = vars(parser.parse_args())

path_ann = args['path_ann']
samples = args['samples']

# Read
obs = pd.read_csv(path_ann, index_col=0)
obs = obs[obs['batch'].isin(samples)]
obs.to_csv(path_ann)
bar2sam = dict()
for i in obs.index:
    sam, bar = i.split('_')
    bar2sam[bar] = sam
# Split into sample files
dirname = os.path.dirname(path_ann)
filenames = glob.glob(os.path.join(dirname, 'GSM5866073_*.gz'))
dict_files = {sample: open(os.path.join(dirname, f'{sample}.frags.tsv'), 'wt') for sample in samples}
for fname in tqdm(filenames):
    with gzip.open(fname, "rt") as i_f:
        for line in i_f:
            chrom, start, end, bar, count = line.rstrip().split('\t')
            bar = bar.split('-')[0]
            if bar in bar2sam:
                sam = bar2sam[bar]
                bar = f'{sam}_{bar}'
                dict_files[sam].write(f'{chrom}\t{start}\t{end}\t{bar}\t{count}\n')
            else:
                continue
# Close all new files
for k in dict_files:
    dict_files[k].close()
# Sort and overwrite
filenames = glob.glob(os.path.join(dirname, '*.frags.tsv'))
for fname in tqdm(filenames):
    bed = pd.read_csv(fname, sep='\t', header=None)
    bed = bed.sort_values([0, 1])
    bed.to_csv(fname, index=False, header=False, sep='\t')
