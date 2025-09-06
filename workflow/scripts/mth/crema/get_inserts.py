import pyranges as pr
import pandas as pd
import os
import sys
from tqdm import tqdm
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f','--path_frgs', required=True, nargs='+')
parser.add_argument('-g','--path_gen', required=True)
parser.add_argument('-d','--path_dir', required=True)
args = vars(parser.parse_args())

path_frgs = args['path_frgs']
path_gen = args['path_gen']
path_dir = args['path_dir']

path_ann = os.path.join(path_dir, 'ann.csv')
ann = pd.read_csv(path_ann, index_col=0)
ann = set(ann.index)
genome = pr.read_bed(path_gen)

def get_inserts(path_frg, ann, genome):
    # Filter by ann
    fr = pr.read_bed(path_frg)
    msk = fr.df['Name'].isin(ann)
    fr = fr[msk]
    # Filter by genome
    fr = fr.overlap(genome)
    # Create inserts
    starts = pd.DataFrame()
    starts['Chromosome'] = fr.Chromosome
    starts['Start'] = fr.Start - 6
    starts['End'] = fr.Start + 5
    starts['Name'] = fr.Name
    starts['Score'] = fr.Score
    ends = pd.DataFrame()
    ends['Chromosome'] = fr.Chromosome
    ends['Start'] = fr.End - 5
    ends['End'] = fr.End + 6
    ends['Name'] = fr.Name
    ends['Score'] = fr.Score
    fr = pd.concat([starts, ends])
    return fr

fr = []
for path_frg in tqdm(path_frgs):
    fr.append(get_inserts(path_frg=path_frg, ann=ann, genome=genome))
fr = pd.concat(fr)
fr = fr.sort_values(['Chromosome', 'Start', 'End'])
# Write
path_out = os.path.join(path_dir, f'inserts.tsv')
fr.to_csv(path_out, sep='\t', header=None, index=False)
