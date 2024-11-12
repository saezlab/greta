import os
import sys
import pandas as pd
from tqdm import tqdm


gdict = pd.read_csv(sys.argv[1]).set_index('id')['symbol'].to_dict()
thr_pval = float(sys.argv[2])
name = os.path.basename(sys.argv[3]).replace('.bed', '')
with open(sys.argv[3], 'w') as f:
    next(sys.stdin)  # skip first line
    for line in tqdm(sys.stdin):
        line = line.strip().split('\t')
        gene, coords, pval = line[1], line[3], float(line[7])
        chrm, start = coords.split('_')[:2]
        start, end = int(start), int(start)
        valid = (pval < thr_pval) and (gene in gdict) and ('_' not in chrm)
        if valid:
            gene = gdict[gene]
            f.write(f'{chrm}\t{start}\t{end}\t{gene}\t{name}\n')
