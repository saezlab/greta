import os
import sys
import pandas as pd
from tqdm import tqdm


tfs = set(pd.read_csv(sys.argv[1], header=None).iloc[:, 0].astype('U'))
mta = pd.read_csv(sys.argv[2], header=None, sep='\t', index_col=0).iloc[:, 0].to_dict()
file_handles = {}
for line in tqdm(sys.stdin):
    if line.startswith('chr'):
        chrm, start, end, tf_ctype = line.strip().split('\t')[:4]
        tf, ctype = tf_ctype.split(':')
        start, end = int(start), int(end)
        if tf in tfs and '_' not in chrm and (end - start) < int(sys.argv[3]):
            ctypes = [mta[c] for c in ctype.split(',') if c in mta]
            if ctypes:
                if tf not in file_handles:
                    file_handles[tf] = open(os.path.join(sys.argv[4], f'{tf}.bed'), 'w')
                file_handles[tf].write(f'{chrm}\t{start}\t{end}\t{tf}\t{",".join(ctypes)}\n')
for tf in file_handles:
    file_handles[tf].close()
