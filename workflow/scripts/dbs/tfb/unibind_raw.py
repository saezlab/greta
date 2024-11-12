import os
import sys
import pandas as pd
from tqdm import tqdm


tfs = set(pd.read_csv(sys.argv[1], header=None).iloc[:, 0].astype('U'))
file_handles = {}
for line in tqdm(sys.stdin):
    chrm, start, end, tmp = line.strip().split('\t')[:4]
    tmp = tmp.split('_')
    if len(tmp) == 4:
        _, ctype, tf, _ = tmp 
        start, end = int(start), int(end)
        ctype = ctype.replace('-', ' ').replace(',', ' ').strip()
        tf = tf.strip()
        valid = (tf in tfs) and ('_' not in chrm) and ((start - end) < int(sys.argv[2]))
        if valid:
            if tf not in file_handles:
                file_handles[tf] = open(os.path.join(sys.argv[3], f'{tf}.bed'), 'w')
            file_handles[tf].write(f'{chrm}\t{start}\t{end}\t{tf}\t{ctype}\n')
for tf in file_handles:
    file_handles[tf].close()