import os
import sys
import pandas as pd


tfs = set(pd.read_csv(sys.argv[1], header=None).values.ravel().astype('U'))
for line in sys.stdin:
    line = line.replace('\n', '').split('\t')[:4]
    chrm, start, end, tmp = line[0], line[1], line[2], line[3]
    tmp = tmp.split('_')
    if len(tmp) == 4:
        _, ctype, tf, _ = tmp
        ctype = ctype.replace('-', ' ').replace(',', ' ').strip()
        tf = tf.strip()
        start, end = int(start), int(end)
        valid = (tf in tfs) and ('_' not in chrm) and ((start - end) < int(sys.argv[2]))
        if valid:
            with open(os.path.join(sys.argv[3], f'{tf}.bed'), 'a') as f:
                f.write(f'{chrm}\t{start}\t{end}\t{tf}\t{ctype}\n')
