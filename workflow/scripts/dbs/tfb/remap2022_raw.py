import os
import sys
import pandas as pd


tfs = set(pd.read_csv(sys.argv[1], header=None).values.ravel().astype('U'))
mta = pd.read_csv(sys.argv[2], header=None, sep='\t').set_index(0)[1].to_dict()
for line in sys.stdin:
    if line.startswith('chr'):
        line = line.replace('\n', '').split('\t')[:4]
        chrm, start, end, tf_ctype = line[0], line[1], line[2], line[3]
        tf, ctype = tf_ctype.split(':')
        start, end = int(start), int(end)
        valid = (tf in tfs) and ('_' not in chrm) and ((start - end) < int(sys.argv[3]))
        if valid:
            ctypes = []
            for c in ctype.split(','):
                if c in mta:
                    ctypes.append(mta[c])
            if len(ctypes) > 0:
                ctypes = ','.join(ctypes)
                with open(os.path.join(sys.argv[4], f'{tf}.bed'), 'a') as f:
                    f.write(f'{chrm}\t{start}\t{end}\t{tf}\t{ctypes}\n')
