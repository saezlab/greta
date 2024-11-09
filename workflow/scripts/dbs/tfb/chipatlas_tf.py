import sys
import os
import re
import pandas as pd


tf = os.path.basename(sys.argv[1]).replace('.bed', '')
meta = pd.read_csv(sys.argv[2], sep='\t', header=None).set_index(0)
pattern = r'ID=(.*?);'
for line in sys.stdin:
    if line.startswith('chr'):
        line = line.replace('\n', '').split('\t')
        chrm, start, end, sample_id = line[0], line[1], line[2], line[3]
        sample_id = re.search(pattern, sample_id).group(1)
        if (sample_id in meta.index) and ('_' not in chrm):
            m_tf = meta.loc[sample_id, 1]
            ctype = meta.loc[sample_id, 2]
            start, end = int(start), int(end)
            if (m_tf == tf) and ((start - end) < int(sys.argv[3])):
                print(f'{chrm}\t{start}\t{end}\t{tf}\t{ctype}')
