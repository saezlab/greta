import sys


meta = pd.read_csv(sys.argv[1], sep='\t', header=None)
meta['id'] = meta[0] + '.' + meta[1]
meta = meta.set_index('id')[2].to_dict()
for line in sys.stdin:
    line = line.replace('\n', '').split('\t')
    chrm, start, end, gene, ctypes = line[0], line[1], line[2], line[3], line[4]
    ctypes = [meta[c] for c in ctypes.split(',')]
    ctypes = ','.join(sorted(set(ctypes)))
    print(f'{chrm}\t{start}\t{end}\t{gene}\t{ctypes}')
