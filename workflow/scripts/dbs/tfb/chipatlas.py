import sys
for line in sys.stdin:
    line = line.replace('\n', '').split('\t')
    chrm, start, end, tf, ctype = line[0], line[1], line[2], line[3], line[4]
    ctype = ','.join(sorted(set(ctype.split(','))))
    print(f'{chrm}\t{start}\t{end}\t{tf}\t{ctype}')
