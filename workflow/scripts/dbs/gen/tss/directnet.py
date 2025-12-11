import pandas as pd
import gzip, re, sys
input_gtf = sys.argv[1]
output_bed = sys.argv[2]
tss = ""
with gzip.open(input_gtf, 'rt') as f:
    for line in f:
        chrom, _, typ, start, end, _, strand, _, body = line.strip().split('\t')
        start, end = int(start) - 1, int(end) - 1
        gene = re.search(r'gene_id "([^"]+)"', body).group(1)
        if typ == 'transcript':
            if strand == '+':
                tss += f'{chrom}\t{start}\t{start}\t{gene}\n'
            elif strand == '-':
                tss += f'{chrom}\t{end}\t{end}\t{gene}\n'
with gzip.open(output_bed, 'rt') as f:
    f.write(tss)
