import os
import sys
import pandas as pd
from tqdm import tqdm
from concurrent.futures import ProcessPoolExecutor


mta = pd.read_csv(sys.argv[1], sep='\t', header=None)
mta['smpl'] = mta[0] + '.' + mta[1]
mta = mta.set_index('smpl')[2].to_dict()

file_data = {}

for line in tqdm(sys.stdin):
    chrm, start, end, gene, smpl = line.strip().split('\t')
    start, end = int(start), int(end)
    ctype = mta[smpl]
    
    if gene not in file_data:
        file_data[gene] = ""
    file_data[gene] += f'{chrm}\t{start}\t{end}\t{gene}\t{ctype}\n'


def write_gene_file(gene, lines, output_dir):
    with open(os.path.join(output_dir, f'{gene}.bed'), 'w') as f:
        f.writelines(lines)


with ProcessPoolExecutor(max_workers=32) as executor:
    futures = {executor.submit(write_gene_file, gene, lines, sys.argv[2]): gene for gene, lines in file_data.items()}
    for future in tqdm(futures, total=len(futures)):
        future.result()