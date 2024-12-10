from tqdm import tqdm
import pandas as pd
import requests
import sys
import re


def get_n_pairs(tf_a, tf_b):
    query = f'{tf_a}[Title/Abstract]+AND+{tf_b}[Title/Abstract]'
    base = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    url = f"{base}esearch.fcgi?db=pubmed&term={query}"
    res = requests.get(url).text
    n = int(re.search(r"<Count>(\d+)</Count>", res).group(1))
    return n


# Read
tfs = pd.read_csv(sys.argv[1], sep='\t', header=None)[0].values.astype('U')

# Find matches
df = []
for i in tqdm(range(tfs.size)):
    tf_a = tfs[i]
    for j in range(i + 1, tfs.size):
        tf_b = tfs[j]
        n = get_n_pairs(tf_a, tf_b)
        if n > 0:
            df.append([tf_a, tf_b, n])
df = pd.DataFrame(df, columns=['tf_a', 'tf_b', 'n'])

# Write
df.to_csv(sys.argv[2], sep='\t', index=False, header=False)
