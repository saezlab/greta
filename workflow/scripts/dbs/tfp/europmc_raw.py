from tqdm import tqdm
import pandas as pd
import requests
import re
import time
import sys


def do_query(query):
    base = 'https://www.ebi.ac.uk/europepmc/webservices/rest/search'
    url = f"{base}?query={query}&format=json"
    res = requests.get(url)
    while res.status_code != 200:
        print(url, flush=True)
        time.sleep(1)
        res = requests.get(url)
    n = int(res.json()['hitCount'])
    return n


def get_n_pairs(tf_a, tf_b):
    query = f'(TITLE:"{tf_a}"+OR+ABSTRACT:"{tf_a}")+AND+(TITLE:"{tf_b}"+OR+ABSTRACT:"{tf_b}")'
    return do_query(query)


def get_n_single(tf):
    query = f'(ABSTRACT:"{tf}"+OR+TITLE:"{tf}")'
    return do_query(query)


# Read args
path_tfs = sys.argv[1]
min_chars = int(sys.argv[2])
min_n = int(sys.argv[3])
path_single = sys.argv[4]
path_pairs = sys.argv[5]

# Open tfs
tfs = pd.read_csv(path_tfs, sep='\t', header=None)[0].values.astype('U')

# Find unique tfs with enough publications (min_n) and characters (min_chars)
single_tfs = []
for tf in tqdm(tfs):
    if len(tf) > min_chars:
        single_tfs.append([tf, get_n_single(tf)])
single_tfs = pd.DataFrame(single_tfs, columns=['tf', 'n']).sort_values('n')
single_tfs = single_tfs[single_tfs['n'] > min_n]
tfs = single_tfs['tf'].sort_values().unique()
single_tfs.to_csv(path_single, index=False)

# Find pairs 
df = []
for i in tqdm(range(tfs.size)):
    tf_a = tfs[i]
    for j in range(i + 1, tfs.size):
        tf_b = tfs[j]
        n = get_n_pairs(tf_a, tf_b)
        if n > 0:
            df.append([tf_a, tf_b, n])
df = pd.DataFrame(df, columns=['tf_a', 'tf_b', 'n'])
df.to_csv(path_pairs, index=False)
