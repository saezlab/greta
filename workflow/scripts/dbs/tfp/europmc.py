import scipy.stats as ss
import pandas as pd
import sys


# Vars
path_single = sys.argv[1]
path_pairs = sys.argv[2]
pval_thr = float(sys.argv[3])
min_odds = float(sys.argv[4])
path_out = sys.argv[5]

# Read
single = pd.read_csv(path_single)
total = single['n'].sum()
single = single.set_index('tf')['n'].to_dict()
pairs = pd.read_csv(path_pairs)

# Compute one-sided Fisher test
df = []
for row in pairs.values:
    tf_a, tf_b, n = row
    only_a = single[tf_a] - n
    only_b = single[tf_b] - n
    backgr = total - (single[tf_a] + single[tf_b])
    s, p = ss.fisher_exact([[n, only_a], [only_b, backgr]], alternative='greater')
    df.append([tf_a, tf_b, s, p])
df = pd.DataFrame(df, columns=['tf_a', 'tf_b', 'stat', 'pval'])
df['padj'] = ss.false_discovery_control(df['pval'])

# Filter
df = df[(df['padj'] < pval_thr) & (df['stat'] > min_odds)].copy()
df['name'] = ['|'.join(sorted([a, b])) for a, b in zip(df['tf_a'], df['tf_b'])]
df[['tf_a', 'tf_b']] = df['name'].str.split('|', expand=True)
df = df.drop(columns='name')

# Save
df.to_csv(path_out, index=False, header=False, sep='\t')
