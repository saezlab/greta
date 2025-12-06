import pandas as pd
import sys


df = pd.read_table(sys.argv[1])
if df.shape[0] > 0:
    df = df[df['regulation'] != 0.]
    #df = df[['TF', 'Region', 'Gene', 'regulation', 'triplet_rank']].groupby(['TF', 'Gene'], as_index=False).mean(numeric_only=True).sort_values('triplet_rank')
    df = df[['TF', 'Region', 'Gene', 'regulation', 'triplet_rank']].sort_values('triplet_rank')
    df = df.reset_index(drop=True).reset_index(names='rank')
    df['score'] = (1 - (df['rank'] / df['rank'].max())) * df['regulation']
    df = df[['TF', 'Region', 'Gene', 'score']]
    df.columns = ['source', 'cre', 'target', 'score']
    df['pval'] = 0.01
else:
    df = pd.DataFrame(columns=['source', 'cre', 'target', 'score', 'pval'])
df.to_csv(sys.argv[2], index=False)
