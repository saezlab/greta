import pandas as pd
import sys

df = pd.read_table(sys.argv[1])
df = df[['TF', 'Region', 'Gene', 'regulation', 'triplet_rank']].groupby(['TF', 'Gene'], as_index=False).mean(numeric_only=True).sort_values('triplet_rank')
df = df.reset_index(drop=True).reset_index(names='rank')
df['score'] = (1 - (df['rank'] / df['rank'].max())) * df['regulation']
df = df[['TF', 'Gene', 'score']]
df.columns = ['source', 'target', 'score']
df['pval'] = 0.01
df.to_csv(sys.argv[2], index=False)
