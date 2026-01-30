import pandas as pd
import sys


df = pd.read_csv(sys.argv[1])

pair = df[df['dts'] == 'Pituitary']
npair = df[df['dts'] == 'Unpaired Pituitary']
fake = df[df['dts'] == 'Synthetic Pituitary']
pvn = (npair.set_index(['class', 'task', 'db', 'name'])['f01'] - pair.set_index(['class', 'task', 'db', 'name'])['f01']).reset_index()
pvf = (fake.set_index(['class', 'task', 'db', 'name'])['f01'] - pair.set_index(['class', 'task', 'db', 'name'])['f01']).reset_index()

mean_pvn = pvn.groupby(['name'])['f01'].mean().sort_values(ascending=True).reset_index()
mean_pvn['type'] = 'Unpaired vs Paired'
mean_pvf = pvf.groupby(['name'])['f01'].mean().sort_values(ascending=True).reset_index()
mean_pvf['type'] = 'Synthethic vs Paired'
mean_df = pd.concat([mean_pvn, mean_pvf])

# Write
mean_df.to_csv(sys.argv[2], index=False)
