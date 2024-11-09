import pandas as pd
import os
import sys


df = pd.read_csv(sys.argv[1], sep='\t', usecols=[0, 1, 2, 3, 4, 5], header=None)
tfs = pd.read_csv(sys.argv[2], header=None).values.ravel()
org = sys.argv[1].split(os.sep)[1]
msk_org = df[1] == org
msk_tfs = df[3].isin(tfs)
msk_unc = ~(df[4] == 'Unclassified')
msk = msk_org & msk_tfs & msk_unc
df = df.loc[msk, :].dropna()
df['ctype'] = df[4] + ',' + df[5]
df = df[[0, 3, 'ctype']]
df.to_csv(sys.argv[1], sep='\t', index=False, header=None)
