import pandas as pd
import sys

expr = pd.read_csv(sys.argv[1], header=0, index_col=0, sep='\t')
blink = pd.read_csv(sys.argv[2], header=0, index_col=0, sep='\t')

namereg = blink.index
nametarget = blink.columns
namet = list(namereg) + [x for x in nametarget if x not in set(namereg)]

expr_sub = expr.loc[namet]
expr_sub = expr_sub.loc[:, expr_sub.sum(axis=0) > 0]
msk = expr.columns.isin(expr_sub.columns)
expr.loc[:, msk].to_csv(sys.argv[3], sep='\t', header=True, index=True)
