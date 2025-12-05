import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import savefigs


path_knocktf = sys.argv[1]
path_hpa = sys.argv[2]
path_out = sys.argv[3]

# knocktf
data = pd.read_csv(os.path.join(path_knocktf, 'diff.csv'), index_col=0)
meta = pd.read_csv(os.path.join(path_knocktf, 'meta.csv'), index_col=0)
fig_ktf, ax = plt.subplots(1, 1, figsize=(4, 3))
sns.histplot(x=meta['logFC'], y=data.std(1), cbar=True, cmap='viridis')
ax.set_xlabel('TF log2FC')
ax.set_ylabel('std log2FC all genes')
ax.axvline(x=-0.5, ls='--', c='gray')

# hpa
df = pd.read_csv(path_hpa, sep='\t', header=None, names=['TFs', 'Contexts'])
df = df.assign(Context=df['Contexts'].str.split(',')).explode('Context')
df['Cell types'] = df['Context'].str.strip()
mat = pd.crosstab(df['Cell types'], df['TFs'])
g = sns.clustermap(mat, cmap='viridis', cbar=False, figsize=(4, 4), dendrogram_ratio=(0.1, 0.1), cbar_pos=None, xticklabels=False, yticklabels=False)
g.fig.suptitle("HPA", y=1.02)
fig_hpa = g.fig

# Write
figs = [fig_ktf, fig_hpa]
savefigs(figs, sys.argv[3])
