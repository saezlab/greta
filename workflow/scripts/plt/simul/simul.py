import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import sys
import os
sys.path.insert(0, os.path.abspath(os.path.join(os.path.dirname(__file__), '..')))
from utils import read_config, savefigs


path_topo = sys.argv[1]
path_eval = sys.argv[2]
path_out = sys.argv[3]

# Topo
df = pd.read_csv(path_topo)
names = df['name_a'].unique()
s_mat = df[['name_a', 'name_b', 'tf_oc']].pivot(index='name_a', columns='name_b', values='tf_oc').loc[names, names]
e_mat = df[['name_a', 'name_b', 'edge_oc']].pivot(index='name_a', columns='name_b', values='edge_oc').loc[names, names]
t_mat = df[['name_a', 'name_b', 'target_oc']].pivot(index='name_a', columns='name_b', values='target_oc').loc[names, names]

fig_topo, axes = plt.subplots(1, 3, figsize=(9, 3), sharey=True)
axes = axes.ravel()

ax = axes[0]
sns.heatmap(s_mat, cmap='Purples', cbar=False, annot=True, fmt=".1f", ax=ax, vmin=0)
ax.set_title('TFs')
ax.set_xlabel('')
ax.set_ylabel('')

ax = axes[1]
sns.heatmap(e_mat, cmap='Purples', cbar=False, annot=True, fmt=".1f", ax=ax, vmin=0)
ax.set_title('Edges')
ax.set_xlabel('')
ax.set_ylabel('')

ax = axes[2]
sns.heatmap(t_mat, cmap='Purples', cbar=False, annot=True, fmt=".1f", ax=ax, vmin=0)
ax.set_title('Genes')
ax.set_xlabel('')
ax.set_ylabel('')

# Eval
df = pd.read_csv(path_eval)
config = read_config()
palette = config['colors']['nets']

fig_eval, ax = plt.subplots(1, 1, figsize=(3, 3))
sns.barplot(data=df, x='f01', y='name', hue='name', palette=palette, ax=ax)
ax.set_xlabel(r'$F_{01}$')
ax.set_ylabel('')

# Write
figs = [fig_topo, fig_eval]
savefigs(figs, path_out)
