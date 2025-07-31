# %%
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt

# %%
input_file = '../../../dts/hpa.tsv'
output_file = '../../../dts/hpa.pdf'

# %%
df = pd.read_csv(input_file, sep='\t', header=None, names=['TF', 'Contexts'])

# %%
df

# %%
# Expand contexts into long format
df_long = df.assign(Context=df['Contexts'].str.split(',')).explode('Context')
df_long['Context'] = df_long['Context'].str.strip()  # Remove extra whitespace

# %%
# Create wide binary matrix
df_wide = pd.crosstab(df_long['TF'], df_long['Context'])

# %%
# Plot clustered binary matrix
clustermap = sns.clustermap(df_wide, cmap='Greys', cbar=False, figsize=(10, 8))

# Save as PDF
clustermap.fig.savefig(output_file, format='pdf', bbox_inches='tight')

plt.close()


