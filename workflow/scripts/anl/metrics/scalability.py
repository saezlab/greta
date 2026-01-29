import pandas as pd
import sys

path_inp = sys.argv[1]  # 'anl/stab/pitupair.ovc.csv'
path_out = sys.argv[2]

uses_gpu = ['dictys', 'scdori', 'scgpt']

mth_dict = {
    'celloracle': 'CellOracle',
    'collectri': 'CollecTRI',
    'crema': 'CREMA',
    'dictys': 'Dictys',
    'directnet': 'DirectNet',
    'dorothea': 'DoRothEA',
    'figr': 'FigR',
    'granie': 'GRaNIE',
    'grnboost': 'GRNBoost2',
    'hummus': 'HuMMuS',
    'inferelator': 'Inferelator',
    'pando': 'Pando',
    'pearson': 'Pearson',
    'random': 'Random',
    'scdori': 'scDORi',
    'scenic': 'SCENIC',
    'scenicplus': 'SCENIC+',
    'scgpt': 'scGPT',
    'scmtni': 'scMTNI',
    'spearman': 'Spearman',
}

df = pd.read_csv(path_inp)
cols = ['mth', 'seed', 'h', 'gb']
df = df[df['cat'] == 'full'][cols].drop_duplicates(['mth', 'h', 'gb']).sort_values(['mth', 'seed'])
df = df.groupby(['mth'])[['h', 'gb']].mean().reset_index()
df['use_gpu'] = df['mth'].isin(uses_gpu)
df['mth'] = [mth_dict[n] for n in df['mth']]
df.to_csv(path_out, index=False)
