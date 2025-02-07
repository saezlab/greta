import pandas as pd
import pyranges as pr
from tqdm import tqdm
import numpy as np
import sys

def set_ocoef(a, b):
    min_s = min(len(a), len(b))
    if min_s == 0:
        return np.nan
    else:
        inter = len(a & b)
        return inter / min_s

# Overlap TFs
tf_dbs = dict()

df = pd.read_csv('dbs/hg38/prt/knocktf/meta.csv', index_col=0)
tf_dbs['knocktf'] = set(df['TF'])

df = pd.read_csv('dbs/hg38/tfm/hpa/hpa.tsv', sep='\t', header=None)
tf_dbs['hpa'] = set(df[0])

df = pd.read_csv('dbs/hg38/tfm/tfmdb/tfmdb.tsv', sep='\t', header=None)
tf_dbs['tfmdb'] = set(df[0])

df = pd.read_csv('dbs/hg38/tfp/europmc/europmc.tsv', sep='\t', header=None)
tf_dbs['europmc'] = set(df[0]) | set(df[1])

df = pd.read_csv('dbs/hg38/tfp/intact/intact.tsv', sep='\t', header=None)
tf_dbs['intact'] = set(df[0]) | set(df[1])

df = pd.read_csv('dbs/hg38/tfb/chipatlas/chipatlas.bed', sep='\t', usecols=[3], header=None)
tf_dbs['chipatlas'] = set(df[3])

df = pd.read_csv('dbs/hg38/tfb/remap2022/remap2022.bed', sep='\t', usecols=[3], header=None)
tf_dbs['remap2022'] = set(df[3])

df = pd.read_csv('dbs/hg38/tfb/unibind/unibind.bed', sep='\t', usecols=[3], header=None)
tf_dbs['unibind'] = set(df[3])

ocf = []
keys = list(tf_dbs.keys())
for i, k_a in enumerate(keys):
    a = tf_dbs[k_a]
    for k_b in keys[i + 1:]:
        b = tf_dbs[k_b]
        ocf.append(['tf', k_a, k_b, set_ocoef(a, b)])

# Overlap Genes
g_dbs = dict()

df = pd.read_csv('dbs/hg38/gst/hall.csv')
g_dbs['hall'] = set(df['target'])

df = pd.read_csv('dbs/hg38/gst/kegg.csv')
g_dbs['kegg'] = set(df['target'])

df = pd.read_csv('dbs/hg38/gst/prog.csv')
g_dbs['prog'] = set(df['target'])

df = pd.read_csv('dbs/hg38/gst/reac.csv')
g_dbs['reac'] = set(df['target'])

df = pd.read_csv('dbs/hg38/prt/knocktf/diff.csv', index_col=0)
g_dbs['knocktf'] = set(df.columns)

df = pd.read_csv('dbs/hg38/c2g/eqtlcatalogue/eqtlcatalogue.bed', sep='\t', header=None, usecols=[3])
g_dbs['eqtlcatalogue'] = set(df[3])

keys = list(g_dbs.keys())
for i, k_a in enumerate(keys):
    a = g_dbs[k_a]
    for k_b in keys[i + 1:]:
        b = g_dbs[k_b]
        ocf.append(['gene', k_a, k_b, set_ocoef(a, b)])

# Overlap bp
bp_dbs = dict()
bp_dbs['chipatlas'] = pr.read_bed('dbs/hg38/tfb/chipatlas/chipatlas.bed')
bp_dbs['remap2022'] = pr.read_bed('dbs/hg38/tfb/remap2022/remap2022.bed')
bp_dbs['unibind'] = pr.read_bed('dbs/hg38/tfb/unibind/unibind.bed')
bp_dbs['blacklist'] = pr.read_bed('dbs/hg38/cre/blacklist/blacklist.bed')
bp_dbs['encode'] = pr.read_bed('dbs/hg38/cre/encode/encode.bed')
bp_dbs['gwascatalogue'] = pr.read_bed('dbs/hg38/cre/gwascatalogue/gwascatalogue.bed')
bp_dbs['phastcons'] = pr.read_bed('dbs/hg38/cre/phastcons/phastcons.bed')
bp_dbs['promoters'] = pr.read_bed('dbs/hg38/cre/promoters/promoters.bed')
bp_dbs['zhang21'] = pr.read_bed('dbs/hg38/cre/zhang21/zhang21.bed')
bp_dbs['eqtlcatalogue'] = pr.read_bed('dbs/hg38/c2g/eqtlcatalogue/eqtlcatalogue.bed')

tfb = ['chipatlas', 'remap2022', 'unibind']
keys = list(bp_dbs.keys())
for i, k_a in enumerate(keys):
    a = bp_dbs[k_a]
    for k_b in keys[i + 1:]:
        b = bp_dbs[k_b]
        if (k_a in tfb) and (k_b in tfb):
            inters = set(bp_dbs[k_a].df['Name']) & set(bp_dbs[k_b].df['Name'])
            mean_oc = []
            for inter in tqdm(inters):
                i_db_a = a[a.df['Name'] == inter]
                i_db_b = b[b.df['Name'] == inter]
                overlap = i_db_a.set_intersect(i_db_b)
                mean_oc.append(overlap.length / np.min([i_db_a.length, i_db_b.length]))
            ocf.append(['bp', k_a, k_b, np.mean(mean_oc)])
        else:
            overlap = a.set_intersect(b)
            ocf.append(['bp', k_a, k_b, overlap.length / np.min([a.length, b.length])])
ocf = pd.DataFrame(ocf, columns=['type', 'db_a', 'db_b', 'ocoeff'])

# Write
ocf.to_csv(sys.argv[1], index=False)
