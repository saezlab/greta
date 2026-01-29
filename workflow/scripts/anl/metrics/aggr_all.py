import pandas as pd
import glob
import sys

path_out = sys.argv[1]

class_dict = {
    'genom': 'Genomic',
    'mech': 'Mechanistic',
    'pred': 'Predictive',
    'prior': 'Prior',
}

db_dict = {
    'eqtlcatalogue': 'eQTL Catalogue',
    'blacklist': 'ENCODE Blacklist',
    'encode': 'ENCODE CREs',
    'gwascatalogue': 'GWAS Catalog',
    'promoters': 'Promoters',
    'zhang21': 'Zhang21',
    'phastcons': 'phastCons',
    'chipatlas': 'ChIP-Atlas',
    'remap2022': 'ReMap 2022',
    'unibind': 'UniBind',
    'knocktf': 'KnockTF',
    'hpa': 'Human Protein Atlas (HPA)',
    'tfmdb': 'TF-Marker',
    'europmc': 'Europe PMC',
    'intact': 'IntAct',
    'hall': 'Hallmarks',
    'reac': 'Reactome',
    'prog': 'PROGENy',
    'kegg': 'KEGG',
    'gtf': 'Gene ~ TFs',
    'cretf': 'CRE ~ TFs',
    'gcre': 'Gene ~ CREs',
    'collectri': 'CollecTRI',
    'sss': 'Boolean rules',
    'knocktf': 'KnockTF',
}

task_dict = {
    'cre': 'CREs',
    'c2g': 'CRE to Gene',
    'tfb': 'TF binding',
    'gsets': 'Gene sets',
    'omics': 'Omics',
    'tfm': 'TF markers',
    'tfp': 'TF pairs',
    'grn': 'Reference GRN',
    'tfa': 'TF scoring',
    'sss': 'Steady state simulation',
    'prt': 'Perturbation forecasting',
}

dts_dict = {
    'lung': 'Lung',
    'embryo': 'Embryo',
    'skin': 'Skin',
    'heart': 'Heart',
    'pbmc10k': 'PBMC',
    'epalate': 'Palate',
    'kidney': 'Kidney',
    'eye': 'Eye',
    'pitupair': 'Pituitary',
    'fakepitupair': 'Synthethic Pituitary',
    'pitunpair': 'Unpaired Pituitary',
    'breast': 'Breast',
    'reprofibro': 'Reprog. Fibro',
    'brain': 'Brain',
}

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
    'scdori': 'scDoRI',
    'scenic': 'SCENIC',
    'scenicplus': 'SCENIC+',
    'scgpt': 'scGPT',
    'scmtni': 'scMTNI',
    'spearman': 'Spearman',
}

lst_paths = glob.glob('anl/metrics/*/*/*/*.*.all.scores.csv')
df = []
for path in lst_paths:
    spath = path.split('/')
    class_metric = spath[2]
    task_metric = spath[3]
    name_db = spath[4]
    name_org = spath[-1].split('.')[0]
    name_dts = spath[-1].split('.')[1]
    tmp = pd.read_csv(path)
    tmp['class'] = class_dict[class_metric]
    tmp['task'] = task_dict[task_metric]
    tmp['db'] = db_dict[name_db]
    tmp['dts'] = dts_dict[name_dts]
    tmp['org'] = name_org
    df.append(tmp)
cols = ['class', 'task', 'db', 'org', 'dts', 'name', 'prc', 'rcl', 'f01']
df = pd.concat(df, axis=0).loc[:, cols].reset_index(drop=True)
df['name'] = [mth_dict[n.split('.')[0].replace('o_', '')] for n in df['name']]

# Write
df.to_csv(path_out, index=False)
