import pandas as pd
import mudata as mu
import pyranges as pr
import json
import os
import decoupler as dc
import re
import sys

def read_config(path_config='config/config.yaml'):
    import yaml
    with open(path_config, 'r') as file:
        config = yaml.safe_load(file)
    return config

def load_cats(dataset, case):
    with open('config/prior_cats.json') as f:
        cats = json.load(f)
    if (dataset == 'pbmc10k'):
        for i in range(4):
            cats[dataset][str(i)] = cats[dataset]['all'].copy()
    cats = cats[dataset][case]
    return cats

def get_size_tfm(path_db, path_mdata):
    db = pd.read_csv(path_db, sep='\t', header=None)
    genes = mu.read(os.path.join(path_mdata, 'mod', 'rna')).var_names.astype('U')
    db.columns = ['gene', 'ctype']
    dts = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(path_mdata))))
    if dts in ['fakepitupair', 'pitunpair']:
        dts = 'pitupair'
    cats = load_cats(dts, 'all')
    cats = cats[db_name]
    cats = [re.escape(c) for c in cats]
    db = db[db['ctype'].str.contains('|'.join(cats))]
    db = db[db['gene'].astype('U').isin(genes)]
    size = db.shape[0]
    return size

def get_size_tfp(path_db, path_mdata):
    genes = mu.read(os.path.join(path_mdata, 'mod', 'rna')).var_names.astype('U')
    db = pd.read_csv(path_db, sep='\t', header=None)
    db = db[db[0].isin(genes) & db[1].isin(genes)]
    size = db.shape[0]
    return size

def get_size_grn(path_db, path_mdata):
    genes = mu.read(os.path.join(path_mdata, 'mod', 'rna')).var_names.astype('U')
    db = pd.read_csv(path_db)
    db = db[db['source'].astype('U').isin(genes) & db['target'].astype('U').isin(genes)]
    size = db.shape[0]
    return size

def get_size_gst(path_db, path_mdata, thr_pval=0.01, thr_prop=0.2):
    rna = mu.read(os.path.join(path_mdata, 'mod', 'rna'))
    db = pd.read_csv(path_db)
    if 'weight' not in db.columns:
        db['weight'] = 1.
    dc.run_ulm(
        mat=rna,
        net=db,
        use_raw=False,
        verbose=False
    )
    pvals = rna.obsm['ulm_pvals'].copy()
    pvals.loc[:, :] = dc.p_adjust_fdr(pvals.values.ravel()).reshape(pvals.shape)
    acts = rna.obsm['ulm_estimate'].copy()
    hits = ((pvals < thr_pval) & (acts > 0)).sum(0).sort_values(ascending=False) / pvals.shape[0]
    hits = hits[hits > thr_prop].index.values.astype('U')
    size = len(hits)
    return size

def get_size_omics(path_mdata, db):
    genes = mu.read(os.path.join(path_mdata, 'mod', 'rna')).var_names.astype('U')
    peaks = mu.read(os.path.join(path_mdata, 'mod', 'atac')).var_names.astype('U')
    if db == 'gtf':
        size = genes.size
    elif db == 'cretf':
        size = peaks.size
    elif db == 'gcre':
        size = genes.size
    return size

def get_size_genom(path_db, path_mdata, filter_name):
    dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(path_mdata))))
    resource_name = os.path.basename(path_db).replace('.bed', '').replace('.gz', '')
    if dataset in ['fakepitupair', 'pitunpair']:
        dataset = 'pitupair'
    cat_resources = ['chipatlas', 'remap2022', 'unibind', 'gwascatalogue', 'eqtlcatalogue']
    db = pr.read_bed(path_db)
    cats = load_cats(dataset, 'all')
    cats = cats[resource_name]
    if resource_name in cat_resources:
        cats_set = set(cats)
        mask = db.df['Score'].apply(lambda x: any(c.strip() in cats_set for c in x.split(',')))
        db = db[mask]
    genes = mu.read(os.path.join(path_mdata, 'mod', 'rna')).var_names.astype('U')
    peaks = mu.read(os.path.join(path_mdata, 'mod', 'atac')).var_names.astype('U')
    peaks = pd.DataFrame(peaks, columns=['cre'])
    peaks[['Chromosome', 'Start', 'End']] = peaks['cre'].str.split('-', n=2, expand=True)
    peaks = pr.PyRanges(peaks[['Chromosome', 'Start', 'End']])
    db = db.overlap(peaks)
    if filter_name:
        db = db[db.df.Name.astype('U').isin(genes)]
    size = db.df.shape[0]
    return size

def get_size_sss(path_mdata):
    mdata = mu.read(path_mdata)
    size = mdata.obs['celltype'].unique().size
    return size

def get_size_prt(path_mdata, org):
    dataset = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(path_mdata))))
    rna = mu.read(os.path.join(path_mdata, 'mod', 'rna'))
    obs = pd.read_csv(os.path.join(f'dbs/{org}/prt/knocktf', 'meta.csv.gz'), index_col=0)
    if dataset in ['fakepitupair', 'pitunpair']:
        dataset = 'pitupair'
    cats = load_cats(dataset, 'all')
    cats = cats['knocktf']
    cats = [re.escape(c) for c in cats]
    msk = obs['Tissue.Type'].isin(cats) & obs['TF'].isin(rna.var_names) & (obs['logFC'] < -0.5)
    obs = obs.loc[msk, :]
    size = obs.shape[0]
    return size

metrics = {
    'hg38': {
        'genom': {
            'cre': ['blacklist', 'zhang21', 'encode', 'promoters', 'phastcons', 'gwascatalogue'],
            'c2g': ['eqtlcatalogue'],
            'tfb': ['chipatlas', 'remap2022', 'unibind'],
        },
        'pred': {
            'gsets': ['reac', 'hall', 'kegg', 'prog'],
            'omics': ['gtf', 'cretf', 'gcre'],
        },
        'prior': {
            'tfm': ['hpa', 'tfmdb'],
            'tfp': ['intact', 'europmc'],
            'grn': ['collectri'],
        },
        'mech': {
            'tfa': ['knocktf'],
            'sss': ['sss'],
            'prt': ['knocktf'],
        },
    },
    'mm10': {
        'genom': {
            'cre': ['blacklist', 'encode', 'promoters', 'phastcons'],
            'tfb': ['chipatlas', 'remap2022', 'unibind'],
        },
        'pred': {
            'gsets': ['reac', 'hall', 'prog'],
            'omics': ['gtf', 'cretf', 'gcre'],
        },
        'prior': {
            'grn': ['collectri'],
        },
        'mech': {
            'tfa': ['knocktf'],
            'sss': ['sss'],
            'prt': ['knocktf'],
        },
    }
}

config = read_config()
dbs_dict = config['dbs_names']
dts_dict = config['dts_names']
tsk_dict = config['task_names']
cls_dict = config['class_names']

df = []
for dts in config['dts']:
    org = config['dts'][dts]['organism']
    path_mdata = f'dts/{org}/{dts}/cases/all/mdata.h5mu'
    for m_class in metrics[org]:
        for m_task in metrics[org][m_class]:
            for db_name in metrics[org][m_class][m_task]:
                print(dts, m_class, m_task, db_name)
                if m_task == 'tfm':
                    path_db = f'dbs/{org}/{m_task}/{db_name}/{db_name}.tsv.gz'
                    size = get_size_tfm(path_db, path_mdata)
                    label = '# TFs'
                elif m_task == 'tfp':
                    path_db = f'dbs/{org}/{m_task}/{db_name}/{db_name}.tsv.gz'
                    size = get_size_tfp(path_db, path_mdata)
                    label = '# TF pairs'
                elif m_task == 'grn':
                    path_db = f'dbs/{org}/gst/{db_name}.csv.gz'
                    size = get_size_grn(path_db, path_mdata)
                    label = '# Edges'
                elif m_task == 'gsets':
                    path_db = f'dbs/{org}/gst/{db_name}.csv.gz'
                    size = get_size_gst(path_db, path_mdata)
                    label = '# Enriched Sets'
                elif m_task == 'omics':
                    size = get_size_omics(path_mdata, db=db_name)
                    label = '# Features'
                elif m_class == 'genom':
                    path_db = f'dbs/{org}/{m_task}/{db_name}/{db_name}.bed.gz'
                    if m_task == 'cre':
                        filter_name = False
                    else:
                        filter_name = True
                    size = get_size_genom(path_db, path_mdata, filter_name=filter_name)
                    label = '# Regions'
                elif m_task == 'sss':
                    size = get_size_sss(path_mdata)
                    label = '# Cell types'
                elif m_task in ['tfa', 'prt']:
                    size = get_size_prt(path_mdata, org=org)
                    label = '# Experiments'
                df.append([org, dts_dict[dts], cls_dict[m_class], tsk_dict[m_task], dbs_dict[db_name], size, label])
df = pd.DataFrame(df, columns=['org', 'dts', 'class', 'task', 'db', 'size', 'label'])
df.to_csv(sys.argv[1], index=False)
