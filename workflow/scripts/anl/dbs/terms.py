import pandas as pd
import numpy as np
from tqdm import tqdm
import os
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-i','--db_paths', required=True, nargs='+')
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

db_paths = args['db_paths']
path_out = args['path_out']

non_term_dbs = ['CollecTRI', 'ENCODE Blacklist', 'ENCODE CREs', 'Promoters', 'Zhang21', 'phastCons']
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
}
df = []
for db_path in db_paths:
    db_name = os.path.basename(os.path.dirname(db_path))
    db_name = db_dict[db_name]
    task = os.path.basename(os.path.dirname(os.path.dirname(db_path)))
    if db_name not in non_term_dbs:
        if task == 'tfb':
            db = pd.read_csv(db_path, header=None, sep='\t', usecols=[4])[4]
            terms = set()
            for r in tqdm(db):
                terms.update(r.split(','))
            terms = sorted(terms)
        elif task == 'tfm':
            db = pd.read_csv(db_path, sep='\t', header=None, usecols=[1])[1]
            terms = set()
            for r in db:
                terms.update(r.split(','))
            terms = sorted(terms)
        elif task == 'prt':
            db = pd.read_csv(db_path)
            terms = np.sort(db['Tissue.Type'].unique())
        elif 'Catalog' in db_name:
            db = pd.read_csv(db_path, header=None, sep='\t', usecols=[4])[4]
            terms = set()
            for r in tqdm(db):
                r = r.split(',')
                if isinstance(r, str):
                    r = [r]
                for s_r in r:
                    terms.update(s_r.split('|'))
            terms = sorted(terms)
        else:
            raise ValueError('db {db} of task {task} has no defined terms'.format(db=db_name, task=task))
        for term in terms:
            df.append([db_name, term])
df = pd.DataFrame(df, columns=['db_name', 'term'])

# Write
df.to_csv(path_out, index=False)
