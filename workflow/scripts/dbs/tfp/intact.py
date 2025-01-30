import pandas as pd
import sys


# Read
db = pd.read_csv(sys.argv[1], sep="\t", usecols=['#ID(s) interactor A', 'ID(s) interactor B', 'Confidence value(s)'])
tfs = pd.read_csv(sys.argv[2], header=None)[0].to_list()
pid = pd.read_csv(sys.argv[3])

# Format
p_to_g = pid.set_index('uniprot_id')['symbol'].to_dict()
db = db.rename(columns={
    '#ID(s) interactor A': 'tf_a',
    'ID(s) interactor B': 'tf_b',
    'Confidence value(s)': 'score',
})
db['tf_a'] = db['tf_a'].str.extract(r'uniprotkb:(\w+)')[0].map(p_to_g)
db['tf_b'] = db['tf_b'].str.extract(r'uniprotkb:(\w+)')[0].map(p_to_g)
db['score'] = db['score'].str.extract(r'intact-miscore:(\d+\.\d+)').astype(float)

# Filter
db = db[db['score'] > 0.75].dropna()
db = db[db['tf_a'].isin(tfs) & db['tf_b'].isin(tfs)]
db = db[db['tf_a'] != db['tf_b']].copy()
db['str'] = ['|'.join(sorted([a, b])) for a, b in zip(db['tf_a'], db['tf_b'])]
db = db.drop_duplicates('str').sort_values('score', ascending=False)
db[['tf_a', 'tf_b']] = db['str'].str.split('|', expand=True)
db = db.drop(columns=['str'])

# Write
db.to_csv(sys.argv[4], index=False, header=False, sep='\t')
