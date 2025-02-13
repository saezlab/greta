import pandas as pd
import numpy as np
import argparse
import os


def read_eval(m_path):
    db_name = os.path.basename(os.path.dirname(m_path))
    task = os.path.basename(os.path.dirname(os.path.dirname(m_path)))
    metric = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname(m_path))))
    case = os.path.basename(m_path).replace('.scores.csv', '')
    df = pd.read_csv('anl/metrics/{0}/{1}/{2}/{3}.scores.csv'.format(metric, task, db_name, case)).sort_values('f01', ascending=False)
    df[['pre', 'p2g', 'tfb', 'mdl']] = df['name'].str.split('.', n=4, expand=True)
    df = df[~df['pre'].str.startswith('o_')]
    df = df.reset_index(drop=True).reset_index(names='rank')
    df['fixed'] = [np.unique(n.split('.')).size == 1 for n in df['name']]
    return metric, task, db_name, case, df


def test_rank(df):
    import decoupler as dc
    steps = ['pre', 'p2g', 'tfb', 'mdl']
    mthds = df['pre'].unique()
    net = []
    sts = []
    for step in steps:
        sts.append(df.groupby([step], as_index=False)['f01'].mean().rename(columns={step: 'name'}).assign(stp=step))
        for mth in mthds:
            for name in df[df[step] == mth]['name']:
                net.append(['{0}.{1}'.format(step, mth), name])
    net = pd.DataFrame(net, columns=['source', 'target'])
    sts = pd.concat(sts)
    res = dc.get_gsea_df(
        df=df.dropna().set_index('name'),
        stat='f01',
        net=net,
        times=1000
    )
    res['padj'] = np.where(res['ES'] > 0, res['FDR p-value'], 1)
    res[['stp', 'name']] = res['Term'].str.split('.', n=2, expand=True)
    res = res[['stp', 'name', 'padj']]
    res = pd.merge(res, sts, how='left', on=['stp', 'name'])
    return res


parser = argparse.ArgumentParser()
parser.add_argument('-m', '--path_mtr', nargs='+', required=True)
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()

# Test each metric-database
df = []
for m_path in args.path_mtr:
    metric, task, db_name, case, m_df = read_eval(m_path)
    m_df = test_rank(m_df)
    m_df[['metric', 'task', 'db', 'case']] = metric, task, db_name, case
    df.append(m_df)
df = pd.concat(df)
df = df[['metric', 'task', 'db', 'stp', 'name', 'case', 'padj', 'f01']]
df = df.sort_values(['metric', 'task', 'db', 'stp', 'name'])

# Write
df.to_csv(args.path_out, index=False)
