from inferelator import inferelator_workflow
from inferelator import MPControl
import glob
import os
import sys
import pandas as pd
import numpy as np

path_tmp = sys.argv[1]
ncores = int(sys.argv[2])
path_out = sys.argv[3]

if os.path.isfile(os.path.join(path_tmp, 'empty.txt')):
    grn = pd.DataFrame(columns=['source', 'cre', 'target', 'score'])
    grn.to_csv(path_out, index=False)
    os._exit(0)

# Run
if __name__ == '__main__':
    MPControl.set_multiprocess_engine("multiprocessing")
    MPControl.client.processes = ncores
    MPControl.connect()
worker = inferelator_workflow(regression="amusr", workflow="multitask")
worker.set_file_paths(
    input_dir=path_tmp,
    output_dir=os.path.join(path_tmp, 'network'),
    tf_names_file='tfs.tsv',
)
worker.set_network_data_flags(use_no_gold_standard=True)
worker.set_file_properties(metadata_handler='nonbranching')
path_adata = os.path.join(path_tmp, 'adatas')
for file in glob.glob(os.path.join(path_adata, '*.h5ad')):
    name = os.path.basename(file).replace('.h5ad', '')
    task = worker.create_task(
        task_name=name,
        input_dir=path_tmp,
        tf_names_file='tfs.tsv',
        priors_file='prior.tsv',
        workflow_type="tfa",
    )
    file = file.replace(path_tmp + '/', '')
    task.set_expression_file(h5ad=file)
worker.set_run_parameters(num_bootstraps=5, random_seed=42, use_numba=True)
worker.run()

# Format grn
path_net = os.path.join(path_tmp, 'network', 'network.tsv.gz')
df = pd.read_csv(path_net, sep='\t')
df = df[df['var.exp.median'] > 0]
df['score'] = np.sign(df['beta.sign.sum']) * df['combined_confidences']
df = df[['regulator', 'target', 'score']].rename(columns={'regulator': 'source', 'target': 'gene'})
df = df[df['score'] != 0.]
p2g = pd.read_csv(os.path.join(path_tmp, 'p2g.csv'))
grn = pd.merge(df, p2g[['cre', 'gene']], on='gene', how='inner')
grn = grn[['source', 'cre', 'gene', 'score']].rename(columns={'gene': 'target'})
grn.to_csv(path_out, index=False)
