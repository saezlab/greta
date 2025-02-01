from pycisTopic.cistopic_class import create_cistopic_object
import scipy.sparse as scs
import pandas as pd
import polars as pl
import mudata as mu
import tempfile
import scipy
import os
import yaml
import pickle
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-b','--path_bl', required=True)
parser.add_argument('-n','--n_topics', required=True)
parser.add_argument('-c','--n_cores', required=True)
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-o','--path_out', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
path_bl = args['path_bl']
n_topics = int(args['n_topics'])
n_cores = int(args['n_cores'])
out_dir = args['out_dir']
path_out = args['path_out']

# Read metadata
mdata = mu.read(path_mdata)
obs = mdata.obs.copy()
obs['celltype'] = obs['celltype'].str.replace(r'\s+', '_', regex=True)

# Write gex
os.makedirs(out_dir, exist_ok=True)
path_rna = os.path.join(out_dir, 'rna.h5ad')
mdata.mod['rna'].write(path_rna)

# Create cistopic object
cistopic_obj = create_cistopic_object(
    fragment_matrix=scs.csr_matrix(mdata.mod['atac'].layers['counts'].T),
    cell_names=mdata.obs_names,
    region_names=mdata.mod['atac'].var_names.str.replace('-', ':', 1),
    path_to_blacklist=path_bl,
    project='object',
    tag_cells=False,
    split_pattern="_",
)

# Write raw data
path_bmat = os.path.join(out_dir, "bool_mat.mtx")
path_cres = os.path.join(out_dir, "cres.txt")
path_brcs = os.path.join(out_dir, "brcs.txt")
scipy.io.mmwrite(path_bmat, scipy.sparse.csr_matrix(cistopic_obj.binary_matrix != 0.))
with open(path_cres, 'w') as f:
    for name in cistopic_obj.region_names:
        f.write(f"{name}\n")

with open(path_brcs, 'w') as f:
    for name in cistopic_obj.cell_names:
        f.write(f"{name}\n")

# Create corpus
command = f"pycistopic topic_modeling mallet create_corpus -i {path_bmat} -o {path_bmat}.mallet"
print(command, flush=True)
os.system(command)

# Compute topics
path_topics = os.path.join(out_dir, "topics")
os.makedirs(path_topics, exist_ok=True)
command = f"pycistopic topic_modeling mallet run \
-i {path_bmat}.mallet \
-o {path_topics}/topic_ \
-t {n_topics} \
-p {n_cores}"
print(command, flush=True)
os.system(command)

# Get topic model
command = f"pycistopic topic_modeling mallet stats \
-i {path_bmat} \
-c {path_brcs} \
-r {path_cres} \
-o {path_topics}/topic_ \
-t {n_topics}"
print(command, flush=True)
os.system(command)

# Load model
path_model = os.path.join(path_topics, f'topic_.{n_topics}_topics.model.pkl')
with open(path_model, 'rb') as pickle_file:
    model = pickle.load(pickle_file)
cistopic_obj.add_LDA_model(model)

# Write object
path_obj = os.path.join(out_dir, "cistopic_obj.pkl")
with open(path_obj, 'wb') as f:
    pickle.dump(cistopic_obj, f)

# Run scenicplus preprocessing
command = f"scenicplus prepare_data prepare_GEX_ACC \
--cisTopic_obj_fname {path_obj} \
--GEX_anndata_fname {path_rna} \
--out_file {path_out} \
--do_not_use_raw_for_GEX_anndata"
print(command, flush=True)
os.system(command)

