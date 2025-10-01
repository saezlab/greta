# for python 3.10
# pip install circe-py pyarrow "git+https://github.com/cantinilab/HuMMuS.git@dask_update#subdirectory=hummuspy" arboreto
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-f', '--path_data', required=True)
parser.add_argument('-r', '--rna_layer', required=True)
parser.add_argument('-a', '--atac_layer', required=True)
parser.add_argument('-o', '--organism', required=True)
parser.add_argument('-w', '--wsize', required=True)
parser.add_argument('-d', '--tmp_dir', required=False, default='/tmp/hummus_tmp/')
parser.add_argument('-c', '--num_workers', required=True)

args = vars(parser.parse_args())

path_data = args['path_data']
path_rna_layer = args['rna_layer']
path_atac_layer = args['atac_layer']
organism = args["organism"]
wsize = int(args['wsize'])
tmp_dir = args['tmp_dir']
n_cpu = int(args['num_workers'])

if organism in ["hg38", "hg19", "GRCh38", "GRCh37"]:
    species = "human"
elif organism in ["mm10", "mm9", "GRCm38", "GRCm39"]:
    species = "mouse"
else:
    species = "unknown"
    print(f"organism {organism} not handled, please use human or mouse references")


os.environ['NUMBA_CACHE_DIR'] = tmp_dir
import numba
import scanpy as sc
import muon as mu
from hummuspy.loader import load_tfs  # import tf names for grnboost2
from arboreto.algo import _prepare_input
from arboreto.core import (EARLY_STOP_WINDOW_LENGTH, RF_KWARGS, SGBM_KWARGS,
                           infer_partial_network, to_tf_matrix)
import joblib
from tqdm import tqdm
import pandas as pd
import pathlib
from typing import List, Literal
import pandas as pd
import circe as ci # atac layer


def compute_rna_network(
        df_exp_mtx: pd.DataFrame,
        tf_names: List[str],
        temp_dir: pathlib.Path,
        method: Literal['GBM', 'RF'] = 'GBM',
        n_cpu: int = 1,
        seed: int = 666) -> pd.DataFrame:
    """
    # Inspired from https://github.com/aertslab/scenicplus/blob/main/src/scenicplus/TF_to_gene.py

    Calculate TF-to-gene relationships using either Gradient Boosting Machine (GBM)
    or Random Forest (RF) regression.

    It is a wrapper around the `infer_partial_network` function from the arboreto package, similarly to GRNBoost2.
    It uses joblib to parallelize the inference of the relationships for each target gene.
    It returns a DataFrame with the TF-to-gene relationships and their importance scores.

    Parameters
    ----------
    df_exp_mtx : pd.DataFrame
        DataFrame of shape (n_cells, n_genes) containing the expression matrix.
        Rows are cells and columns are genes.
    tf_names : List[str]
        List of transcription factor names to consider as potential regulators.
    temp_dir : pathlib.Path
        Path to a temporary directory to store intermediate files during parallel processing.
    method : Literal['GBM', 'RF'], optional
        Method to use for regression. Either 'GBM' for Gradient Boosting Machine or '
        'RF' for Random Forest. Default is 'GBM'.
    n_cpu : int, optional
        Number of CPU cores to use for parallel processing. Default is 1.
    seed : int, optional
        Random seed for reproducibility. Default is 666.

    Returns
    -------
    pd.DataFrame
        DataFrame with columns ['tf', 'target', 'importance'] representing the TF-to-g
        ene relationships and their importance scores.
    
    """

    if(method == 'GBM'):
        method_params = [
            'GBM',      # regressor_type
            SGBM_KWARGS  # regressor_kwargs
        ]
    elif(method == 'RF'):
        method_params = [
            'RF',       # regressor_type
            RF_KWARGS   # regressor_kwargs
        ]

    exp_mtx, gene_names, tf_names = _prepare_input(
        expression_data = df_exp_mtx, gene_names = None, tf_names = tf_names)
    tf_matrix, tf_matrix_gene_names = to_tf_matrix(
        exp_mtx,  gene_names, tf_names)
            
    print('Calculating TF-to-gene importance')
    if temp_dir is not None:
        if type(temp_dir) is str:
            temp_dir = pathlib.Path(temp_dir)
        if not temp_dir.exists():
            Warning(f"{temp_dir} does not exist, creating it.")
            os.makedirs(temp_dir)
        
    TF_to_genes = joblib.Parallel(
        n_jobs = n_cpu,
        temp_folder = temp_dir)(
            joblib.delayed(infer_partial_network)(
                target_gene_name = gene,
                target_gene_expression = exp_mtx[:, gene_names.index(gene)],
                regressor_type = method_params[0],
                regressor_kwargs = method_params[1],
                tf_matrix = tf_matrix,
                tf_matrix_gene_names = tf_matrix_gene_names,
                include_meta = False,
                early_stop_window_length = EARLY_STOP_WINDOW_LENGTH,
                seed = seed)
            for gene in tqdm(
                gene_names, 
                total=len(gene_names), 
                desc=f'Running using {n_cpu} cores'))

    adj = pd.concat(TF_to_genes).sort_values(by='importance', ascending=False)

    return adj


# import tf names list
if species == "human":
    tfs_list = load_tfs("human_tfs_r_hummus")
elif species == "mouse":
    tfs_list = load_tfs("mouse_tfs_r_hummus")
else:
    tfs_list = load_tfs("species not handled")

if __name__ == "__main__":
    # import data - h5mu
    data = mu.read_h5mu(path_data)

    print(data)
    data  = mu.MuData({
        "rna": data["rna"][:,:],
        "atac": data["atac"][:,:]
    })

    print(data["rna"].to_df())
    print(list(data["rna"].var_names[:50]))
    grnboost2_network = compute_rna_network(
        data["rna"][:,:].to_df(),
        tf_names=tfs_list[:],
        temp_dir=pathlib.Path(tmp_dir),
        n_cpu = n_cpu
        )
    grnboost2_network.to_csv(
        path_rna_layer,
        index=False
    )

    # create atac layer
    atac = ci.add_region_infos(data["atac"], sep=(":", "-"))
    #   atac = ci.metacells.compute_metacells(atac)

    ci.compute_atac_network(atac, window_size=wsize, distance_constraint=wsize//2)
    circe_network = ci.extract_atac_links(atac)
    circe_network = circe_network[circe_network["score"] > 0]
    circe_network.to_csv(path_atac_layer, index=False)
