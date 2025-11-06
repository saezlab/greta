## ToDo: Add a wrapping function that is callable from the wrapper

import logging
import torch # type: ignore
import numpy as np
from torch.utils.data import DataLoader, TensorDataset # type: ignore
from sklearn.model_selection import train_test_split
from pathlib import Path
from sklearn.preprocessing import OneHotEncoder


from scdori import ( # type: ignore
    #trainConfig,
    load_scdori_inputs,
    save_model_weights,
    set_seed,
    scDoRI,
    train_scdori_phases,
    train_model_grn,
    initialize_scdori_parameters,
    load_best_model,
)
def wrapper_scdori_training(trainConfig):
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=trainConfig.logging_level)

    logger.info("Starting scDoRI pipeline with integrated GRN.")
    set_seed(trainConfig.random_seed)

    device = torch.device("cuda:0" if torch.cuda.is_available() else "cpu")
    logger.info(f"Using device: {device}")

    rna_metacell, atac_metacell, gene_peak_dist, insilico_act, insilico_rep = (
        load_scdori_inputs(trainConfig)
    )
    gene_peak_fixed = gene_peak_dist.clone()
    gene_peak_fixed[gene_peak_fixed > 0] = 1  # mask for peak-gene links based on distance

    # computing indices of genes which are TFs and setting number of cells per metacell ( set to 1 for single cell data)
    rna_metacell.obs["num_cells"] = 1
    rna_metacell.var["index_int"] = range(rna_metacell.shape[1])
    tf_indices = rna_metacell.var[rna_metacell.var.gene_type == "TF"].index_int.values
    num_cells = rna_metacell.obs.num_cells.values.reshape((-1, 1))

    batch_col = trainConfig.batch_col
    rna_metacell.obs["batch"] = rna_metacell.obs[batch_col].values
    atac_metacell.obs["batch"] = atac_metacell.obs[batch_col].values
    # obtaining onehot encoding for technical batch,

    enc = OneHotEncoder(handle_unknown="ignore")
    enc.fit(rna_metacell.obs["batch"].values.reshape(-1, 1))

    onehot_batch = enc.transform(rna_metacell.obs["batch"].values.reshape(-1, 1)).toarray()
    enc.categories_

    # 2) Make small train/test sets
    n_cells = rna_metacell.n_obs
    # Setting epochs dynamically based on dataset and batch size
     
    indices = np.arange(n_cells)
    train_idx, eval_idx = train_test_split(indices, test_size=0.2, random_state=42)
    train_dataset = TensorDataset(torch.from_numpy(train_idx))
    train_loader = DataLoader(
        train_dataset, batch_size=trainConfig.batch_size_cell, shuffle=True
    )

    eval_dataset = TensorDataset(torch.from_numpy(eval_idx))
    eval_loader = DataLoader(
        eval_dataset, batch_size=trainConfig.batch_size_cell, shuffle=False
    )
    num_genes = rna_metacell.n_vars
    num_peaks = atac_metacell.n_vars

    num_tfs = insilico_act.shape[1]

    num_batches = onehot_batch.shape[1]
    model = scDoRI(
        device=device,
        num_genes=num_genes,
        num_peaks=num_peaks,
        num_tfs=num_tfs,
        num_topics=trainConfig.num_topics,
        num_batches=num_batches,
        dim_encoder1=trainConfig.dim_encoder1,
        dim_encoder2=trainConfig.dim_encoder2,
    ).to(device)

    initialize_scdori_parameters(
        model,
        gene_peak_dist.to(device),
        gene_peak_fixed.to(device),
        insilico_act=insilico_act.to(device),
        insilico_rep=insilico_rep.to(device),
        phase="warmup",
    )
    model = train_scdori_phases(
        model,
        device,
        train_loader,
        eval_loader,
        rna_metacell,
        atac_metacell,
        num_cells,
        tf_indices,
        onehot_batch, 
        trainConfig,
    )

    # saving the model weight correspoinding to final epoch where model stopped training
    save_model_weights(model, Path(trainConfig.weights_folder_scdori), "scdori_final")

    # loading the best checkpoint from Phase 1
    model = load_best_model(
        model, Path(trainConfig.weights_folder_scdori) / "best_scdori_best_eval.pth", device
    )
    initialize_scdori_parameters(
        model,
        gene_peak_dist,
        gene_peak_fixed,
        insilico_act=insilico_act,
        insilico_rep=insilico_rep,
        phase="grn",
    )
    # train Phase 2 of scDoRI model, TF-gene links are learnt in this phase and used to reconstruct gene-expression profiles
    model = train_model_grn(
        model,
        device,
        train_loader,
        eval_loader,
        rna_metacell,
        atac_metacell,
        num_cells,
        tf_indices,
        onehot_batch,
        trainConfig,
    )
    # saving the model weight correspoinding to final epoch where model stopped training
    save_model_weights(model, Path(trainConfig.weights_folder_grn), "scdori_final")

    logger.info("=== Training pipeline completed successfully ===")
    return
