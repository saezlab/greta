import logging
import torch # type: ignore
import numpy as np
from torch.utils.data import DataLoader, TensorDataset # type: ignore
from pathlib import Path
from sklearn.preprocessing import OneHotEncoder
from utils import aggregate_grn_max_val


from scdori import ( # type: ignore
    #trainConfig,
    load_scdori_inputs,
    set_seed,
    scDoRI,
    load_best_model,
    compute_atac_grn_activator_with_significance,
    compute_atac_grn_repressor_with_significance,
    compute_significant_grn,
    get_latent_topics,
    get_tf_expression
)

def wrapper_scdori_grns(trainConfig):
    logger = logging.getLogger(__name__)
    logging.basicConfig(level=trainConfig.logging_level)

    logger.info("Starting scDoRI downstream analysis")
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


    model = load_best_model(
        model, Path(trainConfig.weights_folder_grn) / "best_scdori_best_eval.pth", device
    )
    # creating dataloader for all cells
    n_cells = rna_metacell.n_obs
    indices = np.arange(n_cells)

    all_dataset = TensorDataset(torch.from_numpy(indices))
    all_dataset_loader = DataLoader(
        all_dataset, batch_size=trainConfig.batch_size_cell_prediction, shuffle=False
    )

    # get scDoRI latent embedding (topics)
    scdori_latent = get_latent_topics(
        model,
        device,
        all_dataset_loader,
        rna_metacell,
        atac_metacell,
        num_cells,
        tf_indices,
        onehot_batch,
    )
    np.save(
        Path(trainConfig.weights_folder_scdori) / "scdori_latent.npy",
        scdori_latent
    )
    # adding scDoRI embedding to the anndata object
    rna_metacell.obsm["X_scdori"] = scdori_latent

    grn_act_atac = compute_atac_grn_activator_with_significance(
        model, device, cutoff_val=0.05, outdir=Path(trainConfig.data_dir, "grn"), num_permutations=trainConfig.num_permutations
    )
    # ATAC based GRN for repressors
    grn_rep_atac = compute_atac_grn_repressor_with_significance(
        model, device, cutoff_val=0.05, outdir=Path(trainConfig.data_dir, "grn"), num_permutations=trainConfig.num_permutations
    )
    # calculating TF-expression per topic
    # either from scdori model weights or from true data
    # using true expression here
    tf_normalised = get_tf_expression(
        "True",
        model,
        device,
        all_dataset_loader,
        rna_metacell,
        atac_metacell,
        num_cells,
        tf_indices,
        onehot_batch,
        trainConfig,
    )
    # compute final GRNs which use the significant ATAC based GRNs derived above
    grn_act, grn_rep = compute_significant_grn(
        model,
        device,
        cutoff_val_activator=0.05,
        cutoff_val_repressor=0.05,
        tf_normalised=tf_normalised.detach().cpu().numpy(),
        outdir=Path(trainConfig.data_dir, "grn"),
    )
    # Aggregate GRN activator and repressor matrices
    grn_long = aggregate_grn_max_val(grn_act, grn_rep, rna_metacell)
    grn_long.to_csv(
        trainConfig.grn_file_out, index=False
    )

    logger.info("=== GRN pipeline completed successfully ===")
    return
