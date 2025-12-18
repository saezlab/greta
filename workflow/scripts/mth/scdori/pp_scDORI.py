import logging
from pathlib import Path
import numpy as np
import pandas as pd
import muon as mu
import os
import urllib.request

from scdori.pp import ( # type: ignore
    #ppConfig,
    compute_gene_peak_distance_matrix,
    compute_hvgs_and_tfs,
    compute_in_silico_chipseq,
    compute_motif_scores,
    create_dir_if_not_exists,
    create_extended_gene_bed,
    create_metacells,
    download_genome_references,
    filter_protein_coding_genes,
    intersect_cells,
    keep_promoters_and_select_hv_peaks,
    load_gtf,
    load_motif_database,
    load_anndata,
    remove_mitochondrial_genes,
    run_bedtools_intersect,
    save_processed_datasets,
)
def wrapper_scdori_preprocessing(ppConfig):

    logger = logging.getLogger(__name__)

    logging.getLogger().setLevel(ppConfig.logging_level)

    logger.info("=== Starting multi-ome preprocessing pipeline ===")

    data_dir = Path(ppConfig.data_dir)
    genome_dir = Path(ppConfig.genome_dir)
    motif_dir = Path(ppConfig.motif_directory)
    out_dir = data_dir / Path(ppConfig.output_subdir_name)

    create_dir_if_not_exists(genome_dir)
    create_dir_if_not_exists(motif_dir)
    create_dir_if_not_exists(out_dir)
    create_dir_if_not_exists(ppConfig.weight_dir)

    download_genome_references(
        genome_dir=genome_dir,
        species=ppConfig.species,
        assembly=ppConfig.genome_assembly,
        gtf_url=ppConfig.gtf_url,
        chrom_sizes_url=ppConfig.chrom_sizes_url,
        fasta_url=ppConfig.fasta_url,
    )

    # Checking if mudata object was provided (instead of two separate anndata objects)
    if ppConfig.mudata_file_name.endswith(".h5mu"):
        logger.info("Detected mudata object, using it for preprocessing.")
        mdata = mu.read(ppConfig.mudata_file_name)
        mdata.push_obs()
        data_rna = mdata.mod["rna"].copy()
        data_atac = mdata.mod["atac"].copy()
    elif ppConfig.rna_adata_file_name.endswith(".h5ad"):
        logger.info("Detected separate RNA and ATAC anndata objects.")
        data_rna, data_atac = load_anndata(
            data_dir, ppConfig.rna_adata_file_name, ppConfig.atac_adata_file_name
        )
    else:
        raise ValueError(
            "Please provide either a mudata object or two separate anndata objects."
        )
    data_rna, data_atac = intersect_cells(data_rna, data_atac)

    data_rna = remove_mitochondrial_genes(
        data_rna, mito_prefix=ppConfig.mitochondrial_prefix
    )
    gtf_file = genome_dir / "annotation.gtf"
    gtf_df = load_gtf(gtf_file)
    data_rna = filter_protein_coding_genes(data_rna, gtf_df)

    motif_path = motif_dir / f"{ppConfig.motif_database}_{ppConfig.species}.meme"

    # Downloading custom motif database from scdoris github
    if not os.path.exists(motif_path):
        urllib.request.urlretrieve(ppConfig.motif_url + os.path.basename(motif_path), motif_path)
    
    tf_names_all = []
    with open(motif_path) as f:
        for line in f:
            if line.startswith("MOTIF"):
                parts = line.strip().split()
                if len(parts) >= 3:
                    tf_name = parts[2].split("_")[0].strip("()").strip()
                    tf_names_all.append(tf_name)
    tf_names_all = sorted(list(set(tf_names_all)))


    ## Checking if number of HVGs and TFs to select is not more than available
    if ppConfig.num_genes > data_rna.shape[1]:
        logger.warning(
            f"Requested number of genes ({ppConfig.num_genes}) is more than available genes ({data_rna.shape[1]}). Setting to maximum available."
        )
        ppConfig.num_genes = data_rna.shape[1]
    if ppConfig.num_tfs > len(tf_names_all):
        logger.warning(
            f"Requested number of TFs ({ppConfig.num_tfs}) is more than available TFs ({len(tf_names_all)}). Setting to maximum available."
        )
        ppConfig.num_tfs = len(tf_names_all)

    data_rna, final_genes, final_tfs = compute_hvgs_and_tfs(
        data_rna=data_rna,
        tf_names=tf_names_all,
        user_genes=ppConfig.genes_user,
        user_tfs=ppConfig.tfs_user,
        num_genes=ppConfig.num_genes,
        num_tfs=ppConfig.num_tfs,
        min_cells=ppConfig.min_cells_per_gene,
    )
    ## Check whether any barcode lacks expression of all TFs
    barcode_mask = data_rna[:,data_rna.var['gene_type'] == 'TF'].X.sum(1) != 0
    data_rna = data_rna[barcode_mask].copy()
    data_atac = data_atac[barcode_mask].copy()
    chrom_sizes_path = genome_dir / f"{ppConfig.genome_assembly}.chrom.sizes"
    extended_genes_bed_df = create_extended_gene_bed(
        gtf_df,
        final_genes + final_tfs,  # if we want to include TF genes too
        window_size=ppConfig.window_size,
        chrom_sizes_path=chrom_sizes_path,
    )

    gene_bed_file = out_dir / f"genes_extended_{ppConfig.window_size//1000}kb.bed"
    extended_genes_bed_df.to_csv(gene_bed_file, sep="\t", header=False, index=False)
    logger.info(f"Created extended gene bed => {gene_bed_file}")

    if ppConfig.species.lower() == "human":
        # Fragment coordinates are in different format, chr1-start-end and not chr1:start-end, modified script accordingly
        data_atac.var["chr"] = [x.split("-")[0] for x in data_atac.var.index]
        data_atac.var["start"] = [
            int(x.split("-")[1]) for x in data_atac.var.index
        ]
        data_atac.var["end"] = [
            int(x.split("-")[2]) for x in data_atac.var.index
        ]
        data_atac.var["peak_name"] = data_atac.var.index
        all_peaks_bed = out_dir / "peaks_all.bed"
        data_atac.var[["chr", "start", "end", "peak_name"]].to_csv(
            all_peaks_bed, sep="\t", header=False, index=False
        )
    elif ppConfig.species.lower() == "mouse":

        data_atac.var["chr"] = [x.split(":")[0] for x in data_atac.var.index]
        data_atac.var["start"] = [
            int(x.split(":")[1].split("-")[0]) for x in data_atac.var.index
        ]
        data_atac.var["end"] = [int(x.split(":")[1].split("-")[1]) for x in data_atac.var.index]
        data_atac.var["peak_name"] = data_atac.var.index
        all_peaks_bed = out_dir / "peaks_all.bed"
        data_atac.var[["chr", "start", "end", "peak_name"]].to_csv(
            all_peaks_bed, sep="\t", header=False, index=False
        )
    else:
        raise ValueError("Species not recognized. Please use 'human' or 'mouse'.")

    intersected_bed = out_dir / "peaks_intersected.bed"
    run_bedtools_intersect(
        a_bed=all_peaks_bed, b_bed=gene_bed_file, out_bed=intersected_bed
    )
    peaks_intersected = pd.read_csv(intersected_bed, sep="\t", header=None)
    peaks_intersected.columns = ["chr", "start", "end", "peak_name"]
    windowed_set = set(peaks_intersected["peak_name"])

    # Subset data_atac to these peaks
    data_atac = data_atac[:, list(windowed_set)].copy()
    logger.info(f"After gene-window filtering => shape={data_atac.shape}")

    #
    rna_metacell, atac_metacell = create_metacells(
        data_rna,
        data_atac,
        grouping_key="leiden",
        resolution=ppConfig.leiden_resolution,
        batch_key=ppConfig.batch_key,
    )
    # Copy labels
    data_atac.obs["leiden"] = data_rna.obs["leiden"]

    data_atac = keep_promoters_and_select_hv_peaks(
        data_atac=data_atac,
        total_n_peaks=ppConfig.num_peaks,
        cluster_key="leiden",
        promoter_col=ppConfig.promoter_col,  # column in data_atac.var
    )

    logger.info(f"Final shape after combining promoters + HV => {data_atac.shape}")

    save_processed_datasets(data_rna, data_atac, out_dir)

    peaks_bed = out_dir / "peaks_selected.bed"
    data_atac.var[["chr", "start", "end", "peak_name"]].to_csv(
        peaks_bed, sep="\t", header=False, index=False
    )

    pwms_sub, key_to_tf = load_motif_database(motif_path, final_tfs)
    fasta_path = genome_dir / f"{ppConfig.genome_assembly}.fa"
    df_motif_scores = compute_motif_scores(
        bed_file=peaks_bed,
        fasta_file=fasta_path,
        pwms_sub=pwms_sub,
        key_to_tf=key_to_tf,
        n_peaks=data_atac.shape[1],
        window=500,
        threshold=ppConfig.motif_match_pvalue_threshold,
    )
    df_motif_scores = df_motif_scores[final_tfs]

    df_motif_scores.to_csv(out_dir / "motif_scores.tsv", sep="\t")

    atac_metacell = atac_metacell[:, data_atac.var_names].copy()
    tf_mask = rna_metacell.var["gene_type"] == "TF"
    rna_matrix = rna_metacell.X[:, tf_mask]  # shape=(n_meta, n_tfs)
    atac_matrix = atac_metacell.X  # shape=(n_meta, n_peaks)

    insilico_chipseq_act, insilico_chipseq_rep = compute_in_silico_chipseq(
        atac_matrix=atac_matrix,
        rna_matrix=rna_matrix,
        motif_scores=df_motif_scores,
        percentile=ppConfig.correlation_percentile,
        n_bg=ppConfig.n_bg_peaks_for_corr,
    )
    np.save(out_dir / "insilico_chipseq_act.npy", insilico_chipseq_act)
    np.save(out_dir / "insilico_chipseq_rep.npy", insilico_chipseq_rep)


    # distance is set to 0 if the peak midpoint is within gene-body or promoter (5kb upstream of TSS by default)
    # distance is -1 if peak-gene pairs on different chromosomes

    data_atac.var["index_int"] = range(data_atac.shape[1])
    selected_peak_indices = data_atac.var["index_int"].values

    # Subset GTF to final genes
    gene_info = gtf_df[gtf_df.feature == "gene"].drop_duplicates("gene_name")
    gene_info["gene"] = gene_info["gene_name"].values
    gene_info = gene_info.set_index("gene_name")
    gene_info = gene_info.loc[data_rna.var_names.intersection(gene_info.index)]

    gene_info["chr"] = gene_info["seqname"]  # rename col for consistency
    # Create gene_coordinates_intersect with necessary columns
    gene_info = gene_info[["chr", "start", "end", "strand", "gene"]].copy()
    gene_info.columns = ["chr_gene", "start", "end", "strand", "gene"]

    dist_matrix = compute_gene_peak_distance_matrix(
        data_rna=data_rna, data_atac=data_atac, gene_coordinates_intersect=gene_info
    )
    np.save(out_dir / "gene_peak_distance_raw.npy", dist_matrix)

    dist_matrix[dist_matrix < 0] = 1e8
    dist_matrix = np.exp(
        -1 * dist_matrix.astype(float) / ppConfig.peak_distance_scaling_factor
    )
    dist_matrix = np.where(dist_matrix < ppConfig.peak_distance_min_cutoff, 0, dist_matrix)
    np.save(out_dir / "gene_peak_distance_exp.npy", dist_matrix)

    logger.info("=== Preprocessing pipeline completed successfully ===")
    return
