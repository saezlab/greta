
GRaNIE_singularity_path = "envs/GRaNIE.sif"

def getOutputFiles(dataset, trajectory):
    """"Helper function to retrieve the GRaNIE output file names based on config parameters"""

    filenames =
    expand("resources/{dataset}/{trajectory}/GRaNIE/{modality}_minClusters{minClusters}_{countAggr}.tsv.gz",
        dataset = dataset, trajectory = trajectory,
        modality = ["rna.pseudobulkFromClusters", "atac.pseudobulkFromClusters", "metadata"],
        minClusters = config[dataset]['trajectories'][trajectory]["GRaNIE"]["preprocessing_minClusters"],
        countAggr   = config[dataset]['trajectories'][trajectory]["GRaNIE"]["preprocessing_countAggregation"]
    )

    return filenames

rule preprocess_GRaNIE:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    output:
        lambda w: getOutputFiles(w.dataset, w.trajectory)
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.preprocess_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    params:
        #pseudubulk                       = lambda w: getConfigValue(w.dataset, w.trajectory,  "preprocessing_pseudobulk"),
        preprocessing_countAggregation   = lambda w: getConfigValue(w.dataset, w.trajectory, "preprocessing_countAggregation"),
        preprocessing_minCellsPerCluster = lambda w: getConfigValue(w.dataset, w.trajectory, "preprocessing_minCellsPerCluster"),
        preprocessing_minClusters        = lambda w: getConfigValue(w.dataset, w.trajectory, "preprocessing_minClusters"),
        organism   = lambda w: config[w.dataset]["organism"]
    shell:
        """
        Rscript scripts/GRaNIE/1.preprocess.R \
         --input {input} \
         --output {output} \
         --organism {params.organism} \
         --preprocessing_countAggregation {params.preprocessing_countAggregation} \
         --preprocessing_minCellsPerCluster {params.preprocessing_minCellsPerCluster} \
         --preprocessing_minClusters  {params.preprocessing_minClusters }
        """

def getConfigValue(dataset, trajectory, par):
    """"Helper function to retrieve the GRaNIE config values"""
    return config[dataset]['trajectories'][trajectory]["GRaNIE"][par]


rule run_GRaNIE:
    input:
        data=rules.preprocess_GRaNIE.output
    output:
        GRN          = "resources/{dataset}/{trajectory}/GRaNIE/GRN.qs"
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.run_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    threads: 8
    params:
        organism                = lambda w: config[w.dataset]["organism"],
        name                    = lambda w: w.dataset + "." + w.trajectory,
        TFBS_source             = lambda w: getConfigValue(w.dataset, w.trajectory, "TFBS_source"),
        normalization_peaks     = lambda w: getConfigValue(w.dataset, w.trajectory, "normalization_peaks"),
        normalization_rna       = lambda w: getConfigValue(w.dataset, w.trajectory, "normalization_rna"),
        includeSexChr           = lambda w: getConfigValue(w.dataset, w.trajectory, "includeSexChr"),
        minCV                   = lambda w: getConfigValue(w.dataset, w.trajectory, "minCV"),
        minNormalizedMean_peak  = lambda w: getConfigValue(w.dataset, w.trajectory, "minNormalizedMean_peak"),
        minNormalizedMean_RNA   = lambda w: getConfigValue(w.dataset, w.trajectory, "minNormalizedMean_RNA"),
        minSizePeaks            = lambda w: getConfigValue(w.dataset, w.trajectory, "minSizePeaks"),
        corMethod               = lambda w: getConfigValue(w.dataset, w.trajectory, "corMethod"),
        promoterRange           = lambda w: getConfigValue(w.dataset, w.trajectory, "promoterRange"),
        TF_peak_fdr             = lambda w: getConfigValue(w.dataset, w.trajectory, "TF_peak.fdr.threshold"),
        peak_gene_fdr           = lambda w: getConfigValue(w.dataset, w.trajectory, "peak_gene.fdr.threshold"),
        runTFClassification     = lambda w: getConfigValue(w.dataset, w.trajectory, "runTFClassification"),
        runNetworkAnalyses      = lambda w: getConfigValue(w.dataset, w.trajectory, "runNetworkAnalyses")
    shell:
        """
        Rscript scripts/GRaNIE/2.runGRaNIE.R \
         --input {input} \
         --output {output.data} \
         --threads {threads} \
         --organism {params.organism} \
         --name {params.name} \
         --TBFS_source {params.TFBS_source} \
         --normalization_peaks {params.normalization_peaks} \
         --normalization_rna {params.normalization_rna} \
         --includeSexChr {params.includeSexChr} \
         --minCV {params.minCV} \
         --minNormalizedMean_peak {params.minNormalizedMean_peak} \
         --minNormalizedMean_RNA {params.minNormalizedMean_RNA} \
         --minSizePeaks {params.minSizePeaks} \
         --corMethod {params.corMethod} \
         --promoterRange {params.promoterRange} \
         --TF_peak_fdr {params.TF_peak_fdr} \
         --peak_gene_fdr {params.peak_gene_fdr} \
         --runTFClassification {params.runTFClassification} \
         --runNetworkAnalyses {params.runNetworkAnalyses}
        """

rule postprocess_GRaNIE:
    input:
        GRN = rules.run_GRaNIE.output
    output:
        TF_gene      = "resources/{dataset}/{trajectory}/GRaNIE/grn.csv",
        TF_peak_gene = "resources/{dataset}/{trajectory}/GRaNIE/tri.csv"
    singularity: GRaNIE_singularity_path
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.postprocess_GRaNIE.txt"
    params:
        ranking_TF_gene      = lambda w: getConfigValue(w.dataset, w.trajectory, "ranking_TF_gene"),
        ranking_TF_peak_gene = lambda w: getConfigValue(w.dataset, w.trajectory, "ranking_TF_peak_gene")
    shell:
        """
        Rscript scripts/GRaNIE/3.postprocess.R \
         --input {input.GRN} \
         --output {output} \
         --ranking_TF_gene {params.ranking_TF_gene} \
         --ranking_TF_peak_gene {params.ranking_TF_peak_gene}
        """
