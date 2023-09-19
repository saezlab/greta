
GRaNIE_singularity_path = "envs/GRaNIE.sif"

rule preprocess_GRaNIE:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    output:
        counts_rna  = "resources/{dataset}/{trajectory}/GRaNIE/counts_rna.tsv.gz",
        counts_atac = "resources/{dataset}/{trajectory}/GRaNIE/counts_atac.tsv.gz"
        metadata    = "resources/{dataset}/{trajectory}/GRaNIE/metadata.tsv.gz"
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.preprocess_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    params:
        pseudubulk = lambda w: config[w.dataset]['trajectories'][w.trajectory]["GRaNIE"]["preprocessing_pseudobulk"]
    shell:
        """
        Rscript scripts/GRaNIE/1.preprocess.R \
         --input {input} \
         --output {output} \
         --pseudobulk {params.pseudobulk}
        """

def getConfigValue(dataset, trajectory, par):
    return config[dataset]['trajectories'][trajectory]["GRaNIE"][par]


rule run_GRaNIE:
    input:
        data=rules.preprocess_GRaNIE.output
    output:
        data = "resources/{dataset}/{trajectory}/GRaNIE/connections_filtered.tsv.gz",
        GRN  = "resources/{dataset}/{trajectory}/GRaNIE/GRN.qs"
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.run_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    params:
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
        TF_peak.fdr.threshold   = lambda w: getConfigValue(w.dataset, w.trajectory, "TF_peak.fdr.threshold"),
        peak_gene.fdr.threshold = lambda w: getConfigValue(w.dataset, w.trajectory, "peak_gene.fdr.threshold"),
        runTFClassification     = lambda w: getConfigValue(w.dataset, w.trajectory, "runTFClassification"),
        runNetworkAnalyses      = lambda w: getConfigValue(w.dataset, w.trajectory, "runNetworkAnalyses")
    shell:
        """
        Rscript scripts/GRaNIE/2.runGRaNIE.R \
         --input {input} \
         --output {output} \
         --TBGS_source {params.TFBS_source} \
         --normalization_peaks {params.normalization_peaks} \
         --normalization_rna {params.normalization_rna} \
         --includeSexChr {params.includeSexChr} \
         --minCV {params.minCV} \
         --minNormalizedMean_peak {params.minNormalizedMean_peak} \
         --minNormalizedMean_RNA {params.minNormalizedMean_RNA} \
         --minSizePeaks {params.minSizePeaks} \
         --corMethod {params.corMethod} \
         --promoterRange {params.promoterRange} \
         --TF_peak.fdr.threshold {params.TF_peak.fdr.threshold} \
         --peak_gene.fdr.threshold {params.peak_gene.fdr.threshold} \
         --runTFClassification {params.runTFClassification} \
         --runNetworkAnalyses {params.runNetworkAnalyses}
        """


rule postprocess_GRaNIE:
    input:
        data=rules.run_GRaNIE.output
    output:
        TF_gene      = "resources/{dataset}/{trajectory}/GRaNIE/grn.csv",
        TF_peak_gene = "resources/{dataset}/{trajectory}/GRaNIE/tri.csv"
    singularity: GRaNIE_singularity_path
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.postprocess_GRaNIE.txt"
    params:
    shell:
        """
        Rscript scripts/GRaNIE/3.postprocess.R \
         --input {input} \
         --output {output}
        """
