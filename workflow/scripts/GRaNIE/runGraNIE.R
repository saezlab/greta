library(readr)
library(GRaNIE)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_all_peaks <- args[8]
path_connections <- args[9]

# Read genome
if (organism == 'human'){
  genomeAssembly = "hg38"
} else {
  stop("Not implemented yet")
}

source("workflow/scripts/GRaNIE/functions.R")

runGRaNIE (dir_output = "output_GRaNIE", 
                      datasetName = "undescribed",
                      file_peaks, file_rna, file_metadata,
                      TFBS_source = "custom",
                      TFBS_folder = NULL,
                      TFBS_JASPAR_useSpecificTaxGroup = NULL,
                      genomeAssembly = "hg38",
                      normalization_peaks = "DESeq2_sizeFactors", 
                      idColumn_peaks = "peakID",
                      normalization_rna = "limma_quantile", 
                      idColumn_RNA =  "ENSEMBL",
                      includeSexChr = FALSE,
                      minCV = 0,
                      minNormalizedMean_peaks = 5,
                      minNormalizedMean_RNA = 1,
                      minSizePeaks = 5,
                      corMethod = "pearson",
                      promoterRange = 250000, 
                      useGCCorrection = FALSE,
                      TF_peak.fdr.threshold = 0.2,
                      peak_gene.fdr.threshold = 0.1,
                      runTFClassification = FALSE,
                      runNetworkAnalyses = FALSE, 
                      nCores = 4,
                      forceRerun = TRUE
)

# TODO
# Save
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(path_all_peaks))
write.csv(x = conns, file = file.path(path_connections))


