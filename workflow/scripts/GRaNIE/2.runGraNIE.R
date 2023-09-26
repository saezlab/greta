library(argparse)
library(readr)
library(GRaNIE)

parser <- ArgumentParser(description= 'Run GRaNIE')
parser$add_argument('--input', '-i', help= 'Input files: RNA, ATAC, metadata', nargs = 3,  required= TRUE)
parser$add_argument('--output', '-o', help= 'List of 3 filenames: 1. GRN object, 2. TF_gene table, 3. TF_peak_gene table', nargs= 3, required= TRUE)
parser$add_argument('--threads', help= 'Threads', nargs= 1, required= TRUE)
parser$add_argument('--organism', help= 'Organism', nargs= 1, required= TRUE)
parser$add_argument('--name', help= 'Dataset name', nargs= 1, required= TRUE) 
parser$add_argument('--TBFS_source', help= 'TBFS_source', nargs= 1, required= TRUE) 
parser$add_argument('--normalization_peaks', help= 'normalization_peaks', nargs= 1, required= TRUE) 
parser$add_argument('--normalization_rna', help= 'normalization_rna', nargs= 1, required= TRUE) 
parser$add_argument('--includeSexChr', help= 'includeSexChr', nargs= 1, required= TRUE) 
parser$add_argument('--minCV', help= 'minCV', nargs= 1, required= TRUE) 
parser$add_argument('--minNormalizedMean_peak', help= 'minNormalizedMean_peak', nargs= 1, required= TRUE) 
parser$add_argument('--minNormalizedMean_RNA', help= 'minNormalizedMean_RNA', nargs= 1, required= TRUE) 
parser$add_argument('--minSizePeaks', help= 'minSizePeaks', nargs= 1, required= TRUE) 
parser$add_argument('--corMethod', help= 'corMethod', nargs= 1, required= TRUE) 
parser$add_argument('--promoterRange', help= 'promoterRange', nargs= 1, required= TRUE) 
parser$add_argument('--TF_peak_fdr', help= 'TF_peak.fdr.threshold', nargs= 1, required= TRUE) 
parser$add_argument('--peak_gene_fdr', help= 'peak_gene.fdr.threshold', nargs= 1, required= TRUE) 
parser$add_argument('--runTFClassification', help= 'runTFClassification', nargs= 1, required= TRUE) 
parser$add_argument('--runNetworkAnalyses', help= 'runNetworkAnalyses', nargs= 1, required= TRUE) 


xargs <- parser$parse_args()


# Parse args
outdir <-dirname(xargs$output[1])

# Read genome
if (xargs$organism == 'human'){
  genomeAssembly = "hg38"
} else {
  stop("Not implemented yet")
}

source("workflow/scripts/GRaNIE/functions.R")

GRN = runGRaNIE (dir_output = outdir, 
                      datasetName = xargs$name,
                      file_rna = args$input[1], file_peaks = args$input[2], file_metadata = xargs$input[3],
                      TFBS_source = xargs$TBFS_source,
                      TFBS_folder = NULL,
                      TFBS_JASPAR_useSpecificTaxGroup = NULL,
                      genomeAssembly = genomeAssembly,
                      normalization_peaks = xargs$normalization_peaks, 
                      idColumn_peaks = "peakID",
                      normalization_rna = xargs$normalization_rna, 
                      idColumn_RNA =  "ENSEMBL",
                      includeSexChr = xargs$includeSexChr,
                      minCV = xargs$minCV,
                      minNormalizedMean_peaks = xargs$minNormalizedMean_peak,
                      minNormalizedMean_RNA = xargs$minNormalizedMean_RNA,
                      minSizePeaks = xargs$minSizePeaks,
                      corMethod = xargs$corMethod,
                      promoterRange = xargs$promoterRange, 
                      useGCCorrection = FALSE,
                      TF_peak.fdr.threshold = xargs$TF_peak_fdr,
                      peak_gene.fdr.threshold = xargs$peak_gene_fdr,
                      runTFClassification = xargs$runTFClassification,
                      runNetworkAnalyses = xargs$runNetworkAnalyses, 
                      nCores = xargs$threads,
                      forceRerun = TRUE
)

