library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(Matrix)
library(chromVAR)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
ext <- as.numeric(args[8])
ncres <- as.numeric(args[9])
path_out <- args[10]

nCores <- 32

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`

# ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
h5closeAll()

# Transform atac to sme Object
peaks <- strsplit(rownames(atac_X), "-")
peak_ranges <- GenomicRanges::GRanges(
    seqnames = sapply(peaks, "[[", 1),
    ranges = IRanges::IRanges(start = as.numeric(sapply(peaks, "[[", 2)), end = as.numeric(sapply(peaks, "[[", 3)))
)
atac_X <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=atac_X), rowRanges = peak_ranges)

# Run p2g
cisCorr <- FigR::runGenePeakcorr(
    ATAC.se = atac_X,
    RNAmat = rna_X,
    genome = organism,
    nCores = 6,
    p.cut = NULL,
    n_bg = 100,
    normalizeATACmat = FALSE,
    windowPadSize = round(ext / 2),
)

# Process
p2g <- cisCorr %>%
    filter(pvalZ <= 0.05) %>%
    mutate(cre=stringr::str_replace(PeakRanges, ':', '-')) %>%
    rename(gene=Gene, score=rObs, pval=pvalZ) %>%
    select(cre, gene, score, pval)

# Identify dorcs and keep only their g-p pairs
dorcs <- p2g %>% summarize(counts=n(), .by=gene) %>% filter(counts >= ncres) %>% pull(gene)
p2g <- p2g %>% filter(gene %in% dorcs)

# Write
write.csv(x = p2g, file = path_out, row.names=FALSE)
