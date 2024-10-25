library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
options("optmatch_max_problem_size" = Inf)
optmatch::setMaxProblemSize(size = Inf)


# Parse args
args <- commandArgs(trailingOnly = F)
path_cca <- args[8]
path_barMap_out <- args[10]


# Load Data
CCA_PCs <- readRDS(path_cca)
isATAC <- grepl("^smpl_",rownames(CCA_PCs))
ATAC_PCs <- CCA_PCs[isATAC,]
RNA_PCs <- CCA_PCs[!isATAC,]

# Pair with FigR
pairing <- pairCells(
    ATAC = ATAC_PCs,
    RNA = RNA_PCs,
    keepUnique = TRUE
)

# Filter paired object
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairing <- pairing[!duplicated(pairing$ATAC),]

# Write
write.csv(pairing, path_barMap_out, row.names = FALSE)
