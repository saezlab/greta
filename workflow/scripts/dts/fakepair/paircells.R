library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
options("optmatch_max_problem_size" = Inf)
optmatch::setMaxProblemSize(size = Inf)


# Parse args
args <- commandArgs(trailingOnly = F)
path_cca <- args[6]
path_ctypes <- args[7]
path_barMap_out <- args[8]


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
#euc.dist <- function(x1, x2) sqrt(sum((x1 - x2) ^ 2))
#pairing$dist <- apply(pairing, 1, function(x) { euc.dist(ATAC_PCs[x[1],1:ncol(ATAC_PCs)],RNA_PCs[x[2],1:ncol(RNA_PCs)])})
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairing <- pairing[!duplicated(pairing$ATAC),]
#atac_pairing <- pairing[!duplicated(pairing$ATAC),]
#rna_pairing <- pairing[!duplicated(pairing$RNA),]
#pairing <- merge(atac_pairing, rna_pairing)

# Merge ctype info
ctypes <- read.csv(path_ctypes)
pairing <- merge(pairing, ctypes, by.x='ATAC', by.y='barcode')
pairing['batch'] <- 'smpl'
pairing <- pairing[, c('ATAC', 'RNA', 'batch', 'celltype', 'dist')]

# Write
write.csv(pairing, path_barMap_out, row.names = FALSE)