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
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairing <- pairing[!duplicated(pairing$ATAC),]

# Merge ctype info
ctypes <- read.csv(path_ctypes)
pairing <- merge(pairing, ctypes, by.x='RNA', by.y='X')
pairing['batch'] <- 'smpl'
pairing <- pairing[, c('ATAC', 'RNA', 'batch', 'celltype', 'dist')]
rownames(pairing) <- pairing$ATAC

# Write
write.csv(pairing, path_barMap_out, row.names = FALSE)
