library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)

options("optmatch_max_problem_size" = Inf)


# Parse args
args <- commandArgs(trailingOnly = F)
path_exprMat <- args[6]
path_atacSE <- args[7]
path_cca <- args[8]
path_celltypes <- args[9]
path_barMap_out <- args[10]



# Load Data
RNAmat <- readRDS(path_exprMat)
ATAC.se <- readRDS(path_atacSE)
CCA_PCs <- readRDS(path_cca)


# Identify ATAC and RNA cells
isATAC <- grepl("_2",rownames(CCA_PCs))
table(isATAC) # ATAC vs RNA


# Pair cells with FigR 
ATAC_PCs <- CCA_PCs[isATAC,]
RNA_PCs <- CCA_PCs[!isATAC,]

library('optmatch')
setMaxProblemSize(size = Inf)

pairing <- pairCells(ATAC = ATAC_PCs,
            RNA = RNA_PCs,
            keepUnique = TRUE)


# Filter paired object
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairingClean <- pairing[!duplicated(pairing$ATAC),]
pairingClean <- pairingClean[c("ATAC","RNA")]
pairingClean <- data.frame(lapply(pairingClean, function(x) {gsub("_.*", "", x)}))

# Create annotation
celltypes <- read.csv(path_celltypes)

celltypes <- celltypes[celltypes$X %in% pairingClean$RNA,]
celltypes <- celltypes[order(celltypes$X), ]
pairingClean <- pairingClean[order(pairingClean$RNA),] 
pairingClean$celltype <- celltypes$celltype

# Format and rename
rownames(pairingClean) <- pairingClean$ATAC
pairingClean$batch <- c("smpl")
colnames(pairingClean) <- c("barcodes", "RNA", "celltype", "batch")


write.csv(pairingClean, path_barMap_out, row.names = FALSE)

