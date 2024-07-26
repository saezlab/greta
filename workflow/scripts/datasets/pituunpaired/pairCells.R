library(doParallel)
library(SingleCellExperiment)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)

options("optmatch_max_problem_size" = Inf)

# Parse args
args <- commandArgs(trailingOnly = F)
path_exprMat <- args[1]
path_atacSE <- args[2]
path_cca <- args[3]
path_barMap_out <- args[4]



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

pairing <- pairCells(ATAC = ATAC_PCs,
            RNA = RNA_PCs,
            keepUnique = TRUE)


# Filter paired object
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairingClean <- pairing[!duplicated(pairing$ATAC),]
pairingClean <- pairingClean[c("ATAC","RNA")]
pairingClean <- data.frame(lapply(pairingClean, function(x) {gsub("_.*", "", x)}))
write.csv(pairingClean, path_barMap_out, row.names = FALSE)


# Optional: Filter data based on pairing
# Adapt colnames
#colnames(ATAC.se) <- paste0(colnames(ATAC.se), "_2")
#colnames(RNAmat) <- paste0(colnames(RNAmat), "_1")

# Filter data based on pairing
#ATAC.se.paired <- ATAC.se[,pairingClean$ATAC]
#RNAmat.paired <- RNAmat[,pairingClean$RNA]