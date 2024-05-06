library(rhdf5)
library(dplyr)
library(doParallel)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
nCores <- 16

# Read data
print('Open object')
indata <- H5Fopen(path_data)

# RNA
rna_data <- as.data.frame(indata$mod$rna$X)
colnames(rna_data) <- indata$obs$`_index`
rownames(rna_data) <- indata$mod$rna$var$`_index`

# ATAC
atac_data <- Matrix::sparseMatrix(
    i=indata$mod$atac$layers$counts$indices,
    p=indata$mod$atac$layers$counts$indptr,
    x=as.numeric(indata$mod$atac$layers$counts$data),
    index1 = FALSE
)
colnames(atac_data) <- indata$obs$`_index`
rownames(atac_data) <- indata$mod$atac$var$`_index`

# Dim reduction
#lsi <- indata$obsm$X_spectral
#colnames(lsi) <- colnames(rna_data)
#rownames(lsi) <- paste0("SPL_", 1:nrow(lsi))

# Normalize ATAC data
atac_data <- as.matrix(FigR::centerCounts(atac_data, chunkSize = 100000))
colnames(atac_data) <- as.character(colnames(atac_data))
rownames(atac_data) <- as.character(rownames(atac_data))

# Find KNN based on LSI
#cellknn <- FNN::get.knn(t(lsi), k = k)$nn.index
#rownames(cellknn) <- colnames(rna_data)

# KNN impute
#rna_data <- as.matrix(FigR::smoothScoresNN(NNmat = cellknn, mat = rna_data, nCores = nCores))
#atac_data <- as.matrix(FigR::smoothScoresNN(NNmat = cellknn, mat = atac_data, nCores = nCores))

# Write
#h5write(rna_data, name="mod/rna/X", file=indata)
h5write(atac_data, name="mod/atac/X", file=indata)

# Close
h5closeAll()
