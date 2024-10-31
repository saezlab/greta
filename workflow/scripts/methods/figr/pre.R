library(rhdf5)
library(dplyr)
library(doParallel)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
nCores <- as.numeric(args[7])

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

# Normalize ATAC data
atac_data <- as.matrix(FigR::centerCounts(atac_data, chunkSize = 100000))
colnames(atac_data) <- as.character(colnames(atac_data))
rownames(atac_data) <- as.character(rownames(atac_data))

# Write
h5write(atac_data, name="mod/atac/X", file=indata)

# Close
h5closeAll()
