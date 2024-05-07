library(rhdf5)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]

# Read data
print('Open object')
indata <- H5Fopen(path_data)

# RNA
rna_data <- indata$mod$rna$X
colnames(rna_data) <- indata$obs$`_index`
rownames(rna_data) <- indata$mod$rna$var$`_index`

### ATAC
atac_data <- indata$mod$atac$X
colnames(atac_data) <- indata$obs$`_index`
rownames(atac_data) <- indata$mod$atac$var$`_index`

# Normalize data
norm_data <- function(data, norm){
    if (norm == 'deseq2'){
        data <- DESeq2::DESeqDataSetFromMatrix(
            countData = data,
            colData = data.frame(sampleID = colnames(data)),
            design = stats::as.formula(" ~ 1")
        )
        data <- DESeq2::estimateSizeFactors(data)
        data <- DESeq2::counts(data, normalized = TRUE)
    }
    if (norm == 'limma'){
        data <- limma::normalizeBetweenArrays(
            data,
            method = 'quantile'
        )
    }
    return(data)
}
# Add pseudocounts for sparsity and normalize
rna_data <- norm_data(rna_data, 'limma')
atac_data <- norm_data(atac_data, 'deseq2')

# Write
h5write(rna_data, name="mod/rna/X", file=indata)
h5write(atac_data, name="mod/atac/X", file=indata)

# Close
h5closeAll()
