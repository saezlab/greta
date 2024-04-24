library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
tf_list_file <- args[8]
tf_motif_file <- args[9]

# Read MuData
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- Matrix::sparseMatrix(
  i = as.vector(indata$mod$rna$X[['indices']][] + 1),
  p = as.vector(indata$mod$rna$X[['indptr']][]),
  x = as.vector(indata$mod$rna$X[['data']][])
)

print(rna_X)
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_')
h5closeAll()

print(Matrix::t(rna_X[1:5,1:5]))
# Create HuMMuS object
seurat <- SeuratObject::CreateSeuratObject(
    assay = 'RNA',
    counts = rna_X[,1:2])

seurat[['RNA']]

print(seurat)
hummus <- HuMMuS::Initiate_Hummus_Object(seurat)
rm(seurat)

hummus[['RNA']]

# Set genome
if (organism == 'hg38'){
    hummus@motifs_db <- get_tf2motifs("human")
} else if (organism == 'mm10'){
    hummus@motifs_db <- get_tf2motifs("mouse")
}

# store expressed TFs list
get_tfs(hummus = hummus,
        assay = 'RNA', # Assay to check expression
        store_tfs = TRUE,
        output_file = tf_list_file,
        verbose = 1)

hummus@motifs_db

# store TF motifs
saveRDS(hummus, tf_motif_file)
