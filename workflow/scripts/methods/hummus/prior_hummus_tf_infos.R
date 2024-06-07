library(tidyverse)
library(rhdf5)
library(HuMMuS)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
tf_list_file <- args[8]

# Read MuData
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')

# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$mod$rna$obs$`_index`
rownames(rna_X) <- stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_')
h5closeAll()


# Create HuMMuS object
seurat <- SeuratObject::CreateSeuratObject(
    assay = "RNA",
    counts = rna_X[, 1:2])

hummus <- HuMMuS::Initiate_Hummus_Object(seurat)
rm(seurat)


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
