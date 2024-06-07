library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
granges_hg <- args[8]
granges_mm <- args[9]
rna_network_path <- args[10]
atac_network_path <- args[11]
hummus_object_f <- args[12]
grn_number_edges <- as.numeric(args[13])
tf_layer_method <- args[14]
if (tf_layer_method == "None"){
    tf_layer_method = NULL
}

# Set genome
if (organism == 'hg38'){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    annot <- read.csv(granges_hg)
    genome_annot = get_genome_annotations(
        ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
} else if (organism == 'mm10'){
    annot <- read.csv(granges_mm)
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    genome_annot = get_genome_annotations(
        ensdb_annotations = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
}

# Read peaks and genes
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$mod$rna$obs$`_index`
rownames(rna_X) <- stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_')

#ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$mod$atac$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
h5closeAll()


# Create a HuMMuS object through a Seurat object
seurat_object <- SeuratObject::CreateSeuratObject(rna_X)
rm(rna_X)
seurat_object[['peaks']] <- Signac::CreateChromatinAssay(
  counts = atac_X,
  sep = c("-", "-"))
rm(atac_X)
hummus <- Initiate_Hummus_Object(seurat_object)
rm(seurat_object)
hummus
hummus[['RNA']]
hummus[['peaks']]


# Filter annot by seen genes and add it to peaks matrix
annot <- annot[annot$gene_name %in% intersect(colnames(hummus[['RNA']]), annot$gene_name)]


genome_annot <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
Signac::Annotation(hummus@assays$peaks) <- genome_annot
rm(annot)

rna_network <- read.csv2(rna_network_path, sep = ",")
rna_network <- rna_network[1:grn_number_edges,]
# Add external networks
hummus <- add_network(
    hummus,
    rna_network,
    multiplex_name = "RNA",
    network_name = "GRN_Boost2",
    weighted = TRUE)


atac_network <- read.csv2(atac_network_path, sep = ",")
#atac_network[,'peak1'] <- stringr::str_replace_all(atac_network[,'peak1'], '-', '_')
#atac_network[,'peak2'] <- stringr::str_replace_all(atac_network[,'peak2'], '-', '_')

print(head(atac_network))
hummus <- add_network(
    hummus,
    atac_network,
    multiplex_name = "peaks",
    network_name = "AtacNet",
    weighted = TRUE)

# Connect TF to peaks !!!!!TODO: to move in tf2p
hummus@motifs_db <- get_tf2motifs() # TF motif DB # by default human motifs
hummus <- compute_tf_network(hummus, method = tf_layer_method)


saveRDS(hummus, hummus_object_f)
