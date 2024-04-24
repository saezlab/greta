library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
granges_hg <- args[8]
granges_mm <- args[9]
extend <- as.numeric(args[10])
path_out <- args[11]
rna_network_path <- args[12]
atac_network_path <- args[13]
multilayer_f <- args[14]

# Set genome
if (organism == 'hg38'){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    annot <- read.csv(granges_hg)
} else if (organism == 'mm10'){
    annot <- read.csv(granges_mm)
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
}
print(annot)

# Read peaks and genes
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- as.matrix(indata$mod$rna$X)
#rna_X <- Matrix::sparseMatrix(
#  i = as.vector(indata$mod$rna$X[['indices']][] + 1),
#  p = as.vector(indata$mod$rna$X[['indptr']][]),
#  x = as.vector(indata$mod$rna$X[['data']][])
#)

#print(rna_X)
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_')

#ATAC
atac_X <- as.matrix(indata$mod$atac$X)
#atac_X <- Matrix::sparseMatrix(
#  i = as.vector(indata$mod$atac$X[['indices']][] + 1),
#  p = as.vector(indata$mod$atac$X[['indptr']][]),
#  x = as.vector(indata$mod$atac$X[['data']][])
#)

#print(atac_X)
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- stringr::str_replace_all(indata$mod$atac$var$`_index`, '-', '_')


print(Matrix::t(rna_X[1:5,1:5]))
# Create HuMMuS object
h5closeAll()


# Create a HuMMuS object through a Seurat object
seurat_object <- SeuratObject::CreateSeuratObject(rna_X)
rm(rna_X)
seurat_object[['peaks']] <- Signac::CreateChromatinAssay(counts = atac_X, sep=c('_', '_'))
rm(atac_X)
hummus = Initiate_Hummus_Object(seurat_object)
rm(seurat_object)
hummus
hummus[['RNA']]
hummus[['peaks']]


# Filter annot by seen genes and add it to peaks matrix
annot <- annot[annot$gene_name %in% intersect(colnames(hummus[['RNA']]), annot$gene_name)]


genome_annot = get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
Signac::Annotation(hummus@assays$peaks) = genome_annot
rm(annot)

rna_network = read.csv2(rna_network_path, sep=",")
print(head(rna_network))
# Add external networks
hummus <- add_network(
    hummus, 
    rna_network,
    multiplex_name = "RNA",
    network_name = "GRN_Boost2",
    weighted=TRUE)


atac_network = read.csv2(atac_network_path, sep="\t")
print(head(atac_network))
hummus <- add_network(
    hummus, 
    atac_network,
    multiplex_name = "peaks",
    network_name = "AtacNet",
    weighted=TRUE)

# Connect TF to peaks !!!!!TODO: to move in tf2p
hummus@motifs_db <- get_tf2motifs() # TF motif DB # by default human motifs
hummus <- compute_tf_network(hummus, method=NULL)

hummus <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = NULL, # use to filter TF on only expressed TFs,
                                     # if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_multiplex_name = "TF",
              genome = genome,
              )

hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      upstream = 500,
                      downstream = 500,                      
                      )

save_multilayer(hummus,
                folder_name = multilayer_f)

# Run HuMMuS enhancers search
enhancers <- define_enhancers(
  hummus,
#  gene_list = list("ATF2"),
  multilayer_f = multilayer_f,
  njobs = 1
  )

#get only gene, peaks and score columns
enhancers <- enhancers[, c("gene", "peak", "score")]
colnames(enhancers) <- c("gene", "cre", "score")

# Write
write.csv(
    x = enhancers,
    file = path_out,
    row.names = FALSE,
    quote = FALSE)
