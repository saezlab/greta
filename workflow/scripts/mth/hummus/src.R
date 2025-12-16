library(reticulate)
use_condaenv("env", conda = "/opt/conda/bin/conda", required = TRUE)
hummuspy <- import("hummuspy")

library(HuMMuS)
library(tidyverse)
library(rhdf5)


grn_number_edges <- 50000

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
rna_network_path <- args[7]
atac_network_path <- args[8]
tmp_dir <- args[9]
n_cores <- as.numeric(args[10])
output_file <- args[11]
binding_regions_file <- args[12]

print("Open CIRCE layer")
atac_network <- read.csv(atac_network_path, sep = ",")
atac_network <- atac_network[, c("Peak1", "Peak2", "score")]
atac_network <- atac_network[atac_network[, "score"] > 0, ]


# Load genome information
print('Loading genome information...')
genome_annot <- get_genome_annotations(
  ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
species <- "human"
tf_layer_method <- NULL

# Load the h5mu file as a Seurat object
print('Loading h5mu file...')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
##  RNA
rna_X <- indata$mod$rna$X

colnames(rna_X) <- toupper(indata$mod$rna$obs$`_index`)
rownames(rna_X) <- toupper(stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_'))
rna_X = rna_X[, 1:2]
## ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$mod$atac$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
atac_X = atac_X[,1:2]
h5closeAll()

# Create a HuMMuS object through a Seurat object
## RNA
seurat_object <- SeuratObject::CreateSeuratObject(rna_X)
rm(rna_X)
## ATAC
print('Creating Seurat object with ATAC data...')
seurat_object@assays$peaks <- Signac::CreateChromatinAssay(
  counts = atac_X,
  sep = c(":", "-"))
rm(atac_X)
## Initiate HuMMuS object
hummus <- Initiate_Hummus_Object(seurat_object)
rm(seurat_object)

# Add genome annotations
Signac::Annotation(hummus@assays$peaks) <- genome_annot
write.csv(
  as.data.frame(Signac::Annotation(hummus@assays$peaks)),
  file = file.path(tmp_dir, "genome_annotations.csv"),
  row.names = FALSE,
  quote = FALSE)

rm(genome_annot)

# Add TF motifs
hummus@motifs_db <- get_tf2motifs(species) # TF motif DB # by default human motifs
hummus <- compute_tf_network(hummus, method = tf_layer_method)
print(hummus@assays)


print('Computing the bipartite network between peaks and genes...')
# Compute the bipartite network between peaks and genes
hummus <- bipartite_peaks2genes(
                      hummus_object = hummus,
                      gene_assay = "RNA",
                      peak_assay = "peaks",
                      store_network = FALSE,
                      )

print('Computing the bipartite network between TFs and peaks...')
# Compute the bipartite network between TFs and peaks
hummus <- bipartite_tfs2peaks(
              hummus_object = hummus,
              tf_expr_assay = "RNA", # use to filter TF on only expressed TFs,
                                     # if NULL, all TFs with motifs are used
              peak_assay = "peaks",
              tf_multiplex_name = "TF",
              genome = genome,
              store_network = FALSE,
              )

print("Open GRNBoost2 layer")
rna_network <- read.csv(rna_network_path, sep = ",")
rna_network <- rna_network[, c("TF", "target", "importance")]
rna_network <- rna_network[1:min(grn_number_edges, nrow(rna_network)), ]
# Add external networks
hummus <- add_network(
    hummus,
    rna_network,
    multiplex_name = "RNA",
    network_name = "GRNBoost2",
    weighted = TRUE)


print("Open CIRCE layer")
atac_network <- read.csv(atac_network_path, sep = ",")
atac_network <- atac_network[, c("Peak1", "Peak2", "score")]
atac_network <- atac_network[atac_network[, "score"] > 0, ]
#atac_network[,'Peak1'] <- stringr::str_replace_all(
#  atac_network[,'Peak1'], '-', '_')
atac_network[,'Peak1'] <- stringr::str_replace_all(
  atac_network[,'Peak1'], ':', '-')
#atac_network[,'Peak2'] <- stringr::str_replace_all(
#  atac_network[,'Peak2'], '-', '_')
atac_network[,'Peak2'] <- stringr::str_replace_all(
  atac_network[,'Peak2'], ':', '-')

# Add external networks
hummus <- add_network(
    hummus,
    atac_network,
    multiplex_name = "peaks",
    network_name = "CIRCE",
    weighted = TRUE)

# Save the HuMMuS object
#save rds in tmp dir
hummus_rds_path <- file.path(tmp_dir, "hummus.rds")
print(paste0('Saving the HuMMuS object to ', hummus_rds_path))
saveRDS(hummus, file = hummus_rds_path)
 
# Run HuMMuS exploration
multilayer_f <- file.path(tmp_dir, "multilayer")
print(paste0('Saving the Multilayer folder to ', multilayer_f))
save_multilayer(hummus = hummus,
                folder_name = multilayer_f)

print("Running HuMMuS exploration")
grn = define_target_genes(
  hummus,
  multilayer_f = multilayer_f,
  njobs = n_cores
)
print(grn[1:5,])
#get only gene, peaks and score columns
grn <- grn[, c("gene", "tf", "score")]
colnames(grn) <- c("target", "source",  "score")

# Write
print(paste0('Saving the HuMMuS object to ', output_file))
write.csv(
  x = grn,
  file = output_file,
  row.names = FALSE,
  quote = FALSE)

print('Add enhancers to the GRN...')
# Add enhancers to the GRN
enhancers = define_binding_regions(
  hummus,
  multilayer_f = multilayer_f,
  njobs = n_cores
)[, c("tf", "peak", "score")]

print(paste0('Saving the binding regions to ', binding_regions_file))
write.csv(
  x = enhancers,
  file = binding_regions_file,
  row.names = FALSE,
  quote = FALSE)
