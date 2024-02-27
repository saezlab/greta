library(rhdf5)
library(HuMMuS)

# Parse args
args <- commandArgs(trailingOnly = FALSE)
path_data <- args[6]
genie3_cores <- as.integer(args[7])
organism <- args[8]
genie3_thresh <- as.integer(args[9])
cicero_k_per_pseudocells <- as.integer(args[10])
cicero_thresh <- as.double(args[11])
upstream_gene <- as.integer(args[12])
downstream_gene <- as.integer(args[13])
only_tss <- as.logical(args[14])
path_multilayer <- args[15]
path_grn <- args[16]
path_tri <- args[17]

print(only_tss)

# Set up data
print('Open object')
## Read data
indata <- H5Fopen(path_data)
### RNA
rna_indices <- indata$mod$rna$X$indices
rna_indptr <- indata$mod$rna$X$indptr
rna_data <- as.numeric(indata$mod$rna$X$data)
barcodes <- indata$obs$`_index`
genes <- indata$mod$rna$var$`_index`
rna_X <- Matrix::sparseMatrix(i=rna_indices, p=rna_indptr, x=rna_data, index1 = FALSE)
colnames(rna_X) <- barcodes
row.names(rna_X) <- genes
muo_data <- SeuratObject::CreateSeuratObject(rna_X)
muo_data@assays$RNA@var.features <- genes

### ATAC
atac_indices <- indata$mod$atac$X$indices
atac_indptr <- indata$mod$atac$X$indptr
atac_data <- as.numeric(indata$mod$atac$X$data)
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE)
colnames(atac_X) <- barcodes
row.names(atac_X) <- peaks
muo_data[['peaks']] <- Signac::CreateChromatinAssay(counts = atac_X)
h5closeAll()

## Select annotations and genome
print('Add peak annotations')
if (organism == 'human') {
    genome_annotations <- get_genome_annotations(
        ensdb_annotations = EnsDb.Hsapiens.v86::EnsDb.Hsapiens.v86)
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    omnipath_organism <- "9606"
} else {
    genome_annotations <- get_genome_annotations(
        ensdb_annotations = EnsDb.Mmusculus.v79::EnsDb.Mmusculus.v79)
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    omnipath_organism <- "10090"
}

### Add the annotations to peaks
### get human genome annotation from EndDb data
### wrapper of Signac::GetGRangesFromEnsDb, adapting output to UCSC format

# Create hummus object
print("Create hummus object")
hummus <- as(muo_data, "hummus_object")

# We add annotation to the ATAC/peak assay
# (will be used to define peak-gene links)
Signac::Annotation(hummus@assays$peaks) <- genome_annotations
rm(genome_annotations)

# Add TFs motifs info
hummus@motifs_db <- get_tf2motifs(species = organism)

## Compute gene network
print("Create gene network")
hummus <- compute_gene_network(
            hummus,
            gene_assay = "RNA",
            threshold = genie3_thresh,
            number_cores = genie3_cores
            )

## Compute peak network
print(hummus[['peaks']][1:10, 1:10])
print("Create peak network")
hummus <- compute_atac_peak_network(
            hummus,
            atac_assay = "peaks",
            number_cells_per_clusters = cicero_k_per_pseudocells,
            genome = genome,
            threshold = cicero_thresh
            )

## Create "fake" TF network
print("Create TF network")
hummus <- compute_tf_network(
            hummus = hummus,
            organism = omnipath_organism
            )

## Link TF to peaks
print("Link TF to peaks")
hummus <- bipartite_tfs2peaks(
            hummus,
            tf_expr_assay = "RNA",
            peak_assay = "peaks",
            genome = genome,
            tf_multiplex_name = "TF"
            )

## Link peak to genes
print("Link peaks to genes")
hummus <- bipartite_peaks2genes(
            hummus,
            gene_assay = "RNA",
            peak_assay = "peaks",
            only_tss = only_tss,
            upstream = upstream_gene,
            downstream = downstream_gene
            )

tryCatch(
        {
        print(hummus@multilayer)
        print(hummus@multilayer@bipartites)
        print(hummus@multilayer@multiplex)
        },
        error=function(e) {
            message('An Error Occurred')
            print(e)
        },
        warning=function(w) {
            message('A Warning Occurred')
            print(w)
            return(NA)
        }
    )

# Save multilayer
save_multilayer(hummus,
                path_multilayer)

# Infer GRN and enhancers
grn <- define_grn(hummus,
                  multilayer_f = path_multilayer)
tri <- define_enhancers(hummus,
                        multilayer_f = path_multilayer)

# Write
write.csv(grn, path_grn, row.names = FALSE)
write.csv(tri, path_tri, row.names = FALSE)

