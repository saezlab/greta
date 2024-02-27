library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
library(rhdf5)
library(presto)
options(stringsAsFactors = FALSE)

# Set genome (copied from pando)
if (organism == 'human'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(EnsDb.Hsapiens.v86)
} else {
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(EnsDb.Mmusculus.v79)
}

#For testing
path_data = c("../resources/neurips2021/small/mdata.h5mu")
organism = c('human')

## --- Copied from Hummus - Needs adaption and checking ---
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



# Set up data - checked
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

### Extract celltype code
celltypes <- indata$obs$celltype$codes
Idents(muo_data) <- celltype
muo_data$celltype <- Idents(muo_data) 

### ATAC - checked
atac_indices <- indata$mod$atac$X$indices
atac_indptr <- indata$mod$atac$X$indptr
atac_data <- as.numeric(indata$mod$atac$X$data)
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE)
colnames(atac_X) <- barcodes
row.names(atac_X) <- peaks
muo_data[['ATAC']] <- Signac::CreateChromatinAssay(counts = atac_X)
h5closeAll()



## Annotate peaks - checked
print('Add peak annotations')
if (organism == 'human'){
    gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
} else {
    gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}
seqlevelsStyle(gene.ranges) <- "UCSC"
Signac::Annotation(muo_data[['ATAC']]) <- gene.ranges



## Transform RNA data - checked
DefaultAssay(muo_data) <- "RNA"
muo_data <- SCTransform(muo_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')



## Integrate data 
DefaultAssay(muo_data) <- "ATAC"
muo_data <- RunTFIDF(muo_data)
muo_data <- FindTopFeatures(muo_data, min.cutoff = 'q0')
muo_data <- RunSVD(muo_data)
muo_data <- RunUMAP(muo_data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
# had to modify from DIRECTNET tutorial; default k.nn=20, knn.range=200
muo_data <- FindMultiModalNeighbors(muo_data, k.nn = 20, knn.range = 100, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 1:50))
muo_data <- RunUMAP(muo_data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")


## Find markers 
markers_all <- presto:::wilcoxauc.Seurat(X = muo_data, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_all <- markers_all[which(markers_all$auc > 0.5), , drop = FALSE]
markers <- data.frame(gene = markers_all$feature, group = markers_all$group)
c <- unique(markers$group)
marker_list <- list()
for (i in 1:length(c)) {
  marker1<- markers_all[markers$group == c[i],]
  marker_list[[i]] <- as.character(marker1$feature[marker1$auc > 0.5])
}
markers_groups <- unique(unlist(marker_list))
markers_groups <- lapply(markers_groups, function(x) strsplit(x,"[.]")[[1]][1])
markers_groups <- unique(unlist(markers_groups))

                         
                         
## Infer CREs
muo_data <- Run_DIRECT_NET(muo_data, peakcalling = FALSE, k_neigh = 50, atacbinary = TRUE, max_overlap=0.5, size_factor_normalize = FALSE, genome.info = gene.ranges, focus_markers = markers_groups)
