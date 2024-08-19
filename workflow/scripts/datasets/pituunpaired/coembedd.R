library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)


# Parse args
args <- commandArgs(trailingOnly = F)
path_gex <- args[6]
path_peaks <- args[7]
path_frags <- args[8]
path_annot_out <- args[9]
path_gex_out <- args[10]
path_atac.se_out <- args[11]
path_cca_out <- args[12]

nCores <- 4

print(args)


# Load RNA and ATAC seq matrix
rna <- Read10X_h5(path_gex)
atac <- Read10X_h5(path_peaks)


# Add annotation to peaks 
grange.counts <- StringToGRanges(rownames(atac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac <- atac[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"


# Create Seurat object for RNA-seq data
data.rna <- CreateSeuratObject(counts = rna, project = "RNA", assay = "RNA")

# Create Seurat Object for ATAC-seq data
## Create ChromatinAssay object

frag.file <- path_frags
chrom_assay <- CreateChromatinAssay(
   counts = atac,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = frag.file,
   min.cells = 10,
   annotation = annotations
 )

## Create Seurat object from chromatin assay
data.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", project = "ATAC")



# Perform standard preprocessing of each modality independently
## RNA preprocessing
data.rna <- NormalizeData(data.rna)
data.rna <- FindVariableFeatures(data.rna)
data.rna <- ScaleData(data.rna)
data.rna <- RunPCA(data.rna)
data.rna <- FindNeighbors(data.rna, dims = 1:30)
data.rna <- FindClusters(data.rna, resolution = 0.25)
data.rna <- RunUMAP(data.rna, dims = 1:30)


# Annotate RNA clusters and remove unwanted clusters
## Manual annotation based on markers provided in paper, to be replaced with ground truth annotation as provided by authors (if available)
cells_to_remove <- WhichCells(data.rna, idents = c(5, 6))
data.rna <- subset(data.rna, cells = setdiff(Cells(data.rna), cells_to_remove))

new.cluster.ids <- c("Gonadotropes", "Stem cells", "Somatotropes", "Lactotropes", "Thyrotropes", "Pituicytes", "Stem cells", "Pericytes", "Corticotropes", 
    "Endothelial cells", "Immune cells", "Gonadotropes")

## Rename clusters
names(new.cluster.ids) <- levels(data.rna)
data.rna <- RenameIdents(data.rna, new.cluster.ids)

## Save RNA clusters to metadata
data.rna$celltype <- Idents(object = data.rna)

## Extract gene expression data and save as sparse matrix
exprMat <- GetAssayData(object = data.rna, assay = "RNA", slot = "data")
saveRDS(exprMat, file = path_gex_out)

#------------------------------------------
library(SingleCellExperiment)

# ATAC preprocessing
data.atac <- RunTFIDF(data.atac)
data.atac <- FindTopFeatures(data.atac, min.cutoff = "q0")
data.atac <- RunSVD(data.atac)
data.atac <- RunUMAP(data.atac, reduction = "lsi", dims = 2:30, reduction.name = "umap.atac", reduction.key = "atacUMAP_")

## Extract gene expression data and save as sparse matrix
atac.sce <- as.SingleCellExperiment(data.atac)
rowRanges(atac.sce) <- granges(data.atac)
saveRDS(atac.sce, file = path_atac.se_out)


#------------------------------------------

# Estimating gene activities for CCA
gene.activities <- GeneActivity(data.atac, features = VariableFeatures(data.rna))

## add gene activities as a new assay
data.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)

## normalize gene activities and compute variable features
DefaultAssay(data.atac) <- "ACTIVITY"
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac, features = rownames(data.atac))
data.atac <- FindVariableFeatures(data.atac)


# Run CCA, with default parameters --> features = union of variable features from both modalities
data.cca <- RunCCA(
  data.rna,
  data.atac,
  assay1 = "RNA",
  assay2 = "ACTIVITY",
  num.cc = 50
)

CCA_PCs <- Embeddings(data.cca, reduction = "cca")
saveRDS(CCA_PCs, file = path_cca_out)


#------------------------------------------
# Transfer cell type annotations from RNA to ATAC

## Identify anchors via CCA
transfer.anchors <- FindTransferAnchors(reference = data.rna, query = data.atac, features = VariableFeatures(object = data.rna),
    reference.assay = "RNA", query.assay = "ACTIVITY", reduction = "cca")

## Annotate cells via label transfer
celltype.predictions <- TransferData(anchorset = transfer.anchors, refdata = data.rna$celltype,
    weight.reduction = data.atac[["lsi"]], dims = 2:30)

data.atac <- AddMetaData(data.atac, metadata = celltype.predictions)

## Save ATAC clusters to metadata
data.atac$celltype <- data.atac$predicted.id

## Create annotation df
annot <- data.frame(celltype = data.atac$celltype)
annot$batch <- "smpl"
annot$barcodes <- rownames(annot)

## Save annotation df
write.csv(annot, file = path_annot_out)