library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)
library(SingleCellExperiment)
library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_gex <- args[6]
path_peaks <- args[7]
path_frags <- args[8]
path_cca_out <- args[9]


# Load RNA and ATAC seq matrix

# Process RNA
rna <- Read10X_h5(path_gex)[[1]]
data.rna <- CreateSeuratObject(counts = rna, project = "RNA", assay = "RNA")
data.rna <- NormalizeData(data.rna)
data.rna <- FindVariableFeatures(data.rna)
data.rna <- ScaleData(data.rna)

# Process ATAC
indata <- H5Fopen(path_peaks, flags='H5F_ACC_RDONLY')
indices <- indata$X$indices
indptr <- indata$X$indptr
data <- as.numeric(indata$X$data)
atac <- Matrix::sparseMatrix(i=indices, p=indptr, x=data, index1 = FALSE)
colnames(atac) <- indata$obs$`_index`
rownames(atac) <- indata$var$`_index`
h5closeAll()
grange.counts <- StringToGRanges(rownames(atac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac <- atac[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
chrom_assay <- CreateChromatinAssay(
   counts = atac,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = path_frags,
   annotation = annotations
)
data.atac <- CreateSeuratObject(counts = chrom_assay, assay = "ATAC", project = "ATAC")
data.atac <- RunTFIDF(data.atac)
data.atac <- FindTopFeatures(data.atac, min.cutoff = "q0")
data.atac <- ScaleData(data.atac)

# Infer gene scores
gene.activities <- GeneActivity(data.atac, features = VariableFeatures(data.rna))
data.atac[["ACTIVITY"]] <- CreateAssayObject(counts = gene.activities)
DefaultAssay(data.atac) <- "ACTIVITY"
data.atac <- NormalizeData(data.atac)
data.atac <- ScaleData(data.atac, features = rownames(data.atac))
data.atac <- FindVariableFeatures(data.atac)

# Run CCA
data.cca <- RunCCA(
  data.rna,
  data.atac,
  assay1 = "RNA",
  assay2 = "ACTIVITY",
  num.cc = 50
)

CCA_PCs <- Embeddings(data.cca, reduction = "cca")
saveRDS(CCA_PCs, file = path_cca_out)
