library(Signac)
library(EnsDb.Hsapiens.v86)
library(ggplot2)
library(cowplot)
library(dplyr)
library(Seurat)


# Parse args
args <- commandArgs(trailingOnly = F)
path_gex <- args[6]
path_celltypes <- args[7]
path_peaks <- args[8]
path_frags <- args[9]
path_cca_out <- args[10]


# RNA
rna <- Read10X_h5(path_gex)
data.rna <- CreateSeuratObject(counts = rna, project = "RNA", assay = "RNA")
celltypes <- read.csv(path_celltypes)
cells_to_remove <- Cells(data.rna)[!Cells(data.rna) %in% celltypes$X]
data.rna <- subset(data.rna, cells = setdiff(Cells(data.rna), cells_to_remove))
data.rna <- NormalizeData(data.rna)
data.rna <- FindVariableFeatures(data.rna)
data.rna <- ScaleData(data.rna)

# ATAC
atac <- Read10X_h5(path_peaks)
grange.counts <- StringToGRanges(rownames(atac), sep = c(":", "-"))
grange.use <- seqnames(grange.counts) %in% standardChromosomes(grange.counts)
atac <- atac[as.vector(grange.use), ]
annotations <- GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
seqlevelsStyle(annotations) <- 'UCSC'
genome(annotations) <- "hg38"
colnames(atac) <- gsub("-[0-9]+$", "", colnames(atac))
colnames(atac) <- paste0("smpl_", colnames(atac))
chrom_assay <- CreateChromatinAssay(
   counts = atac,
   sep = c(":", "-"),
   genome = 'hg38',
   fragments = path_frags,
   min.cells = 10,
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

