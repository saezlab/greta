library(tidyverse)
library(rhdf5)
library(Pando)
library(doParallel)
library(BSgenome.Hsapiens.UCSC.hg38)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_gof_plot <- args[8]
path_modules_plot <- args[9]
path_topo_plot <- args[10]

# Set cores
registerDoParallel(parallel::detectCores())

# Set up data

## Read data
indata <- H5Fopen(path_data)

### RNA
rna_indices <- indata$mod$rna$raw$X$indices
rna_indptr <- indata$mod$rna$raw$X$indptr
rna_data <- as.numeric(indata$mod$rna$raw$X$data)
barcodes <- indata$obs$`_index`
genes <- indata$mod$rna$var$`_index`
rna_X <- Matrix::sparseMatrix(i=rna_indices, p=rna_indptr, x=rna_data, index1 = FALSE)
colnames(rna_X) <- barcodes
row.names(rna_X) <- genes
muo_data <- Seurat::CreateSeuratObject(rna_X)

### ATAC
atac_indices <- indata$mod$atac$raw$X$indices
atac_indptr <- indata$mod$atac$raw$X$indptr
atac_data <- as.numeric(indata$mod$atac$raw$X$data)
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE)
colnames(atac_X) <- barcodes
row.names(atac_X) <- peaks
muo_data[['peaks']] <- Signac::CreateChromatinAssay(counts = atac_X, genome = "hg19")

## Annotate peaks
gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v75)
seqlevelsStyle(gene.ranges) <- "UCSC"
Signac::Annotation(muo_data[['peaks']]) <- gene.ranges

# Initiate GRN and filter regions
data('phastConsElements20Mammals.UCSC.hg38')
muo_data <- initiate_grn(
    muo_data,
    rna_assay = 'RNA',
    peak_assay = 'peaks',
    regions = phastConsElements20Mammals.UCSC.hg38 
)

# Scan Motifs
data('motifs')
data('motif2tf')
patterning_genes <- read_tsv('patterning_genes.tsv')
pattern_tfs <- patterning_genes %>% 
    filter(type=='Transcription factor') %>% 
    pull(symbol)
motif2tf_use <- motif2tf %>%
    filter(tf %in% pattern_tfs)
motifs_use <- motifs[unique(motif2tf_use$motif)]
muo_data <- find_motifs(
    muo_data, 
    pfm = motifs_use, 
    motif_tfs = motif2tf_use,
    genome = BSgenome.Hsapiens.UCSC.hg38
)

# Infer GRN
muo_data <- infer_grn(
    muo_data,
    peak_to_gene_method = 'GREAT',
    genes = patterning_genes$symbol,
    parallel = T
)

# Module discovery
muo_data <- find_modules(
    muo_data, 
    p_thresh = 0.1,
    nvar_thresh = 2, 
    min_genes_per_module = 1, 
    rsq_thresh = 0.05
)

# Extract network
modules <- NetworkModules(muo_data) 

# Write
write.csv(modules@meta, path_grn)

# Plotting

## Goodness-of-fit
pdf(file=path_gof_plot, width=8, height=4)
plot_gof(muo_data, point_size=3)
dev.off()

## Distirbution modules
pdf(file=path_modules_plot, width=8, height=4)
plot_module_metrics(muo_data)
dev.off()

## GRN topology
pdf(file=path_topo_plot, width=8, height=4)
plot_network_graph(muo_data)
dev.off()

