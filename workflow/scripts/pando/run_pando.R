library(tidyverse)
library(rhdf5)
library(Pando)
library(doParallel)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
p_thresh = args[8]
rsq_thresh = args[9]
nvar_thresh = args[10]
exclude_exons = as.logical(args[11])
path_grn <- args[12]
path_tri <- args[13]

# Set cores
registerDoParallel(parallel::detectCores())

# Set genome
if (organism == 'human'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    library(EnsDb.Hsapiens.v86)
} else {
    library(BSgenome.Mmusculus.UCSC.mm10)
    library(EnsDb.Mmusculus.v79)
}

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
muo_data <- Seurat::CreateSeuratObject(rna_X)
muo_data@assays$RNA@var.features <- genes

### ATAC
atac_indices <- indata$mod$atac$X$indices
atac_indptr <- indata$mod$atac$X$indptr
atac_data <- as.numeric(indata$mod$atac$X$data)
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE)
colnames(atac_X) <- barcodes
row.names(atac_X) <- peaks
muo_data[['peaks']] <- Signac::CreateChromatinAssay(data = atac_X)
muo_data@assays$peaks@var.features <- peaks
Seurat::DefaultAssay(object = muo_data) <- 'peaks'
h5closeAll()

## Annotate peaks
print('Add peak annotations')
if (organism == 'human'){
    gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
} else {
    gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
}
seqlevelsStyle(gene.ranges) <- "UCSC"
Signac::Annotation(muo_data[['peaks']]) <- gene.ranges

# Initiate GRN and filter regions
print('Init GRN and filter regions based on evoconv')
data('phastConsElements20Mammals.UCSC.hg38')
if (organism == 'human'){
    muo_data <- initiate_grn(
        muo_data,
        rna_assay = 'RNA',
        peak_assay = 'peaks',
        regions = phastConsElements20Mammals.UCSC.hg38,
        exclude_exons = exclude_exons
    )
} else {
    muo_data <- initiate_grn(
        muo_data,
        rna_assay = 'RNA',
        peak_assay = 'peaks',
        regions = phastConsElements20Mammals.UCSC.mm10,
        exclude_exons = exclude_exons
    )
}

# Scan Motifs
print('Scan motifs')
data('motifs')
data('motif2tf')
motif2tf_use <- motif2tf %>%
    dplyr::filter(tf %in% genes)
motifs_use <- motifs[unique(motif2tf_use$motif)]

if (organism == 'human'){
    muo_data <- find_motifs(
            muo_data, 
            pfm = motifs_use, 
            motif_tfs = motif2tf_use,
            genome = BSgenome.Hsapiens.UCSC.hg38
        )
} else {
        muo_data <- find_motifs(
            muo_data, 
            pfm = motifs_use, 
            motif_tfs = motif2tf_use,
            genome = BSgenome.Mmusculus.UCSC.mm10
        )
}

# Infer GRN
print('Infer GRN')
muo_data <- infer_grn(
    muo_data,
    peak_to_gene_method = 'GREAT',
    genes = genes,
    parallel = T
)

# Module discovery
print('Find modules')
muo_data <- find_modules(
    muo_data, 
    p_thresh = p_thresh,
    rsq_thresh = rsq_thresh,
    nvar_thresh = nvar_thresh
)

# Extract and filter
print('Extract and filter')
grn <- NetworkModules(muo_data)@meta %>%
    dplyr::rename(source = tf, weight = estimate, pvals = pval) %>%
    dplyr::select(source, target, weight, pvals, padj) %>%
    dplyr::arrange(source, target)

targets <- unique(grn$target)

tri <- coef(muo_data) %>%
    dplyr::arrange(tf, target, region) %>%
    dplyr::rename(source = tf, weight = estimate, pvals = pval) %>%
    dplyr::filter(padj < p_thresh, target %in% targets) %>%
    dplyr::select(source, target, region, weight, pvals, padj)

# Write
write.csv(grn, path_grn, row.names = FALSE)
write.csv(tri, path_tri, row.names = FALSE)
