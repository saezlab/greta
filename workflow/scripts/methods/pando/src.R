library(tidyverse)
library(rhdf5)
library(Pando)
library(doParallel)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
exclude_exons = as.logical(args[7])
ext = as.numeric(args[8])
thr_corr = as.numeric(args[9])
p_thresh = as.numeric(args[10])
rsq_thresh = as.numeric(args[11])
nvar_thresh = as.numeric(args[12])
min_genes_per_module = as.numeric(args[13])
organism <- args[14]
granges_hg <- args[15]
granges_mm <- args[16]
path_grn <- args[17]

# Set cores
registerDoParallel(32)

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`

### ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
h5closeAll()

# Create object
muo_data <- Seurat::CreateSeuratObject(rna_X)
muo_data[['peaks']] <- Signac::CreateChromatinAssay(data = atac_X)

# Set genome
if (organism == 'hg38'){
    data('phastConsElements20Mammals.UCSC.hg38')
    annot <- read.csv(granges_hg)
    regions <- phastConsElements20Mammals.UCSC.hg38
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38
} else if (organism == 'mm10'){
    annot <- read.csv(granges_mm)
    regions <- phastConsElements20Mammals.UCSC.mm10
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- BSgenome.Mmusculus.UCSC.mm10
}
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'
Signac::Annotation(muo_data[['peaks']]) <- annot

# Init GRN
print('Init GRN and filter regions based on evoconv')
muo_data <- initiate_grn(
    muo_data,
    rna_assay = 'RNA',
    peak_assay = 'peaks',
    regions = regions,
    exclude_exons = exclude_exons
)

# Update object with only peaks that have genes linked
#peaks_near_gene <- find_peaks_near_genes(
#    peaks = NetworkRegions(muo_data)@ranges,
#    genes = annot[annot$gene_name %in% intersect(rownames(muo_data[['RNA']]), annot$gene_name), ],
#    method = 'GREAT',
#    upstream = round(ext / 2),
#    downstream = round(ext / 2),
#)
#peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups=colnames(peaks_near_gene), fun='sum')
#msk <- as.logical(sparseMatrixStats::colMaxs(peaks2gene))
#p2g_peaks <- colnames(peaks2gene)[msk]
#muo_data@grn@regions@ranges <- data.frame(seqnames=p2g_peaks) %>%
#    tidyr::separate(col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE) %>%
#    GenomicRanges::makeGRangesFromDataFrame(.)
#muo_data@grn@regions@peaks <- subjectHits(findOverlaps(muo_data@grn@regions@ranges, Signac::StringToGRanges(rownames(Seurat::GetAssay(muo_data, assay='peaks')))))

# Find motifs
print('Scan motifs')
data('motifs')
data('motif2tf')
motif2tf_use <- motif2tf %>%
    dplyr::filter(tf %in% rownames(muo_data[['RNA']]))
motifs_use <- motifs[unique(motif2tf_use$motif)]
muo_data <- find_motifs(
    muo_data, 
    pfm = motifs_use, 
    motif_tfs = motif2tf_use,
    genome = genome
)

# Infer GRN
print('Infer GRN')
muo_data <- infer_grn(
    muo_data,
    peak_to_gene_method = 'GREAT',
    genes = rownames(muo_data[['RNA']]),
    tf_cor = thr_corr,
    peak_cor = thr_corr,
    upstream = round(ext / 2),
    downstream = round(ext / 2),
    parallel = T
)

# Module discovery
print('Find modules')
muo_data <- find_modules(
    muo_data, 
    p_thresh = p_thresh,
    rsq_thresh = rsq_thresh,
    nvar_thresh = nvar_thresh,
    min_genes_per_module = min_genes_per_module
)

# Extract and filter
print('Extract and filter')
grn <- NetworkModules(muo_data)@meta %>%
    dplyr::select(-pval) %>%
    dplyr::rename(source = tf, score = estimate, pval = padj) %>%
    dplyr::select(source, target, score, pval) %>%
    dplyr::arrange(source, target)

# Write
write.csv(x = grn, file = path_grn, row.names=FALSE)
