library(tidyverse)
library(rhdf5)
library(Pando)
library(GenomicRanges)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
granges_hg <- args[8]
granges_mm <- args[9]
extend <- as.numeric(args[10])
path_out <- args[11]

# Set genome
if (organism == 'hg38'){
    annot <- read.csv(granges_hg)
} else if (organism == 'mm10'){
    annot <- read.csv(granges_mm)
}
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'

# Read peaks and genes
indata <- H5Fopen(path_data)
peaks <- indata$mod$atac$var$`_index`
genes <- indata$mod$rna$var$`_index`
h5closeAll()
peaks <- data.frame(seqnames=peaks)
peaks <- tidyr::separate(data = peaks, col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Filter annot by seen genes
annot <- annot[annot$gene_name %in% intersect(genes, annot$gene_name)]

# Find peak2gene links
peaks_near_gene <- find_peaks_near_genes(
    peaks = peaks,
    genes = annot,
    method = 'GREAT',
    upstream = 100000,
    downstream = 0,
    extend = extend,
)
peaks2gene <- aggregate_matrix(t(peaks_near_gene), groups=colnames(peaks_near_gene), fun='sum')

# Convert from sparse mat to df
sparse <- summary(peaks2gene)
df <- data.frame(
  cre = colnames(peaks2gene)[sparse$j],
  gene = rownames(peaks2gene)[sparse$i],
  score = sparse$x
)
df <- df %>% arrange(cre, desc(score))

# Write
write.csv(x = df, file = path_out, row.names=FALSE)
