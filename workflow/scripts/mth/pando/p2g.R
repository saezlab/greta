library(tidyverse)
library(rhdf5)
library(Pando)
library(GenomicRanges)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_ann <- args[7]
extend <- as.numeric(args[8])
path_out <- args[9]

# Set genome
annot <- read.csv(path_ann)
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'

# Read peaks and genes
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
peaks <- indata$mod$atac$var$`_index`
genes <- indata$mod$rna$var$`_index`
h5closeAll()
peaks <- data.frame(seqnames=peaks)
peaks <- tidyr::separate(data = peaks, col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Filter annot by seen genes
annot <- annot[annot$gene_name %in% intersect(genes, annot$gene_name), ]

# Find peak2gene links
peaks_near_gene <- find_peaks_near_genes(
    peaks = peaks,
    genes = annot,
    method = 'GREAT',
    upstream = round(extend / 2),
    downstream = round(extend / 2),
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
