library(tidyverse)
library(rhdf5)
library(Pando)
library(GenomicRanges)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_ann <- args[7]
exclude_exons <- args[8]
path_cand <- args[9]
path_matches <- args[10]

# Set genome
data('phastConsElements20Mammals.UCSC.hg38')
regions <- phastConsElements20Mammals.UCSC.hg38
annot <- read.csv(path_ann)
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'

# Read peaks
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
peaks <- indata$mod$atac$var$`_index`
h5closeAll()
peaks <- data.frame(seqnames=peaks)
peaks <- tidyr::separate(data = peaks, col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Read exons
exons <- annot[annot$type=='exon', ]
names(exons@ranges) <- NULL
exons <- IRanges::intersect(exons, exons)
exons <- GenomicRanges::GRanges(
    seqnames = exons@seqnames,
    ranges = exons@ranges
)

# Intersect by only shared chromosomes
seqnames <- intersect(intersect(levels(peaks@seqnames), levels(regions@seqnames)), levels(exons@seqnames))
peaks <- keepSeqlevels(peaks, seqnames, pruning.mode = "coarse")
exons <- keepSeqlevels(exons, seqnames, pruning.mode = "coarse")
regions <- keepSeqlevels(regions, seqnames, pruning.mode = "coarse")

# Filter by evo cons regions
hits <- GenomicRanges::findOverlaps(regions, peaks)
cand <- GenomicRanges::pintersect(
    peaks[S4Vectors::subjectHits(hits)],
    regions[S4Vectors::queryHits(hits)]
)

# Substract exons
if (exclude_exons){
    cand <- GenomicRanges::subtract(cand, exons, ignore.strand=TRUE) %>% unlist()
}

# Find matches of new peaks to old peaks
matches <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(cand, peaks))

# Write
write.csv(x = cand, file = path_cand, row.names=FALSE)
write.csv(x = matches, file = path_matches, row.names=FALSE)
