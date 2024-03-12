library(tidyverse)
library(rhdf5)
library(Pando)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
granges_hg <- args[8]
granges_mm <- args[9]
exclude_exons <- args[10]
path_out <- args[11]

# Set genome
if (organism == 'hg38'){
    data('phastConsElements20Mammals.UCSC.hg38')
    annot <- read.csv(granges_hg)
    regions <- phastConsElements20Mammals.UCSC.hg38
} else if (organism == 'mm10'){
    annot <- read.csv(granges_mm)
    regions <- phastConsElements20Mammals.UCSC.mm10
}
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'

# Read peaks
indata <- H5Fopen(path_data)
peaks <- indata$mod$atac$var$`_index`
h5closeAll()
peaks <- data.frame(seqnames=peaks)
peaks <- tidyr::separate(data = peaks, col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Filter by evo cons regions
hits <- GenomicRanges::findOverlaps(regions, peaks)
cand <- GenomicRanges::pintersect(
    peaks[S4Vectors::subjectHits(hits)],
    regions[S4Vectors::queryHits(hits)]
)

# Substract exons
if (exclude_exons){
    exons <- annot[annot$type=='exon', ]
    names(exons@ranges) <- NULL
    exons <- IRanges::intersect(exons, exons)
    exons <- GenomicRanges::GRanges(
        seqnames = exons@seqnames,
        ranges = exons@ranges
    )
    cand <- GenomicRanges::subtract(cand, exons, ignore.strand=TRUE) %>% unlist()
}

# Filter peaks by candidate regions
matches <- S4Vectors::subjectHits(GenomicRanges::findOverlaps(cand, peaks))
peaks <- peaks[unique(matches),]

# Write
write.csv(x = peaks, file = path_out, row.names=FALSE)
