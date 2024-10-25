# Add arguments
args <- commandArgs(trailingOnly = F)
path_out <- args[6]

library(FigR)
library(GenomicRanges)

# Extract TSS annotations
TSSg <- FigR::hg38TSSRanges
chr <- as.character(seqnames(TSSg))
start_pos <- start(TSSg)
end_pos <- end(TSSg)
gene_names <- mcols(TSSg)$gene_name

# Transform it into a data frame
data <- data.frame(Chromosome = chr, Start = start_pos, End = end_pos, Name = gene_names)

write.csv(x = data, file = path_out)





