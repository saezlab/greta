library(FigR)
library(GenomicRanges)


# Add arguments
args <- commandArgs(trailingOnly = F)
path_out <- args[6]

# Extract TSS annotations
TSSg <- FigR::hg38TSSRanges
chr <- as.character(seqnames(TSSg))
start_pos <- start(TSSg)
end_pos <- end(TSSg)
gene_names <- mcols(TSSg)$gene_name

# Transform it into a data frame
data <- data.frame(Chromosome = chr, Start = start_pos - 1, End = end_pos - 1, Name = gene_names)

# Write
write.table(x = data, file = path_out, sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
