library(EnsDb.Hsapiens.v86)
library(dplyr)


# Parse args
args <- commandArgs(trailingOnly = F)
path_out <- args[6]


# Read
gr <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)

# Merge overlaps
merged <- unlist(reduce(split(gr, gr$gene_name)), use.names = TRUE)

# To df
chr_names <- paste0("chr", as.character(seqnames(merged)))
start_pos <- start(merged)
end_pos <- end(merged)
gene_names <- names(merged)
bed <- data.frame(Chromosome = chr_names, Start = start_pos, End = end_pos, Name = gene_names)
bed <- dplyr::arrange(bed, Chromosome, Start, End)

# Write
write.table(x = bed, file = path_out, sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
