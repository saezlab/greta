# Initiate arguments
args <- commandArgs(trailingOnly = F)
path_out <- args[6]


library(HuMMuS)
library(EnsDb.Hsapiens.v86)
library(dplyr)

# Extract TSS
gene_range = get_genome_annotations(EnsDb.Hsapiens.v86)
chr <- as.character(seqnames(gene_range))
start_pos <- start(gene_range)
end_pos <- end(gene_range)
gene_names <- mcols(gene_range)$gene_name
gene_type <- mcols(gene_range)$gene_biotype


# Build dataframe in .csv
data <- data.frame(Chromosome = chr, Start = start_pos, End = end_pos, Name = gene_names, gene.type = gene_type)


# Filter only protein coding genes
data <- data %>% filter(gene.type == "protein_coding")
data <- data %>%
  dplyr::select(Chromosome, Start, End, Name)


write.csv(x = data, file = path_out)



