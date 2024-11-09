library(biomaRt)
library(dplyr)

# Parse args
args <- commandArgs(trailingOnly = F)
window_size <- as.numeric(args[6])
out_path <- args[7]


ensembl <- useMart(
    "ensembl",
    dataset = "hsapiens_gene_ensembl",
    host = "http://www.ensembl.org"
)

gene_data <- getBM(
  attributes = c("ensembl_gene_id", "external_gene_name", "chromosome_name", "transcription_start_site"),
  mart = ensembl
)

gene_data <- gene_data %>%
  mutate(
    promoter_start = transcription_start_site - window_size,
    promoter_end = transcription_start_site + window_size - 1,
    promoter_start = pmax(promoter_start, 0)  # Ensure non-negative values
  )

standard_chromosomes <- c(1:23, "X", "Y")
bed_data <- gene_data %>%
  filter(chromosome_name %in% standard_chromosomes & external_gene_name != "") %>%
  distinct(external_gene_name, .keep_all = TRUE) %>%
  transmute(
    chrom = paste0("chr", chromosome_name),
    chromStart = promoter_start - 1,  # BED format is 0-based
    chromEnd = promoter_end,
    name = external_gene_name
  )

bed_data <- bed_data %>%
  arrange(
    factor(chrom, levels = paste0("chr", c(1:23, "X", "Y"))),
    chromStart
  )

# Write to output file
write.table(bed_data, file = out_path, sep = "\t", quote = FALSE, col.names = FALSE, row.names = FALSE)
