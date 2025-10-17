#!/usr/bin/env Rscript

# Simple script to fetch Ensembl↔MGI↔synonym map for mouse (mm10)
# Uses biomaRt and outputs: ensembl_gene_id, symbol, synonyms

args <- commandArgs(trailingOnly = TRUE)
out_path <- if (length(args) >= 1) args[1] else "mm10_gene_map.tsv"

suppressMessages({
  library(biomaRt)
  library(dplyr)
})

# Connect to Ensembl mouse dataset (works fine for mm10)
mart <- useEnsembl(biomart = "genes", dataset = "mmusculus_gene_ensembl")

# Grab identifiers and synonyms - expanded to include more identifier types
atts <- c("ensembl_gene_id", "mgi_symbol", "external_gene_name", 
          "external_synonym", "entrezgene_id", "description")
x <- getBM(attributes = atts, mart = mart)

# Preferred symbol = MGI if present, otherwise external_gene_name
x$symbol <- ifelse(is.na(x$mgi_symbol) | x$mgi_symbol == "",
                   x$external_gene_name, x$mgi_symbol)

# Add external_gene_name as an additional synonym if it differs from symbol
x$all_synonyms <- sapply(1:nrow(x), function(i) {
  syns <- c()
  if (!is.na(x$external_synonym[i]) && x$external_synonym[i] != "") {
    syns <- c(syns, x$external_synonym[i])
  }
  if (!is.na(x$external_gene_name[i]) && x$external_gene_name[i] != "" && 
      x$external_gene_name[i] != x$symbol[i]) {
    syns <- c(syns, x$external_gene_name[i])
  }
  if (!is.na(x$mgi_symbol[i]) && x$mgi_symbol[i] != "" && 
      x$mgi_symbol[i] != x$symbol[i]) {
    syns <- c(syns, x$mgi_symbol[i])
  }
  paste(unique(syns), collapse = "|")
})

# Collapse synonyms per gene
map <- x %>%
  group_by(ensembl_gene_id, symbol) %>%
  summarise(
    synonyms = paste(unique(na.omit(all_synonyms[all_synonyms != ""])), collapse = "|"),
    .groups = "drop"
  )

dir.create(dirname(out_path), showWarnings = FALSE, recursive = TRUE)
write.table(map, file = out_path, sep = "\t", quote = FALSE, row.names = FALSE)
cat("[DONE] Wrote gene map to", out_path, "\n")
cat("Total genes in map:", nrow(map), "\n")