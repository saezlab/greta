#!/usr/bin/env Rscript

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 2) {
  stop("Usage: rds_to_annot.R IN_RDS OUT_CSV")
}

suppressWarnings(suppressMessages(library(methods)))

obj <- readRDS(args[1])

if (!isS4(obj) || !"meta.data" %in% slotNames(obj)) {
  stop("Expected Seurat object with meta.data")
}

meta <- slot(obj, "meta.data")

if (!"celltype_major" %in% colnames(meta) || !"library" %in% colnames(meta)) {
  stop("Missing celltype_major or library columns")
}

barcodes <- rownames(meta)
celltype <- as.character(meta$celltype_major)
library <- as.character(meta$library)

# Remove dots and underscores from batch names
batch <- gsub("[._]", "", library)

# Prepend batch to barcode (remove -X suffix from barcode)
final_barcodes <- character(length(barcodes))
for (i in seq_along(barcodes)) {
  bc <- barcodes[i]
  parts <- strsplit(bc, "-")[[1]]
  barcode_base <- parts[1]
  final_barcodes[i] <- paste0(batch[i], "_", barcode_base)
}

outdf <- data.frame(
  batch = batch,
  celltype = celltype,
  stringsAsFactors = FALSE
)
rownames(outdf) <- final_barcodes

write.csv(outdf, args[2], row.names = TRUE, quote = FALSE)
cat(sprintf("[rds_to_annot] Wrote %d cells\n", nrow(outdf)), file = stderr())
