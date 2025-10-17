#!/usr/bin/env Rscript

# Extract data from Seurat RDS files for retina dataset
# Usage: Rscript extract_rds.R <atac_rds> <rna_rds> <multiomics_rds> <output_dir>

suppressPackageStartupMessages({
  library(Seurat)
  library(Matrix)
})

args <- commandArgs(trailingOnly = TRUE)
atac_rds <- args[1]
rna_rds <- args[2]
multiomics_rds <- args[3]
output_dir <- args[4]
sample_name <- args[5]

cat("Loading ATAC RDS...\n")
atac_obj <- readRDS(atac_rds)

cat("Loading RNA RDS...\n")
rna_obj <- readRDS(rna_rds)

cat("Loading multiomics RDS...\n")
multiomics_obj <- readRDS(multiomics_rds)

# Extract cell type annotations from RNA object's metadata
# The celltype column in metadata has the actual cell type names, not just cluster numbers
if ("celltype" %in% colnames(rna_obj@meta.data)) {
  celltypes <- as.character(rna_obj@meta.data$celltype)
  names(celltypes) <- rownames(rna_obj@meta.data)
  cat("Extracted cell types from metadata celltype column\n")
  cat("Unique cell types:", length(unique(celltypes)), "\n")
  cat("Cell type examples:", paste(head(unique(celltypes)), collapse=", "), "\n")
} else {
  # Fallback to active.ident if celltype column doesn't exist
  celltypes <- as.character(rna_obj@active.ident)
  names(celltypes) <- names(rna_obj@active.ident)
  cat("Using active.ident as celltype (no celltype column found)\n")
}

# Get barcodes - they should match
atac_barcodes <- colnames(atac_obj)
rna_barcodes <- colnames(rna_obj)

cat("ATAC cells:", length(atac_barcodes), "\n")
cat("RNA cells:", length(rna_barcodes), "\n")

# Check if barcodes match
if (!all(atac_barcodes == rna_barcodes)) {
  stop("ATAC and RNA barcodes do not match!")
}

# Remove the -1 suffix and prepend sample name
clean_barcodes <- gsub("-1$", "", atac_barcodes)
final_barcodes <- paste0(sample_name, "_", clean_barcodes)

# Create annotation dataframe
annot <- data.frame(
  row.names = final_barcodes,
  batch = sample_name,
  celltype = celltypes[atac_barcodes]
)

# Write annotation
annot_path <- file.path(output_dir, "annot.csv")
write.csv(annot, annot_path, quote = FALSE)
cat("Written annotation to:", annot_path, "\n")

# Extract RNA counts from RNA object
rna_counts <- GetAssayData(rna_obj, slot = "counts", assay = "RNA")
rna_genes <- rownames(rna_counts)

# Rename columns to match new barcodes
colnames(rna_counts) <- final_barcodes

# Write RNA counts as sparse matrix
rna_path <- file.path(output_dir, "rna_counts.mtx")
writeMM(rna_counts, rna_path)
cat("Written RNA counts to:", rna_path, "\n")

# Write RNA gene names
genes_path <- file.path(output_dir, "rna_genes.txt")
writeLines(rna_genes, genes_path)
cat("Written RNA genes to:", genes_path, "\n")

# Write RNA barcodes
rna_barcodes_path <- file.path(output_dir, "rna_barcodes.txt")
writeLines(final_barcodes, rna_barcodes_path)
cat("Written RNA barcodes to:", rna_barcodes_path, "\n")

# Extract ATAC peaks from multiomics object (already called peaks)
# The multiomics object has "Peaks" which is the peak count matrix
peaks_matrix <- multiomics_obj$Peaks

# Verify dimensions match expected cells
if (ncol(peaks_matrix) != length(rna_barcodes)) {
  cat("WARNING: Peaks matrix has different number of cells than RNA/ATAC objects\n")
  cat("Peaks cells:", ncol(peaks_matrix), "RNA/ATAC cells:", length(rna_barcodes), "\n")
  # Filter to matching cells
  matching_cells <- intersect(colnames(peaks_matrix), atac_barcodes)
  cat("Matching cells:", length(matching_cells), "\n")
  peaks_matrix <- peaks_matrix[, matching_cells]
  # Update barcodes
  matched_indices <- match(matching_cells, atac_barcodes)
  peak_barcodes <- final_barcodes[matched_indices]
} else {
  peak_barcodes <- final_barcodes
}

# Rename columns to match new barcodes
colnames(peaks_matrix) <- peak_barcodes

# Get peak coordinates
peak_names <- rownames(peaks_matrix)

# Write peaks matrix
peaks_path <- file.path(output_dir, "peaks_counts.mtx")
writeMM(peaks_matrix, peaks_path)
cat("Written peaks counts to:", peaks_path, "\n")

# Write peak names
peaks_names_path <- file.path(output_dir, "peak_names.txt")
writeLines(peak_names, peaks_names_path)
cat("Written peak names to:", peaks_names_path, "\n")

# Write peak barcodes
peaks_barcodes_path <- file.path(output_dir, "peak_barcodes.txt")
writeLines(peak_barcodes, peaks_barcodes_path)
cat("Written peak barcodes to:", peaks_barcodes_path, "\n")

cat("\n=== Extraction complete ===\n")
cat("Cells:", length(final_barcodes), "\n")
cat("Genes:", length(rna_genes), "\n")
cat("Peaks:", nrow(peaks_matrix), "\n")
cat("Cell types:", length(unique(celltypes)), "\n")
