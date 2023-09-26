#!/usr/local/bin/Rscript

library(argparse)
library(rhdf5)
library(Signac)
library(Matrix)
library(Seurat)

# ./1.preprocess.R \
#--input "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/mdata.h5mu" \
#--output "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/counts_RNA_raw.tsv.gz" "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/counts_ATAC_raw.tsv.gz" "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/metadata.tsv.gz" \
# --organism human \
# --preprocessing_countAggregation mean \
# --preprocessing_minCellsPerCluster 25 \
# --preprocessing_minClusters 30


parser <- ArgumentParser(description= 'GRaNIE preprocessing')
parser$add_argument('--input', '-i', help= 'Input file that contains both raw RNA and raw ATAC counts (h5 file)', nargs = 1,  required= TRUE)
parser$add_argument('--output', '-o', help= 'A list of exactly 3 output files: 1. RNA counts, 2. ATAC counts, 3. metadata', nargs= 3, required= TRUE)
parser$add_argument('--organism', help= 'Organism', nargs= 1, required= TRUE)
parser$add_argument('--preprocessing_countAggregation', help= 'Count aggregation (mean or sum)', nargs= 1, required= TRUE)
parser$add_argument('--preprocessing_minCellsPerCluster', help= 'minCellsPerCluster', nargs= 1, required= TRUE)
parser$add_argument('--preprocessing_minClusters', help= 'minClusters', nargs= 1, required= TRUE)


xargs <- parser$parse_args()

outdir <-dirname(xargs$output[1])

# Read genome
if (xargs$organism == 'human'){
    file_RNA_features = "resources/GRaNIE/features_RNA_hg38.tsv.gz"
} else {
    stop("Not implemented yet")
}


cat('Open file ', xargs$input, "\n")
## Read data
indata <- H5Fopen(xargs$input)

### RNA
rna_indices <- indata$mod$rna$layers$counts$indices
rna_indptr <- indata$mod$rna$layers$counts$indptr
rna_data <- as.numeric(indata$mod$rna$layers$counts$data)
barcodes <- indata$obs$`_index`
genes <- indata$mod$rna$var$`_index`
rna_X <- Matrix::sparseMatrix(i=rna_indices, p=rna_indptr, x=rna_data, index1 = FALSE, dimnames = list(genes, barcodes))
seu.s <- Seurat::CreateSeuratObject(rna_X, assay = "RNA")

### ATAC
atac_indices <- indata$mod$atac$layers$counts$indices
atac_indptr <- indata$mod$atac$layers$counts$indptr
atac_data <- as.numeric(indata$mod$atac$layers$counts$data)
barcodes <- indata$mod$atac$obs$`_index`
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE, dimnames = list(peaks, barcodes))
seu.s[['ATAC']] <- Signac::CreateChromatinAssay(data = atac_X)

h5closeAll()

cat("Done with creating Seurat object, now preprocessing...\n")

# Run single-cell preprocessing pipeline
# Some parameters are hard-coded, for now

prepareSeuratData_GRaNIE(seu.s, outputDir = outdir, saveSeuratObject = TRUE,
                                     file_RNA_features = file_RNA_features,
                                     assayName_RNA = "RNA", assayName_ATAC= "ATAC", 
                                     prepareData = TRUE, SCT_nDimensions = 50,  dimensionsToIgnore_LSI_ATAC = 1,
                                     pseudobulk_source = "cluster",
                                     countAggregation = xargs$preprocessing_countAggregation,
                                     clusteringAlgorithm = 3, clusterResolutions = NULL, 
                                     minClusters = xargs$preprocessing_minClusters,
                                     minCellsPerCluster = xargs$preprocessing_minCellsPerCluster,
                                     subsample_percentage = 100, subsample_n = 1,
                                     doDimPlots = TRUE,
                                     run_MultiK = FALSE,
                                     forceRerun = TRUE
) 

