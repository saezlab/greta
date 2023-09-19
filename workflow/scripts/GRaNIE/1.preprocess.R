#!/usr/local/bin/Rscript

library(argparse)
library(rhdf5)
library(Signac)
library(Matrix)
library(Seurat)

parser <- ArgumentParser(description= 'GRaNIE preprocessing')
parser$add_argument('--input', '-i', help= 'Input file that contains both raw RNA and raw ATAC counts (h5 file)',  required= TRUE)
parser$add_argument('--output', '-o', help= 'A list of 3 exactly 3 output files: 1. RNA counts, 2. ATAC counts, 3. metadata', nargs= 3, required= TRUE)
parser$add_argument('--param', '-p', help= 'An optional integer parameter', default= 1, type= 'integer')
xargs <- parser$parse_args()

str(xargs)

# ./1.preprocess.R --input "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/mdata.h5mu" --output "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/counts_RNA_raw.tsv.gz" "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/counts_ATAC_raw.tsv.gz" "/g/zaugg/carnold/Projects/GRETA/src/greta_benchmark/resources/neurips2021/small/GRaNIE/metadata.tsv.gz"

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
muo_data <- Seurat::CreateSeuratObject(rna_X)

### ATAC
atac_indices <- indata$mod$atac$layers$counts$indices
atac_indptr <- indata$mod$atac$layers$counts$indptr
atac_data <- as.numeric(indata$mod$atac$layers$counts$data)
barcodes <- indata$mod$atac$obs$`_index`
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE, dimnames = list(peaks, barcodes))
muo_data[['ATAC']] <- Signac::CreateChromatinAssay(data = atac_X)

h5closeAll()

cat("Done\n")

# TODO: Needed?
## Annotate peaks
# print('Add peak annotations')
# if (organism == 'human'){
#     gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Hsapiens.v86)
# } else {
#     gene.ranges <- Signac::GetGRangesFromEnsDb(ensdb = EnsDb.Mmusculus.v79)
# }
# seqlevelsStyle(gene.ranges) <- "UCSC"
# Signac::Annotation(muo_data[['peaks']]) <- gene.ranges

# Save

# atac_X %>%
#     as.matrix() %>%
#     readr::write_tsv(xargs$output[1])
# 
# 
# all_peaks <- row.names(exprs(input_cds))
# write.csv(x = all_peaks, file = file.path(path_all_peaks))
# write.csv(x = conns, file = file.path(path_connections))
