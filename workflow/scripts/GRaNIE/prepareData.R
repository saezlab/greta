library(cicero)
library(monocle3)
library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_all_peaks <- args[8]
path_connections <- args[9]

# Read genome
if (organism == 'human'){
    genome <- read.table('resources/genome_sizes/human.txt')
} else {
    genome <- read.table('resources/genome_sizes/mouse.txt')
}

# Process mudata
indata <- H5Fopen(path_data)
indices <- indata$mod$atac$layers$counts$indices
indptr <- indata$mod$atac$layers$counts$indptr
data <- as.numeric(indata$mod$atac$layers$counts$data)
barcodes <- indata$mod$atac$obs$`_index`
peaks <- indata$mod$atac$var$`_index`
h5closeAll()

# Build sparse matrix and binarize
indata <- Matrix::sparseMatrix(i=indices, p=indptr, x=data, index1 = FALSE)
indata@x[indata@x > 0] <- 1

# Format cell info
cellinfo <- data.frame(row.names=barcodes, cells=barcodes)

# Format peak info
peakinfo <- data.frame(row.names=peaks, site_name=peaks)
peakinfo <- tidyr::separate(data = peakinfo, col = 'site_name', into = c("chr", "bp1", "bp2"), sep = "-", remove=FALSE)

# Add names
row.names(indata) <- row.names(peakinfo)
colnames(indata) <- row.names(cellinfo)

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)

# Data preprocessing
set.seed(2017)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# Build cicero cds
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Run the main function
conns <- run_cicero(cicero_cds, genome) # Takes a few minutes to run

# Save
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(path_all_peaks))
write.csv(x = conns, file = file.path(path_connections))
