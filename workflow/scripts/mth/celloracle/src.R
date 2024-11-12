library(cicero)
library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_genome <- args[7]
ext <- as.numeric(args[8])
path_all_peaks <- args[9]
path_connections <- args[10]

# Read genome
org <- sub('^dbs/([^/]+)/.*$', '\\1', path_genome)
path_chr_sizes <- file.path(path_genome, org, sprintf('%s.fa.sizes', org))
genome <- read.table(path_chr_sizes)

# Process mudata
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
data <- indata$mod$atac$X
barcodes <- indata$mod$atac$obs$`_index`
peaks <- indata$mod$atac$var$`_index`
h5closeAll()

# Format cell info
cellinfo <- data.frame(row.names=barcodes, cells=barcodes)

# Format peak info
peakinfo <- data.frame(row.names=peaks, site_name=peaks)
peakinfo <- tidyr::separate(data = peakinfo, col = 'site_name', into = c("chr", "bp1", "bp2"), sep = "-", remove=FALSE)

# Add names
row.names(data) <- row.names(peakinfo)
colnames(data) <- row.names(cellinfo)

# Binarize
data[data != 0] <- 1

# Make CDS
input_cds <-  suppressWarnings(
    new_cell_data_set(data,
    cell_metadata = cellinfo,
    gene_metadata = peakinfo)
)

# Data preprocessing
set.seed(2017)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(
    input_cds,
    reduction_method = 'UMAP', 
    preprocess_method = "LSI"
)
umap_coords <- reducedDims(input_cds)$UMAP

# Build "metacells"
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Run cicero
print("Starting Cicero")
print("Calculating distance_parameter value")
distance_parameters <- estimate_distance_parameter(
    input_cds,
    window=ext,
    maxit=100,
    sample_num = 100,
    distance_constraint = round(ext / 2),
    distance_parameter_convergence = 1e-22,
    genomic_coords = genome
)
mean_distance_parameter <- mean(unlist(distance_parameters))
print("Running models")
cicero_out <- generate_cicero_models(
    input_cds,
    distance_parameter = mean_distance_parameter,
    window = ext,
    genomic_coords = genome
)
print("Assembling connections")
conns <- assemble_connections(cicero_out, silent=FALSE)

# Save
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(path_all_peaks))
write.csv(x = conns, file = file.path(path_connections))
