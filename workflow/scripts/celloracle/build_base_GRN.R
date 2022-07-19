library(cicero)
library(monocle3)

max_count <-  15000
min_count <- 2000
organism <- 'human'

# Determine genome
if (organism == 'human'){
    data("human.hg19.genome")
    genome <- human.hg19.genome
} else if (organism == 'mouse'){
    data("mouse.mm9.genome")
    genome <- mouse.mm9.genome
}

# Make CDS
input_cds <-  suppressWarnings(new_cell_data_set(indata,
cell_metadata = cellinfo,
gene_metadata = peakinfo))
input_cds <- monocle3::detect_genes(input_cds)

# Ensure there are no peaks included with zero reads
input_cds <- input_cds[Matrix::rowSums(exprs(input_cds)) != 0,]

# Visualize peak_count_per_cell
plt <- hist(Matrix::colSums(exprs(input_cds)))

# Filter by peak_count
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) >= min_count] 
input_cds <- input_cds[,Matrix::colSums(exprs(input_cds)) <= max_count]

# Data preprocessing
set.seed(2017)
input_cds <- detect_genes(input_cds)
input_cds <- estimate_size_factors(input_cds)
input_cds <- preprocess_cds(input_cds, method = "LSI")

# Dimensional reduction with umap
input_cds <- reduce_dimension(input_cds, reduction_method = 'UMAP', 
                              preprocess_method = "LSI")
umap_coords <- reducedDims(input_cds)$UMAP

# Build cicero cds
cicero_cds <- make_cicero_cds(input_cds, reduced_coordinates = umap_coords)

# Run the main function
conns <- run_cicero(cicero_cds, chromosome_length) # Takes a few minutes to run

# Save
all_peaks <- row.names(exprs(input_cds))
write.csv(x = all_peaks, file = file.path(output_all_peaks))
write.csv(x = conns, file = file.path(output_connections))

