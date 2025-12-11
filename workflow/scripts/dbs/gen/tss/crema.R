library(EnsDb.Hsapiens.v86)

# Add arguments
args <- commandArgs(trailingOnly = F)
path_out <- args[6]

genebody_coords <- keepStandardChromosomes(
    ensembldb::genes(EnsDb.Hsapiens.v86),
    species="Homo_sapiens",
    pruning.mode = 'coarse'
)
seqlevels(genebody_coords) <- paste0("chr", seqlevels(genebody_coords))
genebody_coords <- genebody_coords[sapply(genebody_coords$gene_name, nchar) > 0]
temp_ind <- table(genebody_coords$gene_name)
temp_ind <- names(temp_ind)[temp_ind > 1]
genebody_coords_list <- split(genebody_coords, f = genebody_coords$gene_name)
temp_func <- function(x){
    if ("protein_coding" %in% x$gene_biotype){ return(x[x$gene_biotype == "protein_coding"])
    }else{ return(x) }
}
genebody_coords_list <- c(
    lapply(genebody_coords_list[temp_ind], temp_func),
    as.list(genebody_coords_list[setdiff(names(genebody_coords_list), temp_ind)])
)
genebody_coords <- unlist(as(genebody_coords_list, "GRangesList"))
genebody_coords <- sort(genebody_coords)

df_bed <- data.frame(
  chrom = as.character(seqnames(genebody_coords)),
  start = start(genebody_coords),
  end   = end(genebody_coords),
  name  = genebody_coords$gene_name,
  stringsAsFactors = FALSE
)

# Write
write.table(x = df_bed, file = gzfile(path_out), sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
