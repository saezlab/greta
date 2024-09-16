library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
hummus_object_f <- args[6]
organism <- args[7]
extend <- as.numeric(args[8])
n_cores <- args[9]
n_cores <- if (n_cores == 'NULL') 1 else as.numeric(n_cores)
path_out <- args[10]
multilayer_f <- args[11]
print(extend)
print(multilayer_f)
print("n_cores of type")
print(class(n_cores))

# Set genome
if (organism == 'hg38'){
    genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    specie <- 'human'
} else if (organism == 'mm10'){
    genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    specie <- 'mouse'
}

#load {pre}.hummus_object, containing formatted preprocessed data
hummus <- readRDS(hummus_object_f)
hummus@motifs_db <- get_tf2motifs(species = specie) # TF motif DB # by default human motifs

print("HuMMuS object:")
print(hummus)

# Add bipartites from hummus method
hummus <- bipartite_tfs2peaks(
  hummus_object = hummus,
  tf_expr_assay = NULL, # use to filter TF on only expressed TFs,
                         # if NULL, all TFs with motifs are used
  peak_assay = "peaks",
  tf_multiplex_name = "TF",
  genome = genome,
  )


hummus@assays$peaks
hummus <- bipartite_peaks2genes(
  hummus_object = hummus,
  gene_assay = "RNA",
  peak_assay = "peaks",
  store_network = FALSE,
  upstream = extend,
  downstream = extend
)

# Save the multilayer, necessary since multixrank uses locally saved files
save_multilayer(
  hummus,
  folder_name = multilayer_f
)

print("Calculating peak-to-gene scores")
# Run HuMMuS enhancers search
enhancers <- define_enhancers(
  hummus,
  multilayer_f = multilayer_f,
  njobs = n_cores
)

#get only gene, peaks and score columns
enhancers <- enhancers[, c("gene", "peak", "score")]
colnames(enhancers) <- c("gene", "cre", "score")

print("Saving peak-to-gene results")
# Write
write.csv(
  x = enhancers,
  file = path_out,
  row.names = FALSE,
  quote = FALSE
)
