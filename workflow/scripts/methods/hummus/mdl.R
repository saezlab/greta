library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
hummus_object_f <- args[6]
p2g_f <- args[7]
tfb_f <- args[8]
n_cores <- args[9]
path_out <- args[10]
multilayer_f <- args[11]
p2g <- args[12]
tfb <- args[13]
organism <- args[14]

# Set genome
if (organism == 'hg38') {
  genome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
} else if (organism == 'mm10') {
  genome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
}

#load {pre}.hummus_object, containing formatted preprocessed data
hummus <- readRDS(hummus_object_f)

#load p2g and tfb bipartites, to be added to this hummus_object
if (p2g != "hummus"){
  p2g <- read.csv(p2g_f)
  p2g <- p2g[, c("gene", "cre", "score")]
  hummus@multilayer@bipartites['atac_rna'] <- new(
    "bipartite",
    "network" = p2g,
    "multiplex_left" = "RNA",
    "multiplex_right" = "peaks"
  )
} else {
  hummus <- bipartite_peaks2genes(
    hummus_object = hummus,
    gene_assay = "RNA",
    peak_assay = "peaks",
    store_network = FALSE,
  )
}

if (tfb != "hummus"){
  tfb <- read.csv(tfb_f)
  tfb <- tfb[, c("tf", "cre", "score")]
  hummus@multilayer@bipartites['tf_peak'] <- new(
    "bipartite",
    "network" = tfb,
    "multiplex_left" = "TF",
    "multiplex_right" = "peaks")
} else {
  hummus <- bipartite_tfs2peaks(
    hummus_object = hummus,
    tf_expr_assay = "RNA", # use to filter TF on only expressed TFs,
                           # if NULL, all TFs with motifs are used
    peak_assay = "peaks",
    tf_multiplex_name = "TF",
    genome = genome,
  )
}

# Save the multilayer, necessary since multixrank uses locally saved files
save_multilayer(
  hummus,
  folder_name = multilayer_f
)


# Run HuMMuS enhancers search
grn <- define_grn(
  hummus,
  multilayer_f = multilayer_f,
  njobs = n_cores
)

#get only gene, peaks and score columns
grn <- grn[, c("gene", "tf", "score")]
colnames(grn) <- c("gene", "tf",  "score")

# Write
write.csv(
  x = grn,
  file = path_out,
  row.names = FALSE,
  quote = FALSE)
