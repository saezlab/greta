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
print(multilayer_f)

#load {pre}.hummus_object, containing formatted preprocessed data
hummus <- readRDS(hummus_object_f)

#load p2g and tfb bipartites, to be added to this hummus_object
p2g <- read.csv(p2g_f)
p2g <- p2g[, c("gene", "cre", "score")]

tfb <- read.csv(tfb_f)
tfb <- tfb[, c("tf", "cre", "score")]


hummus@multilayer@bipartites['atac_rna'] <- new("bipartite",
                           "network" = p2g,
                           "multiplex_left" = "RNA",
                           "multiplex_right" = "peaks")

hummus@multilayer@bipartites['tf_peak'] <- new("bipartite",
                           "network" = tfb,
                           "multiplex_left" = "TF",
                           "multiplex_right" = "peaks")


# Save the multilayer, necessary since multixrank uses locally saved files
save_multilayer(
  hummus,
  folder_name = multilayer_f)


# Run HuMMuS enhancers search
grn <- define_grn(
  hummus,
#  gene_list = list("ATF2"),
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
