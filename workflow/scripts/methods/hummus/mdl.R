library(tidyverse)
library(rhdf5)
library(HuMMuS)


# Parse args
args <- commandArgs(trailingOnly = F)
hummus_object_f <- args[6]
p2g_f <- args[7]
tfb_f <- args[8]
path_out <- args[9]
multilayer_f <- args[10]


#load {pre}.hummus_object, containing formatted preprocessed data
hummus <- readRDS(hummus_object_f)

#load p2g and tfb bipartites, to be added to this hummus_object
p2g <- read.csv2(p2g_f)
p2g <- p2g[, c("gene", "cre", "score")]

tfb <- read.csv2(tfb_f)
tfb <- tfb[, c("tf", "cre", "score")]


hummus@multilayer@bipartite['atac_rna'] <- new("bipartite",
                           "network" = p2g,
                           "multiplex_left" = "RNA",
                           "multiplex_right" = "peaks")

hummus@multilayer@bipartite['tf_peak'] <- new("bipartite",
                           "network" = tfb,
                           "multiplex_left" = "TF",
                           "multiplex_right" = "peaks")


# Save the multilayer, necessary since multixrank uses locally saved files
save_multilayer(
  hummus,
  folder_name = multilayer_f)


# Run HuMMuS enhancers search
tfb <- define_grn(
  hummus,
#  gene_list = list("ATF2"),
  multilayer_f = multilayer_f,
  njobs = 1
  )

#get only gene, peaks and score columns
tfb <- tfb[, c("peak", "tf", "score")]
colnames(tfb) <- c("cre", "tf",  "score")

# Write
write.csv(
    x = tfb,
    file = path_out,
    row.names = FALSE,
    quote = FALSE)
