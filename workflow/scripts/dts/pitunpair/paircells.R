library(doParallel)
library(FigR)
library(BSgenome.Hsapiens.UCSC.hg38)
library(SingleCellExperiment)
options("optmatch_max_problem_size" = Inf)
optmatch::setMaxProblemSize(size = Inf)


# Parse args
args <- commandArgs(trailingOnly = F)
path_inp <- args[6]
path_ann <- args[7]
path_out <- args[8]

# Load Data
x <- read.csv(path_inp)
x <- x[!duplicated(x[[1]]), ]
rownames(x) <- x[[1]]
x <- as.matrix(x[ , -1])
x_glue_rna <- x[grepl("^rna_", rownames(x)), ]
x_glue_atac <- x[grepl("^atac_", rownames(x)), ]

# Pair with FigR
pairing <- pairCells(
    ATAC = x_glue_atac,
    RNA = x_glue_rna,
    min_subgraph_size = 100,
    keepUnique = TRUE
)

# Filter paired object
pairing <- pairing[order(pairing$dist, decreasing = FALSE), ]
pairing <- pairing[!duplicated(pairing$ATAC) & !duplicated(pairing$RNA) ,]

# Merge ctype info
ann <- read.csv(path_ann, row.names=1)
ann$X <- paste0("rna_", sapply(strsplit(as.character(ann$X), "-"), `[`, 1))
ann <- ann[!duplicated(ann[[1]]), ]
pairing <- merge(pairing, ann, by.x='RNA', by.y='X')
pairing['batch'] <- 'smpl'
pairing <- pairing[, c('ATAC', 'RNA', 'batch', 'celltype', 'dist')]
pairing$RNA <- gsub("rna_", "smpl_", as.character(pairing$RNA))
pairing$ATAC <- gsub("atac_", "smpl_", as.character(pairing$ATAC))
rownames(pairing) <- pairing$ATAC

# Write
write.csv(pairing, path_out, row.names = FALSE)
