library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_geneids <- args[8]
tmp_dir <- args[9]
path_motifs <- args[10]
path_p2g <- args[11]
path_out <- args[12]

format_peaks <- function(peaks){
    split_strings <- strsplit(peaks, "-")
    chromosomes <- sapply(split_strings, `[`, 1)
    positions <- sapply(split_strings, function(x) paste(x[-1], collapse = "-"))
    formated_peaks <- paste(chromosomes, positions, sep = ":")
    return(formated_peaks)
}

# Read gene table
path_geneids <- file.path(path_geneids, paste0(organism, '.csv'))
gene_table <- read.csv(path_geneids)
gids <- setNames(gene_table$id, gene_table$symbol)
gsym <- setNames(gene_table$symbol, gene_table$id)

# Read data
print('Open object')
indata <- H5Fopen(path_data)

# RNA
rna_data <- as.data.frame(indata$mod$rna$X)
colnames(rna_data) <- indata$obs$`_index`
rna_data[, 'ENSEMBL'] <- unname(gids[indata$mod$rna$var$`_index`])

### ATAC
atac_data <- as.data.frame(indata$mod$atac$X)
colnames(atac_data) <- indata$obs$`_index`
atac_data[, 'peakID'] <- format_peaks(indata$mod$atac$var$`_index`)
h5closeAll()

# Read p2g results
p2g <- read.csv(path_p2g)
p2g$cre <- format_peaks(p2g$cre)
p2g$gene <- unname(gids[p2g$gene])

# Subset data by p2g
rna_data <- dplyr::filter(rna_data, ENSEMBL %in% p2g$gene)
atac_data <- dplyr::filter(atac_data, peakID %in% p2g$cre)
                        
# Init GRN object
GRN <- GRaNIE::initializeGRN(
    objectMetadata = NULL,
    outputFolder = tmp_dir,
    genomeAssembly = organism
)

# Add and annotate genes
GRN = GRaNIE::addData(
    GRN,
    counts_peaks = atac_data,
    normalization_peaks = "none",
    idColumn_peaks = 'peakID',
    counts_rna = rna_data,
    normalization_rna = "none",
    idColumn_RNA = 'ENSEMBL',
    sampleMetadata = NULL,
    genomeAnnotationSource = 'AnnotationHub',
    forceRerun = FALSE
)

# Add TFBS
GRN = GRaNIE::addTFBS(
    GRN,
    motifFolder = file.path(path_motifs, organism),
    TFs = 'all',
    filesTFBSPattern = "_TFBS",
    fileEnding = ".bed.gz",
    forceRerun = TRUE
)

# Overlap with peaks
GRN = GRaNIE::overlapPeaksAndTFBS(
    GRN,
    nCores = parallel::detectCores(),
    forceRerun = TRUE
)

# Extact and process tfb
tfb <- GRN@data$TF$TF_peak_overlap
sparse <- summary(tfb)
tfb <- data.frame(
    cre = rownames(tfb)[sparse$i],
    motif = colnames(tfb)[sparse$j]
)
motif2ens <- setNames(GRN@annotation$TFs$TF.ENSEMBL, GRN@annotation$TFs$TF.ID)
tfb <- dplyr::mutate(tfb,
    tf = gsym[motif2ens[as.character(tfb$motif)]],
    cre = stringr::str_replace_all(cre, ':', '-')
)
tfb <- dplyr::summarize(tfb, score = dplyr::n(), .by=c(cre, tf))
tfb <- dplyr::arrange(tfb, cre, desc(score))

# Write
write.csv(x = tfb, file = path_out, row.names=FALSE)
