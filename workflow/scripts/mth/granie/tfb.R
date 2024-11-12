library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_geneids <- args[7]
tmp_dir <- args[8]
path_motifs <- args[9]
path_p2g <- args[10]
path_out <- args[11]
nCores <- as.numeric(args[12])

format_peaks <- function(peaks){
    split_strings <- strsplit(peaks, "-")
    chromosomes <- sapply(split_strings, `[`, 1)
    positions <- sapply(split_strings, function(x) paste(x[-1], collapse = "-"))
    formated_peaks <- paste(chromosomes, positions, sep = ":")
    return(formated_peaks)
}

# Read gene table
gene_table <- read.csv(path_geneids)
gids <- setNames(gene_table$id, gene_table$symbol)
gsym <- setNames(gene_table$symbol, gene_table$id)
organism <- sub('^dbs/([^/]+)/.*$', '\\1', path_geneids)

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')

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
if (nrow(p2g) == 0){
    tfb <- data.frame(cre=character(), tf=character(), score=numeric())
    write.csv(x = tfb, file = path_out, row.names=FALSE)
    dir.create(tmp_dir)
    quit(save="no")
}
p2g$cre <- format_peaks(p2g$cre)
p2g$gene <- unname(gids[p2g$gene])
                        
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
    motifFolder = path_motifs[grep(organism, path_motifs)],
    TFs = 'all',
    filesTFBSPattern = "_TFBS",
    fileEnding = ".bed.gz",
    forceRerun = TRUE
)

# Overlap with peaks
GRN = GRaNIE::overlapPeaksAndTFBS(
    GRN,
    nCores = nCores,
    forceRerun = TRUE
)

# Calculate tfb corr
GRN@connections$TF_peaks$`0` = GRaNIE:::.computeTF_peak.fdr(
    GRN,
    perm = 0,
    connectionTypes = c("expression"),
    corMethod = 'pearson',
    removeNegativeCorrelation = FALSE,
    maxFDRToStore = 0.3,
    useGCCorrection = FALSE,
    percBackground_size = 75,
    threshold_percentage = 0.05,
    percBackground_resample = TRUE,
    plotDetails = FALSE
)
permIndex = as.character(0)
GRN@connections$TF_peaks$`0`$main = GRaNIE:::.optimizeSpaceGRN(stats::na.omit(GRN@connections$TF_peaks$`0`$main))
GRN@connections$TF_peaks$`1` <- GRN@connections$TF_peaks$`0`

# Extact and process tfb
tfb <- GRN@connections$TF_peaks$`0`$main
min_pval <- min(GRN@connections$TF_peaks$`0`$main$TF_peak.fdr[GRN@connections$TF_peaks$`0`$main$TF_peak.fdr != 0])
tfb <- dplyr::mutate(tfb,
    score = dplyr::if_else(tfb$TF_peak.fdr == 0, min_pval, tfb$TF_peak.fdr)
)

motif2ens <- setNames(GRN@annotation$TFs$TF.ENSEMBL, GRN@annotation$TFs$TF.ID)
tfb <- dplyr::mutate(tfb,
    tf = gsym[motif2ens[as.character(tfb$TF.ID)]],
    cre = stringr::str_replace_all(tfb$peak.ID, ':', '-'),
    score = -log10(score)
)
tfb <- dplyr::select(tfb, tf, cre, score)
tfb <- dplyr::summarize(tfb, score = mean(score), .by=c(cre, tf))

clean_peaks <- function(peaks) {
  parts <- strsplit(peaks, "[:-]")[[1]]
  # Ensure that start and end are not in scientific notation
  chromosome <- parts[1]
  start <- format(as.numeric(parts[2]), scientific = FALSE)
  end <- format(as.numeric(parts[3]), scientific = FALSE)
  fpeaks <- paste0(chromosome, "-", start, "-", end)
  return(fpeaks)
}
tfb$cre <- sapply(tfb$cre, clean_peaks)

# Write
write.csv(x = tfb, file = path_out, row.names=FALSE)
