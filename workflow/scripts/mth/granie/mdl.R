library(rhdf5)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_geneids <- args[7]
tmp_dir <- args[8]
path_p2g <- args[9]
path_tfb <- args[10]
thr_fdr <- as.numeric(args[11])
path_out <- args[12]
nCores <- as.numeric(args[13])
cat("N cores: ", nCores, '\n')

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

# Read p2g and tfb results
p2g <- read.csv(path_p2g)
tfb <- read.csv(path_tfb)
if ((nrow(p2g) == 0) | (nrow(tfb) == 0)){
    mdl <- data.frame(source=character(), target=character(), score=numeric(), pval=numeric())
    write.csv(x = mdl, file = path_out, row.names=FALSE)
    dir.create(tmp_dir)
    quit(save="no")
}
tfb$cre <- format_peaks(tfb$cre)
tfb$tf <- unname(gids[tfb$tf])
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

# Add tf annotation
tf_ids <- unique(tfb$tf)
GRN@annotation$TFs <- data.frame(
    TF.ID = tf_ids,
    TF.name = gsym[tf_ids],
    TF.ENSEMBL = tf_ids,
    row.names = NULL
)

# Add tfb predictions
clean_peaks <- function(peaks) {
  parts <- strsplit(peaks, "[:-]")[[1]]
  # Ensure that start and end are not in scientific notation
  chromosome <- parts[1]
  start <- format(as.numeric(parts[2]), scientific = FALSE)
  end <- format(as.numeric(parts[3]), scientific = FALSE)
  fpeaks <- paste0(chromosome, ":", start, "-", end)
  return(fpeaks)
}
row_idx <- as.factor(tfb$cre)
col_idx <- as.factor(tfb$tf)
tmp_filtered <- Matrix::sparseMatrix(
  i = as.integer(row_idx),
  j = as.integer(col_idx),
  x = 1,
  dims = c(length(levels(row_idx)), length(levels(col_idx))),
  dimnames = list(levels(row_idx), levels(col_idx))
)
tmp_filtered <- cbind(tmp_filtered, isFiltered = 0)
msk <- !(rownames(GRN@data$peaks$counts) %in% rownames(tmp_filtered))
empty_cres <- rownames(GRN@data$peaks$counts)[msk]
tmp_unfiltered <- Matrix::Matrix(0, nrow = length(empty_cres), ncol = length(colnames(tmp_filtered))-1, sparse = TRUE)
rownames(tmp_unfiltered) <- empty_cres
tmp_unfiltered <- cbind(tmp_unfiltered, isFiltered = 0)
colnames(tmp_unfiltered) <- colnames(tmp_filtered)
GRN@data$TFs$TF_peak_overlap <- rbind(tmp_filtered, tmp_unfiltered)
GRN@data$TFs$TF_peak_overlap <- GRN@data$TFs$TF_peak_overlap[rownames(GRN@data$peaks$counts), ]

# Compute tf-cre corrs
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

# Add p2g
# Add peak_gene.r=0.5 and peak_gene.p_raw = 0.01 to bypass filters
p2g_method <- stringr::str_extract(path_p2g, "[^/.]+(?=\\.p2g)")

if (p2g_method != 'granie'){
    p2g <- dplyr::mutate(
        p2g,
        peak_gene.r=0.5,
        peak_gene.p_raw=0.001,
        peak.ID=p2g$cre,
        gene.ENSEMBL=p2g$gene
    )
} else {
    p2g <- dplyr::mutate(
        p2g,
        peak_gene.r=p2g$score,
        peak_gene.p_raw=p2g$pval,
        peak.ID=p2g$cre,
        gene.ENSEMBL=p2g$gene
    )
}
p2g <- dplyr::select(p2g, peak.ID, gene.ENSEMBL, peak_gene.r, peak_gene.p_raw)
GRN@connections$peak_genes$`0` <- p2g
GRN@connections$peak_genes$`1` <- GRN@connections$peak_genes$`0`

# Format final grn
GRN = GRaNIE::filterGRNAndConnectGenes(
    GRN,
    TF_peak.fdr.threshold = thr_fdr,
    peak_gene.fdr.threshold = thr_fdr,
    peak_gene.fdr.method = "BH",
    gene.types = c("all"),
    allowMissingTFs = FALSE,
    allowMissingGenes = FALSE,
    forceRerun = TRUE
)
mdl <- dplyr::select(GRN@connections$all.filtered$`0`, TF.ENSEMBL, TF_peak.r, peak.ID, peak_gene.r, gene.ENSEMBL, TF_peak.fdr, peak_gene.p_adj)
mdl <- dplyr::mutate(dplyr::rowwise(mdl), score = sign(TF_peak.r) * mean(c(abs(TF_peak.r), abs(peak_gene.r))))
mdl <- dplyr::mutate(dplyr::ungroup(mdl), source=gsym[as.character(TF.ENSEMBL)], target=gsym[as.character(gene.ENSEMBL)])
mdl <- dplyr::mutate(dplyr::rowwise(mdl), pval = mean(c(TF_peak.fdr, peak_gene.p_adj)))
mdl <- dplyr::select(dplyr::ungroup(mdl), source, target, score, pval)
mdl <- dplyr::summarize(mdl, score = mean(score), pval=mean(pval), .by=c(source, target))
mdl <- dplyr::arrange(mdl, source, target)

# Write
write.csv(x = mdl, file = path_out, row.names=FALSE)
