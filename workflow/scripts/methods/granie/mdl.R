library(rhdf5)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_geneids <- args[8]
tmp_dir <- args[9]
path_p2g <- args[10]
path_tfb <- args[11]
thr_fdr <- as.numeric(args[12])
path_out <- args[13]


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

# Subset data by p2g and tfb
p2g <- dplyr::filter(p2g, cre %in% tfb$cre)
rna_data <- dplyr::filter(rna_data, (ENSEMBL %in% p2g$gene) | (ENSEMBL %in% tfb$tf))
atac_data <- dplyr::filter(atac_data, peakID %in% tfb$cre)

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
row_idx <- as.factor(tfb$cre)
col_idx <- as.factor(tfb$tf)
GRN@data$TFs$TF_peak_overlap <- Matrix::sparseMatrix(
  i = as.integer(row_idx),
  j = as.integer(col_idx),
  x = 1,
  dims = c(length(levels(row_idx)), length(levels(col_idx))),
  dimnames = list(levels(row_idx), levels(col_idx))
)[rownames(GRN@data$peaks$counts), ]
GRN@data$TFs$TF_peak_overlap <- cbind(GRN@data$TFs$TF_peak_overlap, isFiltered = 0)

# Compute tf-cre corrs
GRN@connections$TF_peaks$`0`  = GRaNIE:::.computeTF_peak.fdr(
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

# Compute p2g corrs
model_p2g <- function(GRN, p2g, nCores=1, chunksize=50000){
    # Format p2g
    overlaps.sub.filt.df <- dplyr::rename(
        p2g,
        peak.ID = cre,
        gene.ENSEMBL = gene
    )
    countsPeaks.clean = GRaNIE::getCounts(GRN, type = "peaks",  permuted = FALSE, includeIDColumn = FALSE)
    countsRNA.clean = GRaNIE::getCounts(GRN, type = "rna", permuted = FALSE, includeIDColumn = FALSE)
    map_peaks = match(overlaps.sub.filt.df$peak.ID,  rownames(countsPeaks.clean))
    map_rna  = match(overlaps.sub.filt.df$gene.ENSEMBL, rownames(countsRNA.clean))
    maxRow = nrow(overlaps.sub.filt.df)
    startIndexMax = ceiling(maxRow / chunksize)
    
    # Run correlation
    res.l = GRaNIE:::.execInParallelGen(
        nCores,
        returnAsList = TRUE,
        listNames = NULL,
        iteration = 0,#0:startIndexMax,
        verbose = TRUE,
        functionName = GRaNIE:::.correlateDataWrapper,
        chunksize = chunksize,
        maxRow = maxRow,
        counts1 = countsPeaks.clean,
        counts2 = countsRNA.clean,
        map1 = map_peaks,
        map2 = map_rna, 
        corMethod = 'pearson'
    )
    res.m  = do.call(rbind, res.l)
    
    # Convert to df and format
    res.df = tibble::as_tibble(res.m)
    res.df = dplyr::mutate(
        res.df,
        peak.ID = rownames(countsPeaks.clean)[map_peaks],
        gene.ENSEMBL = rownames(countsRNA.clean)[map_rna]
    )
    res.df = dplyr::filter(res.df, !is.na(.data$gene.ENSEMBL))
    res.df = dplyr::left_join(res.df, overlaps.sub.filt.df, by = c("gene.ENSEMBL", "peak.ID"), multiple = "all")
    selectColumns = c("peak.ID", "gene.ENSEMBL", "r", "p.raw")
    res.df = dplyr::select(res.df, tidyselect::all_of(selectColumns))
    res.df = dplyr::mutate(
        res.df,
        peak.ID = as.factor(.data$peak.ID),
        gene.ENSEMBL = as.factor(.data$gene.ENSEMBL), 
    )
    res.df = dplyr::rename(res.df, peak_gene.r = "r", peak_gene.p_raw = "p.raw")
    res.df = dplyr::arrange(res.df, peak.ID)
    return(res.df)
}
GRN@connections$peak_genes$`0` <- model_p2g(GRN, p2g, nCores=6, chunksize=50000)
GRN@connections$peak_genes$`1` <- GRN@connections$peak_genes$`0`

# Format final grn
GRN = GRaNIE::filterGRNAndConnectGenes(
    GRN,
    TF_peak.fdr.threshold = thr_fdr,
    peak_gene.fdr.threshold = thr_fdr,
    peak_gene.fdr.method = "BH",
    gene.types = c("protein_coding", "lincRNA"),
    allowMissingTFs = FALSE,
    allowMissingGenes = FALSE,
    forceRerun = TRUE
)
mdl <- dplyr::select(GRN@connections$all.filtered$`0`, TF.ENSEMBL, TF_peak.r, peak.ID, peak_gene.r, gene.ENSEMBL, TF_peak.fdr, peak_gene.p_adj)
mdl <- dplyr::mutate(mdl, score = TF_peak.r * peak_gene.r)
mdl <- dplyr::mutate(mdl, source=gsym[as.character(TF.ENSEMBL)], target=gsym[as.character(gene.ENSEMBL)])
mdl <- dplyr::mutate(dplyr::rowwise(mdl), pval = mean(c(TF_peak.fdr, peak_gene.p_adj)))
mdl <- dplyr::select(dplyr::ungroup(mdl), source, target, score, pval)
mdl <- dplyr::summarize(mdl, score = mean(score), pval=mean(pval), .by=c(source, target))
mdl <- dplyr::arrange(mdl, desc(abs(score)))

# Write
write.csv(x = mdl, file = path_out, row.names=FALSE)
