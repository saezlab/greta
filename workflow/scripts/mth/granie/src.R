library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_geneids <- args[7]
tmp_dir <- args[8]
ext <- as.numeric(args[9])
path_motifs <- args[10]
thr_fdr <- as.numeric(args[11])
path_out <- args[12]
nCores <- as.numeric(args[13])

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

# Add tfb
GRN = GRaNIE::addConnections_TF_peak(
    GRN,
    plotDiagnosticPlots = FALSE,
    connectionTypes = c("expression"),
    corMethod = "pearson",
    forceRerun = TRUE
)

# Add p2g
GRN = GRaNIE::addConnections_peak_gene(
    GRN,
    corMethod = "pearson",
    promoterRange = round(ext / 2),
    TADs = NULL,
    nCores = nCores,
    plotDiagnosticPlots = FALSE,
    plotGeneTypes = list(c("all")),
    forceRerun = TRUE
)

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
