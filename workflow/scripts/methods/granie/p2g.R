library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_geneids <- args[8]
tmp_dir <- args[9]
ext <- as.numeric(args[10])
path_out <- args[11]


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

# P2G 
genomeAssembly = GRN@config$parameters$genomeAssembly
consensusPeaks =  dplyr::filter(GRN@data$peaks$counts_metadata, !.data$isFiltered)
peak.TADs.df = NULL
gene.types = stats::na.omit(unique(GRN@annotation$genes$gene.type))
overlaps.sub.df = GRaNIE:::.calculatePeakGeneOverlaps(
    GRN, allPeaks = consensusPeaks,
    peak.TADs.df,
    neighborhoodSize = ext,
    genomeAssembly = genomeAssembly,
    gene.types = gene.types,
    overlapTypeGene = 'TSS'
)

# Process
p2g = dplyr::mutate(overlaps.sub.df, gene.ENSEMBL = gsub("\\..+", "", .data$gene.ENSEMBL, perl = TRUE))
p2g <- dplyr::select(p2g, peak.ID, gene.ENSEMBL, peak_gene.distance)
p2g <- dplyr::mutate(
    p2g,
    gene = gsym[as.character(gene.ENSEMBL)],
    cre = stringr::str_replace_all(peak.ID, ':', '-'),
    score = -peak_gene.distance
)
p2g <- dplyr::summarise(p2g, score=max(score), .by=c(cre, gene))
p2g <- dplyr::arrange(p2g, desc(score))

# Write
write.csv(x = p2g, file = path_out, row.names=FALSE)
