library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_geneids <- args[7]
tmp_dir <- args[8]
ext <- as.numeric(args[9])
path_out <- args[10]
nCores <- as.numeric(args[11])

print(args)
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

# Identify possible peaks and genes
genomeAssembly = GRN@config$parameters$genomeAssembly
consensusPeaks =  dplyr::filter(GRN@data$peaks$counts_metadata, !.data$isFiltered)
peak.TADs.df = NULL
gene.types = stats::na.omit(unique(GRN@annotation$genes$gene.type))
overlaps.sub.filt.df = GRaNIE:::.calculatePeakGeneOverlaps(
    GRN, allPeaks = consensusPeaks,
    peak.TADs.df,
    neighborhoodSize = round(ext / 2),
    genomeAssembly = genomeAssembly,
    gene.types = gene.types,
    overlapTypeGene = 'TSS'
)
overlaps.sub.filt.df = dplyr::mutate(overlaps.sub.filt.df, gene.ENSEMBL = gsub("\\..+", "", .data$gene.ENSEMBL, perl = TRUE))
overlaps.sub.filt.df <- dplyr::select(overlaps.sub.filt.df, peak.ID, gene.ENSEMBL, peak_gene.distance)

# Compute p2g corrs
model_p2g <- function(GRN, overlaps.sub.filt.df, nCores=nCores, chunksize=50000){
    # Format p2g
    countsPeaks.clean = GRaNIE::getCounts(GRN, type = "peaks",  permuted = FALSE, includeIDColumn = FALSE)
    countsRNA.clean = GRaNIE::getCounts(GRN, type = "rna", permuted = FALSE, includeIDColumn = FALSE)
    map_peaks = match(overlaps.sub.filt.df$peak.ID,  rownames(countsPeaks.clean))
    map_rna  = match(overlaps.sub.filt.df$gene.ENSEMBL, rownames(countsRNA.clean))
    maxRow = nrow(overlaps.sub.filt.df)
    cat('Running p2g for', length(unique(map_peaks)), 'peaks and', length(unique(map_rna)), 'genes.\n')

    # Run correlation
    res.l = GRaNIE:::.execInParallelGen(
        nCores,
        returnAsList = TRUE,
        listNames = NULL,
        iteration = 0,
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
p2g <- model_p2g(GRN, overlaps.sub.filt.df, nCores=nCores, chunksize=1000000)

# Process
p2g <- dplyr::mutate(p2g, pval=peak_gene.p_raw)
p2g <- dplyr::mutate(
    p2g,
    gene = gsym[as.character(gene.ENSEMBL)],
    cre = stringr::str_replace_all(peak.ID, ':', '-'),
    score = peak_gene.r
)
p2g <- dplyr::select(p2g, cre, gene, score, pval)
p2g <- dplyr::arrange(p2g, desc(score))

clean_peaks <- function(peaks) {
  parts <- strsplit(peaks, "[:-]")[[1]]
  # Ensure that start and end are not in scientific notation
  chromosome <- parts[1]
  start <- format(as.numeric(parts[2]), scientific = FALSE)
  end <- format(as.numeric(parts[3]), scientific = FALSE)
  fpeaks <- paste0(chromosome, "-", start, "-", end)
  return(fpeaks)
}
p2g$cre <- sapply(p2g$cre, clean_peaks)

# Write
write.csv(x = p2g, file = path_out, row.names=FALSE)
