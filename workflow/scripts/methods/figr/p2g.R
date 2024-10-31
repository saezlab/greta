library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(Matrix)
library(chromVAR)
library(doParallel)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
ext <- as.numeric(args[8])
thr_pval <- as.numeric(args[9])
ncres <- as.numeric(args[10])
nCores <- as.numeric(args[11])
path_out <- args[12]

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`

# Read normalized ATAC and raw ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`

ATAC.se <- Matrix::sparseMatrix(
    i=indata$mod$atac$layers$counts$indices,
    p=indata$mod$atac$layers$counts$indptr,
    x=as.numeric(indata$mod$atac$layers$counts$data),
    index1 = FALSE
)
colnames(ATAC.se) <- indata$obs$`_index`
rownames(ATAC.se) <- indata$mod$atac$var$`_index`
h5closeAll()

# Transform atac to sme Object
peaks <- strsplit(rownames(ATAC.se), "-")
peak_ranges <- GenomicRanges::GRanges(
    seqnames = sapply(peaks, "[[", 1),
    ranges = IRanges::IRanges(start = as.numeric(sapply(peaks, "[[", 2)), end = as.numeric(sapply(peaks, "[[", 3)))
)
ATAC.se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.se), rowRanges = peak_ranges)

# Run p2g

runGenePeakcorr <- function(ATACmat,
                            ATAC.se, # SummarizedExperiment object of scATAC data
                            RNAmat, # Paired normalized scRNA-seq data, with gene names as rownames
                            genome, # Must be one of "hg19", "mm10", or "hg38"
                            geneList=NULL, # 2 or more valid gene symbols (if only running on subset of genes)
                            windowPadSize=50000, # base pairs padded on either side of gene TSS
                            normalizeATACmat=TRUE, # Whether or not to normalize scATAC counts (default is yes, assumes raw counts)
                            nCores=4, # Number of cores if parallelization support
                            keepPosCorOnly=TRUE,
                            keepMultiMappingPeaks=FALSE,
                            n_bg=100, # Number of background peaks to use
                            p.cut=NULL # Optional, if specified, will only return sig peak-gene hits
) {

  stopifnot(inherits(ATAC.se,"RangedSummarizedExperiment"))
  stopifnot(inherits(RNAmat,c("Matrix","matrix")))

  if(!all.equal(ncol(ATAC.se),ncol(RNAmat)))
    stop("Input ATAC and RNA objects must have same number of cells")

  message("Assuming paired scATAC/scRNA-seq data ..")

  peakRanges.OG <- granges(ATAC.se) # Peak ranges in reference input SE (pre-filtering)

  # Function needs rownames for both matrices or gives error
  rownames(ATAC.se) <- paste0("Peak",1:nrow(ATAC.se))
  rownames(ATACmat) <- rownames(ATAC.se)

  if(is.null(rownames(RNAmat)))
    stop("RNA matrix must have gene names as rownames")

  # Check for peaks/genes with 0 accessibility/expression

  if(any(Matrix::rowSums(assay(ATAC.se))==0)){
    message("Peaks with 0 accessibility across cells exist ..")
    message("Removing these peaks prior to running correlations ..")
    peaksToKeep <- Matrix::rowSums(assay(ATAC.se))!=0
    ATAC.se <- ATAC.se[peaksToKeep,] # Subset ranges
    ATACmat <- ATACmat[peaksToKeep,]
    message("Important: peak indices in returned gene-peak maps are relative to original input SE")
  }


  peakRanges <- granges(ATAC.se) # Peak ranges

  if(any(Matrix::rowSums(RNAmat)==0)){
    message("Genes with 0 expression across cells exist ..")
    message("Removing these genes prior to running correlations ..")
    genesToKeep <- Matrix::rowSums(RNAmat)!=0
    RNAmat <- RNAmat[genesToKeep,]
  }

  cat("Number of peaks in ATAC data:",nrow(ATACmat),"\n")
  cat("Number of genes in RNA data:",nrow(RNAmat),"\n")


  if (!genome %in% c("hg19", "hg38", "mm10"))
    stop("You must specify one of hg19, hg38 or mm10 as a genome build for currently supported TSS annotations..\n")
  switch(genome, hg19 = {
    TSSg <- FigR::hg19TSSRanges
  }, hg38 = {
    TSSg <- FigR::hg38TSSRanges
  }, mm10 = {
    TSSg <- FigR::mm10TSSRanges
  })

  # Keep genes that have annotation and are in RNA matrix
  names(TSSg) <- as.character(TSSg$gene_name)

  if(!is.null(geneList)){
    if(length(geneList)==1)
      stop("Please specify more than 1 valid gene symbol")

    if(any(!geneList %in% names(TSSg))){
      cat("One or more of the gene names supplied is not present in the TSS annotation specified: \n")
      cat(geneList[!geneList %in% names(TSSg)], sep = ", ")
      cat("\n")
      stop()
    }

    TSSg <- TSSg[geneList]
  }

  # Checking in case some genes in RNA don't overlap our TSS annotations
  genesToKeep <- intersect(names(TSSg),rownames(RNAmat))
  cat("\nNum genes overlapping TSS annotation and RNA matrix being considered: ",length(genesToKeep),"\n")

  # Match gene order in RNA matrix and TSS ranges
  RNAmat <- RNAmat[genesToKeep,]
  TSSg <- TSSg[genesToKeep]

  # Pad TSS by this much *either side*
  TSSflank <- GenomicRanges::flank(TSSg,
                                   width = windowPadSize,
                                   both = TRUE)

  # Get peak summit
  cat("\nTaking peak summits from peak windows ..\n")
  peakSummits <- resize(peakRanges,width = 1,fix = "center")

  # Find overlap of all peaks to all genes given window
  # Subject is Peaks, query is Gene
  cat("Finding overlapping peak-gene pairs ..\n")
  genePeakOv <- findOverlaps(query = TSSflank,subject = peakSummits)
  numPairs <- length(genePeakOv)

  cat("Found ",numPairs,"total gene-peak pairs for given TSS window ..\n")

  cat("Number of peak summits that overlap any gene TSS window: ",length(unique(subjectHits(genePeakOv))),"\n")
  cat("Number of gene TSS windows that overlap any peak summit: ",length(unique(queryHits(genePeakOv))),"\n\n")

  # For each gene, determine observed correlation of each overlapping peak to its associated gene (gene expression)

  # For each of those genes, also determine correlation based on background peaks (run in parallel) and save
  # Do this in parallel, and get data frame of gene-peak-pearson values
  # Fetch background peaks for each peak tested (i.e. that has overlap in window with gene)
  set.seed(123)
  cat("Determining background peaks ..\n")

  if(is.null(rowData(ATAC.se)$bias)){
    if(genome %in% "hg19")
      myGenome <- BSgenome.Hsapiens.UCSC.hg19::BSgenome.Hsapiens.UCSC.hg19
    if(genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if(genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38

    ATAC.se <- chromVAR::addGCBias(ATAC.se,genome=myGenome) }

  cat("Using ",n_bg," iterations ..\n\n")

  set.seed(123)
  bg <- chromVAR::getBackgroundPeaks(ATAC.se,niterations=n_bg)

  cat("Computing gene-peak correlations ..\n")

  pairsPerChunk <- 500

  # This defines the outer (larger chunks)
  largeChunkSize <- 5000

  startingPoint <- 1 # If for any reason the workers fail, resume from where it failed by specifying the starting point here
  chunkStarts <- seq(startingPoint, numPairs, largeChunkSize)
  chunkEnds <- chunkStarts + largeChunkSize -1
  chunkEnds[length(chunkEnds)] <- numPairs

  library(doParallel)

  dorcList <- list()
  for(i in 1:length(chunkStarts)){
    cat("Running pairs: ",chunkStarts[i], "to",chunkEnds[i],"\n")
    # This fill further chunk this up and run in parallel, saving the merged output ObsCor
    ObsCor <- FigR::PeakGeneCor(ATAC = ATACmat,
                                RNA = RNAmat,
                                OV = genePeakOv[chunkStarts[i]:chunkEnds[i]],
                                chunkSize = pairsPerChunk,
                                ncores = nCores,
                                bg = bg)
    gc()

    dorcList[[i]] <- ObsCor
  }

  cat("\nMerging results ..\n")
  dorcTab <- bind_rows(dorcList)

  cat("Performing Z-test for correlation significance ..\n")
  permCols <- 4:(ncol(bg)+3)


  if (keepPosCorOnly){
    # Filter to positive correlations
    cat("Only considering positive correlations ..\n")
    dorcTab <- dorcTab %>% dplyr::filter(rObs > 0)
  }

  if(!keepMultiMappingPeaks){
  # Remove multi-mapping peaks (force 1-1 mapping)
  cat("Keeping max correlation for multi-mapping peaks ..\n")
  dorcTab <- dorcTab %>% dplyr::group_by(Peak) %>% dplyr::filter(rObs==max(rObs))
  }

  # Swap gene number for gene symbol from TSS annotation lookup
  dorcTab$Gene <- as.character(TSSg$gene_name)[dorcTab$Gene]

  # Swap peak numbers to match reference input peak numbers
  # This only changes if some peaks had zero accessibility and were filtered out internally
  # Use rownames from reference matching
  dorcTab$Peak <- as.numeric(splitAndFetch(rownames(ATACmat)[dorcTab$Peak],"Peak",2))

  # # Z test pval
  dorcTab$rBgSD <- matrixStats::rowSds(as.matrix(dorcTab[,permCols]))
  dorcTab$rBgMean <- rowMeans(dorcTab[,permCols])
  dorcTab$pvalZ <- 1-stats::pnorm(q = dorcTab$rObs,mean = dorcTab$rBgMean,sd = dorcTab$rBgSD)


  cat("\nFinished!\n")

  if(!is.null(p.cut)){
    cat("Using significance cut-off of ",p.cut," to subset to resulting associations\n")
    dorcTab <- dorcTab[dorcTab$pvalZ <= p.cut,] # Subset to significant correlations only
  }

  # Add peak ranges to final data frame output
  dorcTab$PeakRanges <- paste(as.character(seqnames(peakRanges.OG[dorcTab$Peak])),paste(start(peakRanges.OG[dorcTab$Peak]),end(peakRanges.OG[dorcTab$Peak]),sep="-"),sep=":")

  return(as.data.frame(dorcTab[,c("Peak","PeakRanges","Gene","rObs","pvalZ")],stringsAsFactors=FALSE))
}

cisCorr <- runGenePeakcorr(
    ATACmat = atac_X, 
    ATAC.se = ATAC.se,
    RNAmat = rna_X,
    genome = organism,
    nCores = nCores,
    normalizeATACmat = FALSE,
    windowPadSize = round(ext / 2),
)

# Process
p2g <- cisCorr %>%
    filter(pvalZ <= thr_pval) %>%
    mutate(cre=stringr::str_replace(PeakRanges, ':', '-')) %>%
    rename(gene=Gene, score=rObs, pval=pvalZ) %>%
    select(cre, gene, score, pval)

# Identify dorcs and keep only their g-p pairs
dorcs <- p2g %>% summarize(counts=n(), .by=gene) %>% filter(counts >= ncres) %>% pull(gene)
p2g <- p2g %>% filter(gene %in% dorcs)

# Write
write.csv(x = p2g, file = path_out, row.names=FALSE)
