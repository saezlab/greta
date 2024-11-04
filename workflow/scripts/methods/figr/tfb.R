library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(doParallel)
library(FigR)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
genome <- args[7]
path_p2g <- args[8]
cellK <- as.numeric(args[9])
dorcK <- as.numeric(args[10])
nCores <- as.numeric(args[11])
path_out <- args[12]
n_bg <- 50

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`

# ATAC (read only raw counts)
ATAC.se <- Matrix::sparseMatrix(
    i=indata$mod$atac$layers$counts$indices,
    p=indata$mod$atac$layers$counts$indptr,
    x=as.numeric(indata$mod$atac$layers$counts$data),
    index1 = FALSE
)
colnames(ATAC.se) <- indata$obs$`_index`
rownames(ATAC.se) <- indata$mod$atac$var$`_index`

# Dim reduction
lsi <- indata$obsm$X_spectral
if (!is.null(lsi)){
    colnames(lsi) <- colnames(rna_X)
    rownames(lsi) <- paste0("SPL_", 1:nrow(lsi))
}
h5closeAll()

# Transform atac to sme Object
peaks <- strsplit(rownames(ATAC.se), "-")
peak_ranges <- GenomicRanges::GRanges(
    seqnames = sapply(peaks, "[[", 1),
    ranges = IRanges::IRanges(start = as.numeric(sapply(peaks, "[[", 2)), end = as.numeric(sapply(peaks, "[[", 3)))
)
ATAC.se <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=ATAC.se), rowRanges = peak_ranges)

# Read p2g
p2g <- read.csv(path_p2g)
if (nrow(p2g) == 0){
    tfb <- data.frame(cre=character(), tf=character(), score=numeric())
    write.csv(x = tfb, file = path_out, row.names=FALSE)
    quit(save="no")
}
p2g <- p2g %>%
    mutate(Peak=match(cre,  rownames(ATAC.se))) %>%
    rename(PeakRanges=cre, Gene=gene) %>%
    mutate(PeakRanges=sub("-", ":", PeakRanges, fixed = TRUE)) %>%
    select(Peak, PeakRanges, Gene)

# Compute sum of peaks per gene (DORCs)
dorcMat <- FigR::getDORCScores(
    ATAC.se = ATAC.se,
    dorcTab = p2g,
    normalizeATACmat=TRUE,
    nCores = nCores
)

# Smooth data
set.seed(123)
if (!is.null(lsi)){
    cellkNN <- FNN::get.knn(t(lsi), k = cellK)$nn.index
    rownames(cellkNN) <- as.character(colnames(rna_X))
    dorcMat <- FigR::smoothScoresNN(NNmat = cellkNN, mat = dorcMat, nCores = nCores)
}

get_TFenrich <- function(
    ATAC.se,
    dorcMat,
    rna_X,
    dorcK,
    genome,
    n_bg,
    p2g
    ){
    # Process
    dorcGenes <- rownames(dorcMat)
    DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
    rownames(DORC.knn) <- rownames(dorcMat)
    
    if (is.null(SummarizedExperiment::rowData(ATAC.se)$bias)) {
        if (genome %in% "mm10")
          myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
        if (genome %in% "hg38")
          myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
        ATAC.se <- chromVAR::addGCBias(ATAC.se, genome = myGenome)
    }
    
    # Set data subfolder path
    packagePath <- find.package("FigR", lib.loc=NULL, quiet = TRUE)
    
    # Read pwm matrices
    if(grepl("hg", genome)){
        pwm <- readRDS(paste0(packagePath,"/data/cisBP_human_pfms_2021.rds"))
    } else {
        pwm <- readRDS(paste0(packagePath,"/data/cisBP_mouse_pfms_2021.rds"))
    }
    if(all(grepl("_",names(pwm),fixed = TRUE)))
        names(pwm) <- FigR::extractTFNames(names(pwm))
    
    # Modify gene names
    #myGeneNames <- gsub(x = rownames(rna_X),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
    #rownames(rna_X) <- myGeneNames
    myGeneNames <- rownames(rna_X)
    
    # Only non-zero expression TFs (also found in rna_X)
    motifsToKeep <- intersect(names(pwm), myGeneNames)
    
    # Find motifs in peaks
    cat("Getting peak x motif matches ..\n")
    motif_ix <- motifmatchr::matchMotifs(subject = ATAC.se, pwms = pwm[motifsToKeep], genome=genome)
    motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix))!=0]
    
    # Select background peaks
    set.seed(123)
    bg <- chromVAR::getBackgroundPeaks(ATAC.se, niterations = n_bg)
    
    # Find enriched motifs
    #library(doParallel)
    opts <- list()
    pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
    progress <- function(n) setTxtProgressBar(pb, n)
    opts <- list(progress = progress)
    time_elapsed <- Sys.time()
    cl <- parallel::makeCluster(nCores)
    parallel::clusterEvalQ(cl, .libPaths())
    doSNOW::registerDoSNOW(cl)
    mZtest.list <- foreach(
        g=dorcGenes,
        .options.snow = opts,
        .packages = c("FigR", "dplyr","Matrix","Rmpfr")
    ) %dopar% {
        # Take peaks associated with gene and its k neighbors
        # Pool and use union for motif enrichment
        DORCNNpeaks <- unique(p2g$Peak[p2g$Gene %in% c(g,rownames(dorcMat)[DORC.knn[g,]])])
        if(length(DORCNNpeaks) <= 1){
            return()
        }
    
        # Z test
        mZ <- FigR::motifPeakZtest(
            peakSet = DORCNNpeaks,
            bgPeaks = bg,
            tfMat = assay(motif_ix)
        )
        # Process
        mZ <- mZ[,c("gene","z_test")]
        colnames(mZ)[1] <- "Motif"
        colnames(mZ)[2] <- "Enrichment.Z"
        mZ$Enrichment.P <- 2*pnorm(abs(mZ$Enrichment.Z),lower.tail = FALSE) # One-tailed
        mZ <- cbind("DORC"=g,mZ)
        mZ <- mZ %>% mutate(Enrichment.P=dplyr::if_else(Enrichment.Z < 0, 1, Enrichment.P)) # Negative enrichment penalize p-value
        return(mZ)
    }
    TFenrich.d <- do.call('rbind',mZtest.list)
    dim(TFenrich.d)
    rownames(TFenrich.d) <- NULL
    return(TFenrich.d)
}

TFenrich.d <- get_TFenrich(
    ATAC.se,
    dorcMat,
    rna_X,
    dorcK,
    genome,
    n_bg,
    p2g
)

# Process tfb
tfb <- left_join(
    rename(TFenrich.d, Gene=DORC),
    p2g %>% select(Gene, PeakRanges),
    relationship = "many-to-many",
    by='Gene') %>%
    summarize(Enrichment.P=mean(Enrichment.P), .by=c(Motif, PeakRanges)) %>% 
    mutate(score=-log10(Enrichment.P)) %>%
    rename(tf=Motif, cre=PeakRanges) %>%
    mutate(cre=sub(":", "-", cre, fixed = TRUE)) %>%
    select(cre, tf, score) %>%
    arrange(cre, desc(score))

if (any(is.infinite(tfb$score))) {
  max_finite <- max(tfb$score[is.finite(tfb$score)])
  tfb$score[is.infinite(tfb$score)] <- max_finite
}

# Write
write.csv(x = tfb, file = path_out, row.names=FALSE)
