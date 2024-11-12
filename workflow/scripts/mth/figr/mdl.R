library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)
library(doParallel)
library(FigR)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_p2g <- args[7]
path_tfb <- args[8]
cellK <- as.numeric(args[9])
score_thr <- as.numeric(args[10])
nCores <- as.numeric(args[11])
path_out <- args[12]

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- as.character(indata$obs$`_index`)
rownames(rna_X) <- as.character(indata$mod$rna$var$`_index`)

# ATAC (read only raw counts)
ATAC.se <- Matrix::sparseMatrix(
    i=indata$mod$atac$layers$counts$indices,
    p=indata$mod$atac$layers$counts$indptr,
    x=as.numeric(indata$mod$atac$layers$counts$data),
    index1 = FALSE
)
colnames(ATAC.se) <- as.character(indata$obs$`_index`)
rownames(ATAC.se) <- as.character(indata$mod$atac$var$`_index`)

# Dim reduction
lsi <- indata$obsm$X_spectral
if (!is.null(lsi)){
    colnames(lsi) <- as.character(colnames(rna_X))
    rownames(lsi) <- as.character(paste0("SPL_", 1:nrow(lsi)))
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
p2g <- read.csv(path_p2g) %>%
    mutate(Peak=match(cre,  rownames(ATAC.se))) %>%
    rename(PeakRanges=cre, Gene=gene) %>%
    mutate(PeakRanges=sub("-", ":", PeakRanges, fixed = TRUE)) %>%
    select(Peak, PeakRanges, Gene)
tfb <- read.csv(path_tfb) %>%
    mutate(cre=sub("-", ":", cre, fixed = TRUE)) %>%
    rename(PeakRanges=cre)
if ((nrow(p2g) == 0) | (nrow(tfb) == 0)){
    mdl <- data.frame(source=character(), target=character(), score=numeric(), pval=numeric())
    write.csv(x = mdl, file = path_out, row.names=FALSE)
    quit(save="no")
}

# Make negative scores to 0
tfb$score[tfb$score < 0] <- 0

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
    rna_X <- FigR::smoothScoresNN(NNmat = cellkNN, mat = rna_X, nCores = nCores)
}
# Process
dorcGenes <- as.character(rownames(dorcMat))

# Find enriched motifs
opts <- list()
pb <- txtProgressBar(min = 0, max = length(dorcGenes), style = 3)
progress <- function(n) setTxtProgressBar(pb, n)
opts <- list(progress = progress)
time_elapsed <- Sys.time()
cl <- parallel::makeCluster(nCores)
parallel::clusterEvalQ(cl, .libPaths())
doSNOW::registerDoSNOW(cl)

get_pval <- function(x, corr){
    # Number of observations
    n <- length(x)
    # Standard error of the Spearman correlation coefficient
    standard_error <- 1 / sqrt(n - 3)
    # Z-score
    z <- 0.5 * log((1 + corr) / (1 - corr))
    z_score <- z / standard_error
    # Calculate p-value
    p_value <- 2 * (1 - pnorm(abs(z_score)))
    return(p_value[1,])
}

mZtest.list <- foreach(
    g=dorcGenes,
    .options.snow = opts,
    .packages = c("FigR", "dplyr","Matrix","Rmpfr")
) %dopar% {
    mZ <- inner_join(p2g[p2g$Gene == g, ], tfb, by = 'PeakRanges') %>%
    summarize(tf_pval=mean(score), .by=c(tf)) %>%
    mutate(tf_pval = 10 ** -tf_pval)
    if(nrow(mZ) <= 1){
        return()
    }
    
    corr.r <- cor(dorcMat[g,], t(as.matrix(rna_X[mZ$tf,])), method = "spearman")
    stopifnot(all(colnames(corr.r) == mZ$tf))
    
    mZ$Corr <- corr.r[1,] # Correlation coefficient
    mZ$Corr.Z <- scale(mZ$Corr,center = TRUE,scale = TRUE)[,1] # Z-score among all TF correlations
    mZ$Corr.P <- 2*pnorm(abs(mZ$Corr.Z),lower.tail = FALSE) # One-tailed
    
    mZ <- mZ %>%
    mutate(score=sign(Corr)*as.numeric(-log10(1-(1-Rmpfr::mpfr(tf_pval,100))*(1-Rmpfr::mpfr(Corr.P,100))))) %>%
    rowwise() %>%
    mutate(target = g, pval=mean(c(tf_pval, Corr.P))) %>%
    ungroup() %>%
    select(tf, target, score, pval)
    return(mZ)
}
mdl <- do.call('rbind', mZtest.list) %>%
    filter(abs(score) > score_thr) %>%
    rename(source=tf) %>%
    arrange(source, target)

# Write
write.csv(x = mdl, file = path_out, row.names=FALSE)
