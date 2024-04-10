library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_p2g <- args[7]
path_tfb <- args[8]
score_thr <- as.numeric(args[9])
path_out <- args[10]
nCores <- 1

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`

# ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
h5closeAll()

# Transform atac to sme Object
peaks <- strsplit(rownames(atac_X), "-")
peak_ranges <- GenomicRanges::GRanges(
    seqnames = sapply(peaks, "[[", 1),
    ranges = IRanges::IRanges(start = as.numeric(sapply(peaks, "[[", 2)), end = as.numeric(sapply(peaks, "[[", 3)))
)
atac_X <- SummarizedExperiment::SummarizedExperiment(assays = list(counts=atac_X), rowRanges = peak_ranges)

# Read p2g
p2g <- read.csv(path_p2g) %>%
    mutate(Peak=match(cre,  rownames(atac_X))) %>%
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

# Compute sum of peaks per gene (DORCs)
dorcMat <- FigR::getDORCScores(
    ATAC.se = atac_X,
    dorcTab = p2g,
    normalizeATACmat=FALSE,
    nCores = nCores
)

# Process
dorcGenes <- rownames(dorcMat)

# Find enriched motifs
library(doParallel)
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
    rename(source=tf)

# Write
write.csv(x = mdl, file = path_out, row.names=FALSE)
