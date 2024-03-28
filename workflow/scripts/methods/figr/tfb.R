library(rhdf5)
library(GenomicRanges)
library(SummarizedExperiment)
library(dplyr)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
genome <- args[7]
path_p2g <- args[8]
dorcK <- as.numeric(args[9])
path_out <- args[10]
nCores <- 1
n_bg <- 50

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

# Compute sum of peaks per gene (DORCs)
dorcMat <- FigR::getDORCScores(
    ATAC.se = atac_X,
    dorcTab = p2g,
    normalizeATACmat=FALSE,
    nCores = nCores
)

# Process
dorcGenes <- rownames(dorcMat)
DORC.knn <- FNN::get.knn(data = t(scale(Matrix::t(dorcMat))),k = dorcK)$nn.index # Scaled
rownames(DORC.knn) <- rownames(dorcMat)

if (is.null(SummarizedExperiment::rowData(atac_X)$bias)) {
    if (genome %in% "mm10")
      myGenome <- BSgenome.Mmusculus.UCSC.mm10::BSgenome.Mmusculus.UCSC.mm10
    if (genome %in% "hg38")
      myGenome <- BSgenome.Hsapiens.UCSC.hg38::BSgenome.Hsapiens.UCSC.hg38
    atac_X <- chromVAR::addGCBias(atac_X, genome = myGenome)
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
myGeneNames <- gsub(x = rownames(rna_X),pattern = "-",replacement = "") # NKX2-1 to NKX21 (e.g.)
rownames(rna_X) <- myGeneNames

# Only non-zero expression TFs (also found in rna_X)
motifsToKeep <- intersect(names(pwm),myGeneNames)

# Find motifs in peaks
cat("Getting peak x motif matches ..\n")
motif_ix <- motifmatchr::matchMotifs(subject = atac_X, pwms = pwm[motifsToKeep], genome=genome)
motif_ix <- motif_ix[, Matrix::colSums(assay(motif_ix))!=0]

# Select background peaks
set.seed(123)
bg <- chromVAR::getBackgroundPeaks(atac_X, niterations = n_bg)

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
    #mZ <- mZ %>% filter(Enrichment.P < 0.05)
    return(mZ)
}
TFenrich.d <- do.call('rbind',mZtest.list)
dim(TFenrich.d)
rownames(TFenrich.d) <- NULL

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

# Write
write.csv(x = tfb, file = path_out, row.names=FALSE)
