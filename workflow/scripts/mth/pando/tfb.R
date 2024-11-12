library(tidyverse)
library(rhdf5)
library(Pando)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
organism <- args[7]
path_p2g <- args[8]
path_out <- args[9]

# Read genome
if (organism == 'hg38'){
    library(BSgenome.Hsapiens.UCSC.hg38)
    genome <- BSgenome.Hsapiens.UCSC.hg38
} else {
    library(BSgenome.Mmusculus.UCSC.mm10)
    genome <- BSgenome.Mmusculus.UCSC.mm10
}

# Read p2g
p2g <- read.csv(path_p2g)
if (nrow(p2g) == 0){
    tfb <- data.frame(cre=character(), tf=character(), score=numeric())
    write.csv(x = tfb, file = path_out, row.names=FALSE)
    quit(save="no")
}

# Read genes
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
genes <- indata$mod$rna$var$`_index`
h5closeAll()

# Transform motif2tf to mat
data('motif2tf')
motif2tf <- motif2tf %>% select('motif'=1,'tf'=2) %>%
    distinct() %>% mutate(val=1) %>%
    tidyr::pivot_wider(names_from = 'tf', values_from=val, values_fill=0) %>%
    tibble::column_to_rownames('motif') %>%
    as.matrix() %>% Matrix::Matrix(sparse=TRUE)

# Subset motifs to tfs in data
data('motifs')
motif2tf <- motif2tf[, intersect(genes, colnames(motif2tf))]
motifs <- motifs[rownames(motif2tf)[Matrix::rowSums(motif2tf) != 0 ]]

# Transform peaks to Granger
peaks <- data.frame(seqnames=p2g$cre) %>% distinct()
peaks <- tidyr::separate(data = peaks, col = 'seqnames', into = c("seqnames", "start", "end"), sep = "-", remove=FALSE)
peaks <- GenomicRanges::makeGRangesFromDataFrame(peaks)

# Run motif enrichment using motifmatcher (MOODS)
peak_motifs <- Signac::CreateMotifMatrix(
    features = peaks,
    pwm = motifs,
    genome = genome,
    use.counts = FALSE,
    score=TRUE
)

# Extact list of motifs to tfs
sparse <- summary(motif2tf)
motif2tf_lst <- data.frame(
    motif = rownames(motif2tf)[sparse$i],
    tf = colnames(motif2tf)[sparse$j]
) %>%
group_by(motif) %>%
summarize(values = list(tf)) %>%
deframe()

# Convert from sparse mat to df
sparse <- summary(peak_motifs)
df <- data.frame(
    cre = rownames(peak_motifs)[sparse$i],
    tf = colnames(peak_motifs)[sparse$j],
    score = sparse$x
) %>%
mutate(tf=motif2tf_lst[tf]) %>%
unnest(tf) %>%
summarize(score = max(score), .by=c(cre, tf)) %>%
mutate(score = ifelse(score < 0, 0, score)) %>%  # Sometimes MOODs returns negative values
arrange(cre, desc(score))

# Write
write.csv(x = df, file = path_out, row.names=FALSE)
