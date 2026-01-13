library(Seurat)
library(Signac) 
library(Matrix)
library(TFBSTools)
library(motifmatchr)
library(rhdf5)
library(EnsDb.Hsapiens.v86)
library(CREMA)
library(parallel)
library(BSgenome.Hsapiens.UCSC.hg38)


# Parse args
args <- commandArgs(trailingOnly = F)
print(args)
path_data <- args[6]
path_inserts <- args[7]
path_motifdb <- args[8]
ext <- as.numeric(args[9])
site_extension <- as.numeric(args[10])
thr_fdr <- as.numeric(args[11])
path_out <- args[12]
threads <- as.numeric(args[13])
path_hvg <- args[14]

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
rna <- Matrix::sparseMatrix(
    i=indata$mod$rna$layers$counts$indices,
    p=indata$mod$rna$layers$counts$indptr,
    x=as.numeric(indata$mod$rna$layers$counts$data),
    index1 = FALSE
)
colnames(rna) <- indata$obs$`_index`
rownames(rna) <- indata$mod$rna$var$`_index`
h5closeAll()

# Filter hvg
hvg <- readLines(gzfile(path_hvg))
rna <- rna[hvg, ]

# Format
rna <- CreateSeuratObject(
  counts = rna,
  assay = "RNA"
)
DefaultAssay(rna) <- "RNA"
rna <- SCTransform(rna)

# Read inserts
frag_inserts_object <- CreateFragmentObject(path=path_inserts, cells=colnames(rna))

# Format gene coords
print('Format gene coords')
genebody_coords <- keepStandardChromosomes(
    ensembldb::genes(EnsDb.Hsapiens.v86),
    species="Homo_sapiens",
    pruning.mode = 'coarse'
)
seqlevels(genebody_coords) <- paste0("chr", seqlevels(genebody_coords))
genebody_coords <- genebody_coords[sapply(genebody_coords$gene_name, nchar) > 0]
temp_ind <- table(genebody_coords$gene_name)
temp_ind <- names(temp_ind)[temp_ind > 1]
genebody_coords_list <- split(genebody_coords, f = genebody_coords$gene_name)
temp_func <- function(x){
    if ("protein_coding" %in% x$gene_biotype){ return(x[x$gene_biotype == "protein_coding"])
    }else{ return(x) }
}
genebody_coords_list <- c(
    lapply(genebody_coords_list[temp_ind], temp_func),
    as.list(genebody_coords_list[setdiff(names(genebody_coords_list), temp_ind)])
)
genebody_coords <- unlist(as(genebody_coords_list, "GRangesList"))
genebody_coords <- sort(genebody_coords)
genebody_coords <- trim(genebody_coords)

# Load TF motifs
motifs <- readRDS(path_motifdb)
motif_tfs <- unique(sapply(strsplit(names(motifs), split = "_"), function(x){x[1]}))

# Filter
print('Filter assay')
exp_mtx <- GetAssayData(rna, assay = "RNA", slot = "counts")
exp_mtx <- filter_exp_mtx_genes(exp_mtx, gene_names_remove_pattern = "^MT-", proportion_cells_detected = 1e-3)
TFs_select <- intersect(motif_tfs, row.names(exp_mtx))
genes_select <- row.names(exp_mtx)
chromosomes <- paste0("chr", c(seq(1,22), "X", "Y"))
genes_select <- genes_select[genes_select %in% genebody_coords$symbol[as.character(seqnames(genebody_coords)) %in% chromosomes]]
print(paste("num of TFs selected:", length(TFs_select)))
print(paste("num of genes selected:", length(genes_select)))
exp_mtx <- GetAssayData(rna, assay = "SCT", slot = "counts")
# Filter just by TFs else it cannot scale
TFs_select <- intersect(TFs_select, row.names(exp_mtx))
#genes_select <- intersect(genes_select, row.names(exp_mtx))
#exp_mtx <- as.matrix(exp_mtx[union(TFs_select, genes_select), ])
genes_select <- TFs_select
exp_mtx <- as.matrix(exp_mtx[TFs_select, ])
genebody_coords <- genebody_coords[genebody_coords$gene_name %in% genes_select, ]

# Open windows
print('Open windows')
genebody_coords <- genebody_coords[genebody_coords$gene_name %in% genes_select]
tss <- ifelse(
    test = (GenomicRanges::strand(genebody_coords) == "+" | GenomicRanges::strand(genebody_coords) == "*"),
    yes = GenomicRanges::start(genebody_coords), no = GenomicRanges::end(genebody_coords)
)
msk <- !is.na(tss)
genebody_coords <- genebody_coords[msk, ]

crema_regions <- select_proximal_regions(
    genes = genes_select,
    gene_body_gr = genebody_coords, 
    window_up = ext,
    window_down = ext,
)
crema_regions_str <- lapply(crema_regions, Signac::GRangesToString)

# Run for all genes
test_genes <- rownames(exp_mtx)
print(paste("Running on TFs selected:", length(test_genes)))
if (length(test_genes) == 0) {
    writeLines("source,cre,target,score", path_out)
    quit(save = "no", status = 0)
}
print('Running across genes')
grn <- mclapply(test_genes, function(test_gene) {
    res_gene <- ATAC_weighted_tf_model_highres(
        test_gene,
        TFs = TFs_select,
        regions_str = crema_regions_str[[test_gene]],
        exp_mtx = exp_mtx,
        motifs = motifs,
        fragment_object = frag_inserts_object,
        cells = colnames(rna),
        genome = BSgenome.Hsapiens.UCSC.hg38,
        return_val = "df",
        regression_method = "ols",
        site_extension = site_extension,
    )
    if (is.data.frame(res_gene)) {
        res_gene$target <- test_gene
        res_gene <- res_gene[p.adjust(res_gene$`Pr(>|t|)`, method = "fdr") < thr_fdr, ]
    }
    res_gene
}, mc.cores = threads)
print('Merging')
grn <- do.call(rbind, grn[!is.na(grn)])
grn <- grn[, c('Estimate', 'Pr(>|t|)', 'target')]
colnames(grn) <- c('score', 'pval', 'target')
grn$source <- sapply(strsplit(rownames(grn), "_"), `[`, 1)
grn$cre <- sapply(strsplit(rownames(grn), "_"), `[`, 2)
rownames(grn) <- NULL
grn <- grn[, c('source', 'cre', 'target', 'score', 'pval')]
grn <- grn[order(grn$source, grn$target, grn$cre), ]

# Write
write.csv(x = grn, file = path_out, row.names=FALSE)
