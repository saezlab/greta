library(DIRECTNET)
library(Seurat)
library(Signac)
library(EnsDb.Hsapiens.v86)
library(patchwork)
library(dplyr)
library(ggplot2)
options(stringsAsFactors = FALSE)
library(rhdf5)


# Parse args
args <- commandArgs(trailingOnly = F)
print(args)
path_data <- args[6]
path_ann = args[7]
path_tss = args[8]
ext = as.numeric(args[9])
path_grn <- args[10]

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
rna_X <- Matrix::sparseMatrix(
    i=indata$mod$rna$layers$counts$indices,
    p=indata$mod$rna$layers$counts$indptr,
    x=as.numeric(indata$mod$rna$layers$counts$data),
    index1 = FALSE
)
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- indata$mod$rna$var$`_index`
atac_X <- Matrix::sparseMatrix(
    i=indata$mod$atac$layers$counts$indices,
    p=indata$mod$atac$layers$counts$indptr,
    x=as.numeric(indata$mod$atac$layers$counts$data),
    index1 = FALSE
)
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- indata$mod$atac$var$`_index`
obs = data.frame(row.names=indata$obs$`_index`)
obs$celltype = factor(indata$obs$celltype$categories[indata$obs$celltype$codes + 1], levels = indata$obs$celltype$categories)
h5closeAll()

# Create object
muo_data <- Seurat::CreateSeuratObject(counts=rna_X, meta.data=obs)
muo_data[['ATAC']] <- Signac::CreateChromatinAssay(counts = atac_X)
annot <- read.csv(path_ann)
library(BSgenome.Hsapiens.UCSC.hg38)
genome <- BSgenome.Hsapiens.UCSC.hg38
annot <- GenomicRanges::makeGRangesFromDataFrame(annot, keep.extra.columns=TRUE)
GenomeInfoDb::seqlevelsStyle(annot) <- 'UCSC'
Signac::Annotation(muo_data[['ATAC']]) <- annot
genome.info <- read.csv(path_tss, sep='\t', header=FALSE)
names(genome.info) <- c("Chrom","Starts","Ends","genes")
genome.info <- genome.info[!grepl("_", genome.info$Chrom), ]
genes <- lapply(genome.info$genes, function(x) strsplit(x,"[|]")[[1]][1])
genes <- lapply(genes, function(x) strsplit(x,"[.]")[[1]][1])
genes <- unlist(genes)
genome.info$genes <- genes
unik <- !duplicated(genes)
genome.info <- genome.info[unik,]

# Process
DefaultAssay(muo_data) <- "RNA"
muo_data <- SCTransform(muo_data, verbose = FALSE) %>% RunPCA() %>% RunUMAP(dims = 1:50, reduction.name = 'umap.rna', reduction.key = 'rnaUMAP_')
DefaultAssay(muo_data) <- "ATAC"
muo_data <- RunTFIDF(muo_data)
muo_data <- FindTopFeatures(muo_data, min.cutoff = 'q0')
muo_data <- RunSVD(muo_data)
muo_data <- RunUMAP(muo_data, reduction = 'lsi', dims = 2:50, reduction.name = "umap.atac", reduction.key = "atacUMAP_")
muo_data <- FindMultiModalNeighbors(muo_data, reduction.list = list("pca", "lsi"), dims.list = list(1:50, 2:50))
muo_data <- RunUMAP(muo_data, nn.name = "weighted.nn", reduction.name = "wnn.umap", reduction.key = "wnnUMAP_")
Idents(muo_data) <- muo_data@meta.data$celltype

# Find marker genes and peaks
library(presto)
markers_all <- presto:::wilcoxauc.Seurat(X = muo_data, group_by = 'celltype', assay = 'data', seurat_assay = 'SCT')
markers_all <- markers_all[which(markers_all$auc > 0.5), , drop = FALSE]
markers <- data.frame(gene = markers_all$feature, group = markers_all$group)
groups <- unique(markers$group)
marker_list <- list()
for (i in 1:length(groups)) {
  marker1<- markers_all[markers$group == groups[i],]
  marker_list[[i]] <- as.character(marker1$feature[marker1$auc > 0.5])
}
markers_groups <- unique(unlist(marker_list))
markers_groups <- lapply(markers_groups, function(x) strsplit(x,"[.]")[[1]][1])
markers_groups <- unique(unlist(markers_groups))
groups <- levels(Idents(muo_data))
da_peaks_list <- list()
for (i in 1:length(groups)) {
  print(i)
  da_peaks <- FindMarkers(
    object = muo_data,
    min.pct = 0.2,
    logfc.threshold = 0.6,
    ident.1 = groups[i],
    group.by = "celltype",
    test.use = 'LR',
    only.pos = TRUE
  )
  da_peaks_list[[i]] <- da_peaks
}

# Run GRN
k_neigh <- min(table(muo_data@meta.data$celltype)) - 1
muo_data <- Run_DIRECT_NET(muo_data, k_neigh = k_neigh, atacbinary=TRUE, max_overlap=0.5, size_factor_normalize=TRUE,
                           genome.info=genome.info, focus_markers=markers_groups, window = ext)
direct.net_result <- Misc(muo_data, slot = 'direct.net')
direct.net_result <- as.data.frame(do.call(cbind,direct.net_result))
CREs_Gene <- generate_CRE_Gene_links(direct.net_result, markers = markers)
Focused_CREs <- generate_CRE(L_G_record = CREs_Gene$distal, P_L_G_record = CREs_Gene$promoter, da_peaks_list=da_peaks_list)
L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$distal, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = markers)
if (length(L_TF_record) == 0) {
    writeLines("source,cre,target,score", path_grn)
    quit(save = "no", status = 0)
}
P_L_TF_record <- generate_peak_TF_links(peaks_bed_list = Focused_CREs$promoter, species="Homo sapiens", genome = BSgenome.Hsapiens.UCSC.hg38, markers = markers)
msk <- !sapply(P_L_TF_record, is.null)
network_links <- generate_links_for_Cytoscape(L_G_record = Focused_CREs$L_G_record[msk], L_TF_record[msk], P_L_G_record = Focused_CREs$P_L_G_record[msk], P_L_TF_record[msk], groups=groups[msk])

# Format GRN
grn <- network_links[, c('TF', 'Gene')]
colnames(grn) <- c('source', 'gene')
p2g <- c(Focused_CREs$L_G_record[msk],
                 Focused_CREs$P_L_G_record[msk])
p2g <- do.call(rbind, p2g)
p2g <- p2g[!duplicated(p2g[c("loci", "gene")]), ]
rownames(p2g) <- NULL
colnames(p2g) <- c('cre', 'gene')
grn <- merge(grn, p2g, by = "gene", all = FALSE)[, c('source', 'cre', 'gene')]
colnames(grn) <- c('source', 'cre', 'target')
grn$score <- 1
grn$cre <- gsub("_", "-", grn$cre)

# Write
write.csv(x = grn, file = path_grn, row.names=FALSE)
