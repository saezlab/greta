library(tidyverse)
library(rhdf5)
library(Pando)
library(doParallel)
registerDoParallel(detectCores())


# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_p2g <- args[7]
path_tfb = args[8]
path_out = args[9]

# Read dfs
p2g <- read.csv(path_p2g)[, c('cre', 'gene')]
tfb <- read.csv(path_tfb)[, c('cre', 'tf')]

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_indices <- indata$mod$rna$X$indices
rna_indptr <- indata$mod$rna$X$indptr
rna_data <- as.numeric(indata$mod$rna$X$data)
barcodes <- indata$obs$`_index`
genes <- indata$mod$rna$var$`_index`
rna_X <- Matrix::sparseMatrix(i=rna_indices, p=rna_indptr, x=rna_data, index1 = FALSE)
colnames(rna_X) <- barcodes
row.names(rna_X) <- stringr::str_replace_all(genes, '-', '_')

### ATAC
atac_indices <- indata$mod$atac$X$indices
atac_indptr <- indata$mod$atac$X$indptr
atac_data <- as.numeric(indata$mod$atac$X$data)
peaks <- indata$mod$atac$var$`_index`
atac_X <- Matrix::sparseMatrix(i=atac_indices, p=atac_indptr, x=atac_data, index1 = FALSE)
colnames(atac_X) <- barcodes
row.names(atac_X) <- stringr::str_replace_all(peaks, '-', '_')
h5closeAll()

# Run per feature
p2g$gene = stringr::str_replace_all(p2g$gene, '-', '_')
features <- unique(p2g$gene)
model_fits <- Pando::map_par(features, function(g){
    # Subset scaffold
    g_scaff_grn <- inner_join(p2g[p2g$gene == g, ], tfb, by = 'cre') %>%
        filter(tf != gene) # Remove self regulation
    if (nrow(g_scaff_grn) == 0){
        return()
    }
    # Update strings
    g_scaff_grn$cre = stringr::str_replace_all(g_scaff_grn$cre, '-', '_')
    g_scaff_grn$tf = stringr::str_replace_all(g_scaff_grn$tf, '-', '_')
    
    # Extract data for given gene
    g_gex <- rna_X[c(g, unique(g_scaff_grn$tf)), , drop=FALSE]
    g_acc <- atac_X[unique(g_scaff_grn$cre), , drop=FALSE]
    model_mat <- as.data.frame(t(as.matrix(rbind(g_gex, g_acc))))

    # Generate formula
    formula <- g_scaff_grn %>%
    summarize(cres =
        paste(paste(cre, ':', paste('`', tf, '`', sep=''),sep=''), collapse = ' + '), .by=c(gene)
    ) %>%
    pull(cres)
    formula <- as.formula(paste(g, '~', formula))

    # Run model
    result <- fit_model(
        formula,
        data = model_mat,
        method = 'glm',
    )
    result$gof$nvariables <- nrow(g_scaff_grn)
    return(result)
}, parallel=TRUE)

# Format results
gof <- map_dfr(model_fits, function(x) x$gof, .id='target')
gof$target <- features[as.numeric(gof$target)]
coefs <- map_dfr(model_fits, function(x) x$coefs, .id='target')
coefs <- format_coefs(coefs, term=':', adjust_method='fdr')
coefs$target <- features[as.numeric(coefs$target)]

# Build network object from results
params <- list()
params[['method']] <- 'glm'
network_obj <- new(
    Class = 'Network',
    features = features,
    coefs = coefs,
    fit = gof,
    params = params
)

# Extract GRN
grn <- find_modules(
    network_obj,
)@modules@meta %>%
select(tf, target, estimate, padj) %>%
rename(source=tf, score=estimate, pval=padj) %>%
mutate(
    source=stringr::str_replace_all(source, '_', '-'),
    target=stringr::str_replace_all(target, '_', '-'),
)

# Write
write.csv(x = grn, file = path_out, row.names=FALSE)
