library(tidyverse)
library(rhdf5)
library(Pando)
library(doParallel)

nCores <- 64
cat("N cores: ", nCores, '\n')
registerDoParallel(nCores)

# Parse args
args <- commandArgs(trailingOnly = F)
path_data <- args[6]
path_p2g <- args[7]
path_tfb = args[8]
thr_cor = 0.05
thr_rsq = 0.05
path_out = args[9]

# Read dfs
p2g <- read.csv(path_p2g)[, c('cre', 'gene')]
tfb <- read.csv(path_tfb)[, c('cre', 'tf')]
if ((nrow(p2g) == 0) | (nrow(tfb) == 0)){
    mdl <- data.frame(source=character(), target=character(), score=numeric(), pval=numeric())
    write.csv(x = mdl, file = path_out, row.names=FALSE)
    quit(save="no")
}

# Read data
print('Open object')
indata <- H5Fopen(path_data, flags='H5F_ACC_RDONLY')
# RNA
rna_X <- indata$mod$rna$X
colnames(rna_X) <- indata$obs$`_index`
rownames(rna_X) <- stringr::str_replace_all(indata$mod$rna$var$`_index`, '-', '_')

### ATAC
atac_X <- indata$mod$atac$X
colnames(atac_X) <- indata$obs$`_index`
rownames(atac_X) <- stringr::str_replace_all(indata$mod$atac$var$`_index`, '-', '_')
h5closeAll()

# Run per feature
p2g$gene = stringr::str_replace_all(p2g$gene, '-', '_')
features <- unique(p2g$gene)
cat("Number of target genes to fit: ", length(features), '\n')
model_fits <- Pando::map_par(features, function(g){
    # Subset scaffold
    print(g)
    g_scaff_grn <- inner_join(p2g[p2g$gene == g, ], tfb, by = 'cre') %>%
        filter(tf != gene) # Remove self regulation

    # Update strings
    g_scaff_grn$cre = stringr::str_replace_all(g_scaff_grn$cre, '-', '_')
    g_scaff_grn$tf = stringr::str_replace_all(g_scaff_grn$tf, '-', '_')
    
    # Extract data for given gene
    g_gex <- t(rna_X[c(g, unique(g_scaff_grn$tf)), , drop=FALSE])
    g_acc <- t(atac_X[unique(g_scaff_grn$cre), , drop=FALSE])
    
    # Filter by min corr
    pks <- cor(g_acc, g_gex[, g, drop=FALSE]) %>%
        as.data.frame() %>%
        arrange(!!as.symbol(g)) %>%
        filter(abs(!!as.symbol(g)) > thr_cor) %>%
        rownames()
    gns <- cor(g_gex, g_gex[, g, drop=FALSE]) %>%
        as.data.frame() %>%
        arrange(!!as.symbol(g)) %>%
        filter(abs(!!as.symbol(g)) > thr_cor) %>%
        rownames()
    
    g_scaff_grn <- g_scaff_grn %>% filter(cre %in% pks, tf %in% gns)
    n_covs = nrow(g_scaff_grn)
    if ((n_covs == 0) | ((n_covs + 2) > ncol(rna_X))){
        return()
    }
    model_mat <- as.data.frame(cbind(g_gex[, gns, drop=FALSE], g_acc[, pks, drop=FALSE]))

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

# Filter features by failed models
msk <- unlist(map(model_fits, function(x){!is.null(x)}))
features <- features[msk]

# Format results
gof <- map_dfr(model_fits, function(x) x$gof, .id='target')
if (nrow(gof) == 0){
    mdl <- data.frame(source=character(), target=character(), score=numeric(), pval=numeric())
    write.csv(x = mdl, file = path_out, row.names=FALSE)
    quit(save="no")
}
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
    rsq_thresh = thr_rsq,
    p_thresh = 0.1,
    nvar_thresh = 2,
    min_genes_per_module = 3
)@modules@meta %>%
select(tf, target, estimate, padj) %>%
rename(source=tf, score=estimate, pval=padj) %>%
mutate(
    source=stringr::str_replace_all(source, '_', '-'),
    target=stringr::str_replace_all(target, '_', '-'),
)

# Write
write.csv(x = grn, file = path_out, row.names=FALSE)
