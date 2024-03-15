library(biomaRt)

# Parse args
args <- commandArgs(trailingOnly = F)
path_hg <- args[6]
path_mm <- args[7]

get_gene_table <- function(dataset){
    # Connect to the Ensembl database
    ensembl <- useEnsembl(
        biomart = 'genes',
        dataset = 'hsapiens_gene_ensembl',
        version = 111
    )
    # Specify the attributes to retrieve
    attributes <- c("ensembl_gene_id", "external_gene_name")
    # Retrieve the data
    gene_data <- getBM(
        attributes = attributes,
        mart = ensembl,
        useCache=FALSE
    )
    colnames(gene_data) <- c('id', 'symbol')
    return(gene_data)
}

hg38 <- get_gene_table('hsapiens_gene_ensembl')
mm10 <- get_gene_table('mmusculus_gene_ensembl')

# Write
write.csv(x = hg38, file = path_hg, row.names=FALSE)
write.csv(x = mm10, file = path_mm, row.names=FALSE)
