library(biomaRt)

# Parse args
orgms <- commandArgs(trailingOnly = TRUE)

get_gene_table <- function(dataset){
    # Connect to the Ensembl database
    ensembl <- useEnsembl(
        biomart = 'genes',
        dataset = dataset,
        version = 111
    )
    # Specify the attributes to retrieve
    attributes <- c("uniprotswissprot", "external_gene_name")
    # Retrieve the data
    gene_data <- getBM(
        attributes = attributes,
        mart = ensembl,
        useCache=FALSE,
        verbose=FALSE
    )
    colnames(gene_data) <- c('uniprot_id', 'symbol')
    return(gene_data)
}

org_table <- list(
    'hg38'='hsapiens_gene_ensembl',
    'mm10'='mmusculus_gene_ensembl'
)

for (path_org in orgms) {
    org <- sub('^dbs/([^/]+)/.*$', '\\1', path_org)
    org <- org_table[org]
    gid <- get_gene_table(org)
    gid <- gid[(gid$symbol != "") & (gid$uniprot_id != ""), ]  # Exclude rows with empty gene symbols
    write.csv(x = gid, file = path_org, row.names=FALSE, quote=FALSE)
}
