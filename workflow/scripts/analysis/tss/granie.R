library(AnnotationHub)


# Initiate args
args <- commandArgs(trailingOnly = F)
path_out <- args[6]


# Load db
ah <- AnnotationHub()

# Get the newest version of annotation
results = AnnotationHub::query(ah, c("EnsDb", "Homo sapiens"))
annotationDatasets <- as.data.frame(mcols(results))
newestAnno.title = tail(annotationDatasets$title, 1)
newestAnno.ID = tail(rownames(annotationDatasets), 1)
ensdb.newest <- ah[[newestAnno.ID]]

# Read
gr <- ensembldb::genes(ensdb.newest)

# Merge overlaps
merged <- unlist(reduce(split(gr, gr$gene_name)), use.names = TRUE)

# To df
chr_names <- paste0("chr", as.character(seqnames(merged)))
start_pos <- start(merged) - 1
end_pos <- end(merged) - 1
gene_names <- names(merged)
bed <- data.frame(Chromosome = chr_names, Start = start_pos, End = end_pos, Name = gene_names)

# Filter empty names
bed <- bed[bed$Name != '', ]

# Sort
bed <- bed[order(bed$Chromosome, bed$Start, bed$End), ]

# Write
write.table(x = bed, file = path_out, sep = '\t', row.names = FALSE, quote = FALSE, col.names = FALSE)
