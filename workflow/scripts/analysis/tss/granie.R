# Initiate args
args <- commandArgs(trailingOnly = F)
path_out <- args[6]


# They use the newest version of the homo sapiens
library(AnnotationHub)
ah <- AnnotationHub()

# Get the newest version of annotation
results = AnnotationHub::query(ah, c("EnsDb", "Homo sapiens"))
annotationDatasets <- as.data.frame(mcols(results))
newestAnno.title = tail(annotationDatasets$title, 1)
newestAnno.ID = tail(rownames(annotationDatasets), 1)
ensdb.newest <- ah[[newestAnno.ID]]


# Select columns of interest
genes.df = as.data.frame(suppressWarnings(ensembldb::genes(ensdb.newest))) %>%
  tibble::as_tibble() %>%
  dplyr::mutate(gene.chr = paste0("chr", .data$seqnames)) %>%
  dplyr::select(-"seqnames") %>%
  dplyr::rename(gene.ENSEMBL = "gene_id", gene.start = "start", gene.end = "end",
                gene.strand = "strand", gene.name = "gene_name", gene.type = "gene_biotype") %>%
  dplyr::select("gene.chr", "gene.start", "gene.end", "gene.strand", "gene.ENSEMBL", "gene.type", "gene.name") %>%
  tidyr::replace_na(list(gene.type = "unknown")) %>%
  dplyr::mutate(gene.strand = factor(.data$gene.strand, levels = c(1,-1,0), labels = c("+", "-", "*"))) %>%
  dplyr::mutate_if(is.character, as.factor) %>%
  dplyr::mutate(gene.type = dplyr::recode_factor(.data$gene.type, lncRNA = "lincRNA")) 

# Filter only protein coding
genes.df <- genes.df %>% dplyr::filter(gene.type == "protein_coding")


# Keep interest columns
genes.df <- genes.df %>%
  dplyr::rename(Chromosome = gene.chr, 
         Start = gene.start, 
         End = gene.end, 
         Name = gene.name) %>%
  dplyr::select(Chromosome, Start, End, Name) 


write.csv(x = genes.df, file = path_out)
