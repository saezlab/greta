library(argparse)
library(readr)
library(qs)
library(checkmate)

parser <- ArgumentParser(description= 'Postprocess GRaNIE results')
parser$add_argument('--input', '-i', help= 'Output file from GRaNIE (GRN object)', nargs = 1,  required= TRUE)
parser$add_argument('--output', '-o', help= 'Filenames for processed output files: 1. TF-gene, 2. TF-peak-gene', nargs= 2, required= TRUE)
parser$add_argument('--ranking_TF_gene', help= 'Ranking of TF-genes', nargs= 1, required= TRUE)
parser$add_argument('--ranking_TF_peak_gene', help= 'Ranking of TF-peak-genes', nargs= 1, required= TRUE)


xargs <- parser$parse_args()


checkmate::assertSubset(xargs$ranking_TF_gene, c("TF_gene.p_raw", "TF_gene.r", "TF_peak.fdr", "TF_gene.p_raw", "peak_gene.p_adj", "none"))
checkmate::assertSubset(xargs$ranking_TF_peak_gene, c("TF_gene.p_raw", "TF_gene.r", "TF_peak.fdr", "TF_gene.p_raw", "peak_gene.p_adj", "none"))

# Read in GRN object, because we have to calculate TF-gene correlations as well
GRN = qs::qread(xargs$input[1])

con.df = GRaNIE::getGRNConnections(GRN, include_TF_gene_correlations = TRUE)

checkmate::assertSubset(xargs$ranking_TF_gene, colnames(con.df))
checkmate::assertSubset(xargs$ranking_TF_peak_gene, colnames(con.df))

con.sel.df = con.df %>%
    dplyr::rename(source = "TF.ID", target = "gene.name", region = "peak.ID") %>%
    dplyr::mutate(weight = 1) %>%
    dplyr::select("source", "target", "region", "weight", tidyselect::everything())

# Select the user-chosen ranking column for both tables
TF_gene.df = con.sel.df %>%
    dplyr::select("source", "target", "weight", one_of(xargs$ranking_TF_gene))

TF_peak_gene.df = con.sel.df %>%
    dplyr::select("source", "target", "region", "weight", one_of(xargs$ranking_TF_peak_gene))


readr::write_csv(TF_gene.df, xargs$output[1])
readr::write_csv(TF_peak_gene.df, xargs$output[2])

    

    