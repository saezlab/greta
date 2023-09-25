library(argparse)
library(readr)
library(qs)

CURRENTLY NOT USED ANYMORE

parser <- ArgumentParser(description= 'Postprocess GRaNIE results')
parser$add_argument('--input', '-i', help= 'Output file from GRaNIE (filtered connections)', nargs = 2,  required= TRUE)
parser$add_argument('--output', '-o', help= 'Filenames for processed output files: 1. TF-gene, 2. TF-peak-gene', nargs= 2, required= TRUE)

xargs <- parser$parse_args()

# Read in GRN object, because we have to calculate TF-gene correlations as well
GRN = qs::qread(xargs$input[2])


    

    