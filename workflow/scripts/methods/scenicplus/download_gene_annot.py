import pybiomart as pbm
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-j', '--hsapiens', required=True)
parser.add_argument('-m', '--mmusculus', required=True)
args = parser.parse_args()

organisms = {"hsapiens": args.hsapiens, "mmusculus": args.mmusculus}


for organism in organisms:
    if organism == 'hsapiens':
        dataset = pbm.Dataset(name='hsapiens_gene_ensembl',  host='http://www.ensembl.org')
    if organism == 'mmusculus':
        dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host='http://www.ensembl.org')

    annot = dataset.query(
        attributes=[
            'chromosome_name',
            'transcription_start_site',
            'strand', 'external_gene_name',
            'transcript_biotype'
            ])
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
    filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot = annot[annot.Transcript_type == 'protein_coding']
    annot["Strand"] = annot["Strand"].replace({1: "+", -1: "-"})
    annot.Start = annot.Start.astype(np.int32)
    annot.to_csv(organisms[organism], sep="\t", index=False)
