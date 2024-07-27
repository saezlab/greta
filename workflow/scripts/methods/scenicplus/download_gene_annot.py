import pybiomart as pbm
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('-h', '--hsapiens', required=True)
parser.add_argument('-m', '--mmusculus', required=True)
args = parser.parse_args()
organisms = {"hsapiens": args.hspaiens, "mmusculus": args.mmusculus}


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
        annot.to_csv(organisms[organism], overwrite=True)
