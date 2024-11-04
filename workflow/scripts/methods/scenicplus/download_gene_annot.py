from scenicplus.cli.commands import download_gene_annotation_chromsizes
import os
import argparse
import pandas as pd

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--gann', required=True, nargs='+')
parser.add_argument('-c', '--cist', required=True, nargs='+')
parser.add_argument('-s', '--csiz', required=True, nargs='+')
args = parser.parse_args()


for gann, cist, csiz in zip(args.gann, args.cist, args.csiz):
    org = os.path.basename(gann).split('_')[0]
    if org == 'hg38':
        species = "hsapiens"
    elif org == 'mm10':
        species = 'mmusculus'
    else:
        raise ValueError(f'Organism {org} not implemented.')

    download_gene_annotation_chromsizes(
        species=species,
        genome_annotation_out_fname=gann,
        chromsizes_out_fname=csiz,
        biomart_host="http://ensembl.org/",
        use_ucsc_chromosome_style=True
    )
    df = pd.read_csv(gann, sep="\t")
    df['Start'] = df['Start'] - 1
    df['Score'] = '.'
    df = df[['Chromosome', 'Start', 'End', 'Gene', 'Score', 'Strand', 'Transcript_type']]
    with open(cist, "w") as f:
        f.write("# " + "\t".join(df.columns) + "\n")
        for i, row in df.iterrows():
            text = '{c}\t{s}\t{e}\t{g}\t{r}\t{t}\n'.format(
                c=row['Chromosome'],
                s=row['Start'],
                e=row['End'],
                g=row['Gene'],
                r=row['Score'],
                t=row['Strand']
            )
            f.write(text)
