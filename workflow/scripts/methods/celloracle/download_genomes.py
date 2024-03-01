from genomepy import install_genome
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--gdir', required=True)
args = vars(parser.parse_args())

# Get dir
gdir = args['gdir']

# Install genomes
install_genome(name='hg38', genomes_dir=gdir, provider="UCSC")
install_genome(name='mm10', genomes_dir=gdir, provider="UCSC")
