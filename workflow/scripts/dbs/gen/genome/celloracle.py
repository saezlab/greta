from genomepy import install_genome
import os
import re
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-o','--orgms', required=True, nargs='+')
args = vars(parser.parse_args())

# Get dir
orgms = args['orgms']

# Install genomes
for path_org in orgms:
    org = re.search(r'^dbs/([^/]+)/.*$', path_org).group(1)
    install_genome(name=org, genomes_dir=path_org, provider="UCSC")
