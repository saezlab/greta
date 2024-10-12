import subprocess
import os
import argparse

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d','--gdir', required=True)
args = vars(parser.parse_args())

to_ = args['gdir']


# Find Homer Executable path
result = subprocess.run("which homer", shell=True, text=True, stdout=subprocess.PIPE, stderr=subprocess.PIPE)
homer_path = result.stdout.strip()
from_ = os.path.dirname(os.path.dirname(os.path.realpath(homer_path)))

os.system(f'mkdir -p {to_}')
os.system(f'cp -R {from_}/data/genomes/hg38/* {to_}')
