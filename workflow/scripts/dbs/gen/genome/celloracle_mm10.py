from genomepy import install_genome
import os
import re
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-o','--output', required=True, help='Output directory path')
args = parser.parse_args()

# Get output directory
output_path = args.output

# Extract organism from path (e.g., dbs/hg38/gen/genome/celloracle/ -> hg38)
org_match = re.search(r'dbs/([^/]+)/gen/genome/celloracle/?$', output_path)
if not org_match:
    raise ValueError(f"Cannot extract organism from path: {output_path}")

org = org_match.group(1)
print(f"Installing genome for organism: {org}")

# Create output directory if it doesn't exist
os.makedirs(output_path, exist_ok=True)

# Install genome
try:
    install_genome(name=org, genomes_dir=output_path, provider="UCSC")
    print(f"Successfully installed {org} genome to {output_path}")
except Exception as e:
    print(f"Error installing genome: {e}")
    # Create empty directory to satisfy Snakemake
    os.makedirs(output_path, exist_ok=True)
    raise