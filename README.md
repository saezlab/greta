# GRETA benchmark
Benchmark of GRN using the GRETA package

This pipeline uses `Snakemake` to ensure reproducibility.

## Installation
Clone repo:
```
git clone git@github.com:saezlab/greta_benchmark.git
cd greta_benchmark
```

Install `mamba` (this might take a while) to install packages faster:
```
conda install -n base -c conda-forge mamba
```

Then create a new enviroment specific for `Snakemake`:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
conda activate snakemake
```

