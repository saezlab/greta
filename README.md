# Gene Regulatory nETwork Analysis (GRETA) <img src="https://drive.google.com/uc?id=1DFGeAuSp8w1kDlMaS4zyeXfepKVW14Ym" align="right" width="120" class="no-scaled-link" alt='GRETA logo' />

<!-- badges: start -->
[![GitHub issues](https://img.shields.io/github/issues/saezlab/greta.svg)](https://github.com/saezlab/greta/issues/)
<!-- badges: end -->

Cells regulate their functions through gene expression, driven by a complex interplay of transcription factors and other regulatory mechanisms that together can be modeled as gene regulatory networks (GRNs).
The emergence of single-cell multi-omics technologies has driven the development of several methods that integrate transcriptomics and chromatin accessibility data to infer GRNs.
**Gene Regulatory nETwork Analysis (GRETA)** is a `Snakemake` pipeline that implements state-of-the-art multimodal GRN inference methods. It organizes the steps of these methods into a modular framework, enabling users to infer, compare, and benchmark GRN approaches.

<div align="center">
   <img src="https://drive.google.com/uc?id=1HpJx1deKivG2DRv3uXp_xLz90R0YfOwU" alt="GRETA graphical abstract" width="500" style="pointer-events: none;" />
</div>

## Installation
Clone repo:
```
git clone git@github.com:saezlab/greta.git
cd greta
```

Then create a new enviroment specific for `Snakemake`:
```
mamba create -c conda-forge -c bioconda -n snakemake snakemake
mamba activate snakemake
```

## Overview
Due to the magnitude of datasets and analyses, the repository is organized as a reproducible [Snakemake](https://snakemake.readthedocs.io/en/stable/) pipeline and uses [singularity](https://docs.sylabs.io/guides/3.5/user-guide/introduction.html) images to handle dependencies:
```
greta/
├── config/
│   ├── slurm/            # Cluster configuration (assumes Slurm architecture)
│   ├── config.yaml       # Specifies methods, datasets, and databases
│   └── prior_cats.json   # Specifies database labels for each dataset
└── workflow/
    ├── envs/             # Singularity definition (.def) and image (.sif) files
    ├── rules/            # Snakemake rules for:
    │   ├── anl              # analyses
    │   ├── dbs              # databases
    │   ├── dts              # datasets
    │   ├── mth              # methods
    │   └── plt              # plots
    ├── scripts/          # Helper scripts for:
    │   ├── anl              # analyses
    │   ├── dbs              # databases
    │   ├── dts              # datasets
    │   ├── mth              # methods
    │   └── plt              # plots
    └── Snakefile         # Main Snakemake file
```

Here are some lines to generate important intermediate outputs:
```
# Downloads and processes a dataset, for example pbmc10k
snakemake --profile config/slurm/ dts/pbmc10k/cases/all/mdata.h5mu

# Computes Pando's preprocessing step on the pbmc10k dataset
snakemake --profile config/slurm/ dts/pbmc10k/cases/all/runs/pando.pre.h5mu

# Computes GRaNIE's p2g step on Pando's pre
snakemake --profile config/slurm/ dts/pbmc10k/cases/all/runs/pando.granie.p2g.csv

# Computes CellOracles's tfb step on GRaNIE's p2g
snakemake --profile config/slurm/ dts/pbmc10k/cases/all/runs/pando.granie.celloracle.tfb.csv

# Computes Dictys's mdl step on the previous results
snakemake --profile config/slurm/ dts/pbmc10k/cases/all/runs/pando.granie.celloracle.dictys.mdl.csv

# Runs all possible method combinations, baselines and original implementations
snakemake --profile config/slurm/ anl/topo/pbmc10k.all.sims_mult.csv

# Downloads and processess all databases
snakemake --profile config/slurm/ anl/dbs/stats.csv

# Runs the mechanistic metric forecasting (perturbation) for all method combinations
snakemake --profile config/slurm/ anl/metrics/mech/prt/knocktf/pbmc10k.all.scores.csv

# Runs the benchmark for all databases and metrics
snakemake --profile config/slurm/ anl/metrics/pbmc10k.all.csv
```

## How to
- [add methods](docs/mth.md)
- [add datasets](docs/dts.md)
- [add databases](docs/dbs.md)

## Citation
Badia-i-Mompel et al. Comparison and evaluation of methods to infer gene regulatory networks from multimodal single-cell data. *bioRxiv* (2024) [doi:10.1101/2024.12.20.629764](https://doi.org/10.1101/2024.12.20.629764)
