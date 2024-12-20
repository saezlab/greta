# Gene Regulatory nETwork Analysis (GRETA) <img src="https://drive.google.com/uc?id=1nx66ibhwjjnooCRGaxDvX27A3IT7EaPg" align="right" width="120" class="no-scaled-link" alt='GRETA logo' />

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
