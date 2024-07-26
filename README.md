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
## Adding new GRN inference methods
### 1. Download test data
The first thing is to download the toy dataset. It is a downsampled `MuData` object of an erythrocyte trajectory from the NeurIPS2021 competition. You can do this by running:
```
mkdir -p resources/neurips2021/small/
wget -O "resources/neurips2021/small/mdata.h5mu" "https://docs.google.com/uc?export=download&id=1ens89w191zDX1y444cf2TFzBZUOPNycu"
```

### 2. Install dependancies and run toy dataset
Then, you should install the necessary dependencies to run the new method. Here I would recommend creating a new conda enviroment for testing:
```
mamba create -n=test_method -c conda-forge -c bioconda \
python=3.10 \
# add packages here, if you need R also add r-base==4.2
```

Once installed, try running the method according to their vignette. For R, check [Pando's code](https://github.com/saezlab/greta_benchmark/blob/main/workflow/scripts/pando/run_pando.R#L33) on how to read `MuData` objects in R.
For the final output, we expect 2 things whenever possible:
1) A `grn.csv` data frame containing the classic TF-Gene networks with the following format:
   ```
   source target weight pval (other)
      TF1     G1    2.5  0.02       X
      TF1     G2   -1.7  0.05       X
      TF2     G1    1.2  0.01       X
      ...    ...    ...   ...     ...
   ```
   It needs to have at least a `source` (the TF name) and a `target` (the gene name, use name and not ID if possible), but other attributes might be included depending on the method. For ranking, we either need a `weight` column and/or a column for statistical significance (either`pval` or `padj`). Please stick to these names so that they are unified across methods. If your method doesn't use any weights, simply put all values to 1 or skip the column altogether. 
2) A `tri.csv` data frame containing the TF-CRE-Gene triplets. Some methods do not directly provide this information but it can be extracted for most. Check [Pando's code](https://github.com/saezlab/greta_benchmark/blob/main/workflow/scripts/pando/run_pando.R#L145) for an example, note that we filter triplets based on if they match TF-Genes from `grn.csv`.
   ```
   source target   region weight pval (other)
      TF1     G1 ch1-1-10    2.5  0.02       X
      TF1     G2 ch3-3-27   -1.7  0.05       X
      TF2     G1 ch9-8-42    1.2  0.01       X
      ...    ...      ...    ...   ...     ...
   ```
   Here, the same structure is needed as summarized above, with an additional `region` column (the CRE). Make sure the `region` columns uses only `-` as separator and not `:` in addition.

Once you manage to make the method run and generate these two outputs, write a script in `workflow/scripts/` that accepts command line parameters reproducing the results and name it after the method. Try to identify all the relevant parameters of the method and add them as arguments.

### 3. Integration into the `Snakemake` pipeline
To add your code to the pipeline you will have to write a `Snakemake` rule file (`.smk`) in `workflow/rules/` with the name of the method. If this is unfamiliar to you, check [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) and [Pando's example](https://github.com/saezlab/greta_benchmark/blob/main/workflow/rules/pando.smk).
Also add all the relevant parameters to the `small` dataset in the `config/config.yaml` [file](https://github.com/saezlab/greta_benchmark/blob/main/config/config.yaml) following the same format shown in `pando` and `celloracle`.
For the rule definition, you can skip the `singularity` argument until you have not generated the required container, more info in the next section.
To check that it works run `Snakemake`, for example for Pando it would be:
```
snakemake resources/neurips2021/small/pando/tri.csv
```

### 4. Setting up a singularity container
To handle multiple dependencies across methods we need to build a singularity image container per method. Here are the instructions on [how to install singularity](https://apptainer.org/docs/user/latest/) in case you don't have it.

To create a singularity image file (`.sif`), first you need to manually write a definition file (`.def`) that contains all the instructions to install the necessary dependencies. Please do it inside `workflow/envs` and name it `namemethod.def`.

Here are two examples found in `workflow/envs/`, for [CellOracle](https://github.com/saezlab/greta_benchmark/blob/main/workflow/envs/celloracle.def) (which uses R and Python) and for [Pando](https://github.com/saezlab/greta_benchmark/blob/main/workflow/envs/pando.def) (which only uses R).

For a new method, check what are their dependencies and add them in the `.def` file. Try to add as many packages inside the `mamba` call, it is way faster than installing with `pip` or `devtools`. `R` packages are also available in conda, just keep in mind they might have a different name than their usual, a good strategy is to google this `r_package_name conda r`, and copy the designated name.

Once you think you have all possible dependencies, run the following command:

```
sudo singularity build workflow/envs/namemethod.sif workflow/envs/namemethod.def 
```

This will create the method's singularity container ready to be used. Now you can add this image to the `singularity` rule of your method.
