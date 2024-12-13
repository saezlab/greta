# GRETA benchmark
Benchmark of GRN using the GRETA pipeline.

This pipeline uses `Snakemake` to ensure reproducibility.

<div align="center">
   <img src="https://drive.google.com/uc?id=1HpJx1deKivG2DRv3uXp_xLz90R0YfOwU" alt="GRETA graphical abstract" width="500" />
</div>

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

### 1. Install dependancies and run on any dataset
Then, you should install the necessary dependencies to run the new method. Here I would recommend creating a new conda enviroment for testing:
```
mamba create -n=test_method -c conda-forge -c bioconda \
python=3.10 \
# add packages here, if you need R also add r-base
```

Once you manage to make the method run and generate these two outputs, write a script in `workflow/scripts/` that accepts command line parameters reproducing the results and name it after the method. Try to identify all the relevant parameters of the method and add them as arguments.

### 2. Integration into the `Snakemake` pipeline
To add your code to the pipeline you will have to write a `Snakemake` rule file (`.smk`) in `workflow/rules/` with the name of the method. If this is unfamiliar to you, check [Snakemake's documentation](https://snakemake.readthedocs.io/en/stable/snakefiles/rules.html) and [`pando`'s example](https://github.com/saezlab/greta_benchmark/blob/main/workflow/rules/pando.smk).
For the rule definition, you can skip the `singularity` argument until you have not generated the required container, more info in the next section.
To check that it works run `Snakemake`, for example for `pando` it would be:
```
snakemake --profile config/slurm/ /path/to/output
```

### 3. Setting up a singularity container
To handle multiple dependencies across methods we need to build a singularity image container per method. Here are the instructions on [how to install singularity](https://apptainer.org/docs/user/latest/) in case you don't have it.

To create a singularity image file (`.sif`), first you need to manually write a definition file (`.def`) that contains all the instructions to install the necessary dependencies. Please do it inside `workflow/envs` and name it `namemethod.def`.

Here are two examples found in `workflow/envs/`, for [CellOracle](https://github.com/saezlab/greta_benchmark/blob/main/workflow/envs/celloracle.def) (which uses R and Python) and for [Pando](https://github.com/saezlab/greta_benchmark/blob/main/workflow/envs/pando.def) (which only uses R).

For a new method, check what are their dependencies and add them in the `.def` file. Try to add as many packages inside the `mamba` call, it is way faster than installing with `pip` or `devtools`. `R` packages are also available in conda, just keep in mind they might have a different name than their usual, a good strategy is to google this `r_package_name conda r`, and copy the designated name.

Once you think you have all possible dependencies, run the following command:

```
sudo singularity build workflow/envs/namemethod.sif workflow/envs/namemethod.def 
```

This will create the method's singularity container ready to be used. Now you can add this image to the `singularity` rule of your method.
