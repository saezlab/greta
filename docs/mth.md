## Adding methods
A new method can be added to `greta` in two ways:
- Decoupled: The method is implemented to perform the `pre`, `p2g`, `tfb`, and `mdl` steps separately, allowing the effect of each step to be evaluated individually. This approach results in several GRNs.
- Baseline: Alternatively, the method can be added as a baseline without decoupling its inference steps. This approach results in a single GRN.

Assuming the name of the method is `newmethod`, a rule file for the method must be added in `workflow/rules/mth/newmethod.smk`. This file should include at least the following rules in a specified format.
For each rule, scripts running the specific steps should be placed in `workflow/scripts/mth/{method_name}/{step}.py`.
For reference, other `.smk` files and associated scripts can be consulted for the already implemented methods.
To ensure reproducibility, it is encouraged to methods' code inside a method-specific [singularity](https://docs.sylabs.io/guides/3.0/user-guide/quick_start.html) image, stored in `workflow/envs/{method_name}.sif`. Note that `.sif` files should never be pushed to a repository, instead only commit their `.def` file. Conda enviroment are also accepted, but only when singularity containers are not an option (for example as in `dictys`).
### Decoupled version
#### `pre` - preprocessing
```
rule pre_newmethod:  # Rule to preprocess the data
    threads: 16  # Change to a desired number of CPUs
    singularity: 'workflow/envs/newmethod.sif'  # Path to singularity image
    input: rules.extract_case.output.mdata   # Automated path to any input dataset
    output:
        out='dts/{dat}/cases/{case}/runs/mymethod.pre.h5mu'
    params:
        a=config['methods']['mymethod']['a'],  # Method parameters if needed
        b=config['methods']['mymethod']['b'],
    resources:  # Handles resource allocation in the cluster
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/mth/mymethod/pre.py ...  # Note the script's path
        """
```
The input is a `h5mu` dataset with specific properties (see Adding datasets). The output should be a `h5mu` object containing the processed dataset after applying the method's preprocessing strategy.

#### `p2g` - CRE (peak) to gene
```
rule p2g_newmethod:  # Rule to find the CRE-Gene links
    threads: 16
    singularity: 'workflow/envs/newmethod.sif'
    input:
        pre=lambda w: map_rules('pre', w.pre),  # Automated path to any other method pre step
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.newmethod.p2g.csv',
    # Here params if needed
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/mymethod/p2g.py ... '
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,gene,score,pval" }}' > {output.out}
        fi
        """
```
Input is the result of any other methods `pre` step. Output should be a `.csv` file with the following format:
```
cre,gene,score
chr1-100039944-100040444,SLC35A3,24
...
```
To ensure that the pipeline runs in a reasonable amount of time, if an individual job lasts more than `runtime` an empty result is returned.

#### `tfb` - TF to CRE binding prediction
```
rule tfb_newmethod: # Rule to find TF-CRE links
    threads: 16
    singularity: 'workflow/envs/newmethod.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),  # Automated path to any other method p2g step
    output: out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.newmethod.tfb.csv'
    # Here params if needed
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/mymethod/tfb.py ... '
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "cre,tf,score" }}' > {output.out}
        fi
        """
```
Inputs are the result of any other methods `pre` and `p2g` steps. Output should be a `.csv` file with the following format:
```
cre,tf,score
chr1-100039944-100040444,POU1F1,20.76
...
```
Note: here tfb scores should reflect some degree of significance, for example, most methods represent this score as `-log10(p-value)` of the event happening. This format is relevant later for the `mdl` step of some of the methods. 

#### `mdl` - Modeling
```
rule mdl_newmethod: # Rule to model final TF-Gene links
    threads: 16
    singularity: 'workflow/envs/newmethod.sif'
    input:
        pre=lambda wildcards: map_rules('pre', wildcards.pre),
        p2g=lambda wildcards: map_rules('p2g', wildcards.p2g),
        tfb=lambda wildcards: map_rules('tfb', wildcards.tfb),  # Automated path to any other method tfb step
    output:
        out='dts/{dat}/cases/{case}/runs/{pre}.{p2g}.{tfb}.newmethod.mdl.csv'
    # Here params if needed
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/mymethod/mdl.py ... '
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
```
Inputs are the result of any other methods `pre`, `p2g` and `tfb` steps. Output should be a `.csv` file with the following format:
```
source,target,score,pval
JUND,JUNB,3.24,0.01
...
```
Note: Here `pval` is optional but `score` is requiered.
Note: Benchmarking pipeline assumes that GRN methods perform their own edge prunning strategy, it should be added in this step if it is not available originaly.

#### `mdl_o_` - Original implementation
This rule implements the original method as it is in its source code, without decoupling it. This is requiered to make sure that the decoupled version is similar to the original one.
```
rule mdl_o_newmethod:  # Rule to infer the GRN as in the original implementation of the method
    threads: 16
    singularity: 'workflow/envs/newmethod.sif'
    input:
        mdata=rules.extract_case.output.mdata,
    output:
        out='dts/{dat}/cases/{case}/runs/o_newmethod.o_newmethod.o_newmethod.o_newmethod.mdl.csv',
    # Here params if needed
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,  # Increased time as source implementations can be slower
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/mymethod/src.py ... '
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
```
The input is a the same for the `pre` rule. Output should be a `.csv` file with the same format as `mdl`.

### Baseline
When adding a method as a baseline, only one rule is needed:
```
rule mdl_newmethod:
    threads: 16
    singularity: 'workflow/envs/newmethod.sif'
    input:
        mdata=rules.extract_case.output.mdata,
        ...
    output:
        out='dts/{dat}/cases/{case}/runs/newmethod.newmethod.newmethod.newmethod.mdl.csv'  # Important to repeat the name 4 times here
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        set +e
        timeout $(({resources.runtime}-20))m bash -c \
        'python workflow/scripts/mth/mymethod/mymethod.py ... '
        if [ $? -eq 124 ]; then
            awk 'BEGIN {{ print "source,target,score,pval" }}' > {output.out}
        fi
        """
```
Input should be a `h5mu` object. Output should be in the same format as in `mdl`.
