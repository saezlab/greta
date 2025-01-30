localrules: download_intact, tfp_intact, tfp_europmc


rule download_intact:
    output:
        temp("dbs/hg38/tfp/intact/raw/intact.txt")
    params:
        url=config['dbs']['hg38']['tfp']['intact']
    shell:
        """
        wget --no-verbose {params.url} -O {output}.zip && \
        unzip -o {output}.zip -d $( dirname {output} ) && \
        rm {output}.zip
        """


rule tfp_intact:
    input:
        inc=rules.download_intact.output,
        lmb=rules.gen_tfs_lambert.output,
        pid=rules.gen_pid_uniprot.output,
    output: 'dbs/hg38/tfp/intact/intact.tsv'
    shell:
        """
        python workflow/scripts/dbs/tfp/intact.py \
        {input.inc} {input.lmb} {input.pid} {output}
        """


rule tfp_europmc_raw:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_tfs_lambert.output,
    output:
        single='dbs/hg38/tfp/europmc/raw/single.csv',
        pairs='dbs/hg38/tfp/europmc/raw/pairs.csv'
    params:
        min_chars=2,
        min_n=49
    resources:
        runtime=config['max_mins_per_step'] * 2,
    shell:
        """
        python workflow/scripts/dbs/tfp/europmc_raw.py \
        {input} {params.min_chars} {params.min_n} {output.single} {output.pairs}
        """


rule tfp_europmc:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        single=rules.tfp_europmc_raw.output.single,
        pairs=rules.tfp_europmc_raw.output.pairs,
    output: 'dbs/hg38/tfp/europmc/europmc.tsv'
    params:
        pval_thr=2.2e-16,
        min_odds=5,
    shell:
        """
        python workflow/scripts/dbs/tfp/europmc.py \
        {input.single} {input.pairs} {params.pval_thr} {params.min_odds} {output}
        """
