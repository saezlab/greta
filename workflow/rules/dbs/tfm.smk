localrules: tfm_hpa


rule tfm_hpa:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.lambert.output[0]
    output: 'dbs/hg38/tfm/hpa/hpa.tsv'
    params:
        url=config['dbs']['hg38']['tfm']['hpa']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.zip && \
        python workflow/scripts/dbs/tfm/hpa.py \
        -i {output}.zip \
        -t {input} \
        -o {output}
        """
