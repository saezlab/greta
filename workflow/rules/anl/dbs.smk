rule dbs_stats:
    threads: 1
    input:
        paths_prt=expand('dbs/hg38/prt/{prt}/meta.csv', prt=config['dbs']['hg38']['prt'].keys()),
        paths_gst=expand('dbs/hg38/gst/{gst}.csv', gst=config['dbs']['hg38']['gst'].keys()),
        paths_tfm=expand('dbs/hg38/tfm/{tfm}/{tfm}.tsv', tfm=config['dbs']['hg38']['tfm'].keys()),
        paths_tfp=expand('dbs/hg38/tfp/{tfp}/{tfp}.tsv', tfp=config['dbs']['hg38']['tfp'].keys()),
        paths_tfb=expand('dbs/hg38/tfb/{tfb}/{tfb}.bed', tfb=config['dbs']['hg38']['tfb'].keys()),
        paths_cre=expand('dbs/hg38/cre/{cre}/{cre}.bed', cre=config['dbs']['hg38']['cre'].keys()),
        paths_c2g=expand('dbs/hg38/c2g/{c2g}/{c2g}.bed', c2g=config['dbs']['hg38']['c2g'].keys()),
    output: 'anl/dbs/stats.csv'
    resources:
        mem_mb=32000
    shell: 
        """
        python workflow/scripts/anl/dbs/stats.py \
        -p {input.paths_prt} \
        -g {input.paths_gst} \
        -m {input.paths_tfm} \
        -t {input.paths_tfp} \
        -b {input.paths_tfb} \
        -c {input.paths_cre} \
        -e {input.paths_c2g} \
        -o {output}
        """


rule dbs_terms:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 
        paths_prt=expand('dbs/hg38/prt/{prt}/meta.csv', prt=config['dbs']['hg38']['prt'].keys()),
        paths_tfm=expand('dbs/hg38/tfm/{tfm}/{tfm}.tsv', tfm=config['dbs']['hg38']['tfm'].keys()),
        paths_tfb=expand('dbs/hg38/tfb/{tfb}/{tfb}.bed', tfb=config['dbs']['hg38']['tfb'].keys()),
        paths_cre=expand('dbs/hg38/cre/{cre}/{cre}.bed', cre=config['dbs']['hg38']['cre'].keys()),
        paths_c2g=expand('dbs/hg38/c2g/{c2g}/{c2g}.bed', c2g=config['dbs']['hg38']['c2g'].keys()),
    output: 'anl/dbs/terms.csv'
    resources:
        mem_mb=64000
    shell:
        """
        python workflow/scripts/anl/dbs/terms.py -i {input} -o {output}
        """


rule dbs_ocoef:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'anl/dbs/stats.csv',
    output: 'anl/dbs/ocoef.csv',
    shell:
        """
        python workflow/scripts/anl/dbs/ocoef.py {output}
        """
