localrules: prt_knocktf


rule prt_knocktf:
    threads: 1
    output: 
        meta='dbs/hg38/prt/knocktf/meta.csv',
        diff='dbs/hg38/prt/knocktf/diff.csv',
    params:
        url_m=config['dbs']['hg38']['prt']['knocktf']['meta'],
        url_d=config['dbs']['hg38']['prt']['knocktf']['diff'],
    shell:
        """
        wget --no-verbose '{params.url_m}' -O {output.meta} && \
        wget --no-verbose '{params.url_d}' -O {output.diff}
        """
