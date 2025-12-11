localrules: prt_knocktf, prt_knocktf_mm10


rule prt_knocktf_old:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 
        meta='dbs/hg38/prt/knocktf/meta.csv',
        diff='dbs/hg38/prt/knocktf/diff.csv',
        dir=directory('dbs/hg38/prt/knocktf/')
    params:
        url_m=config['dbs']['hg38']['prt']['knocktf']['meta'],
        url_d=config['dbs']['hg38']['prt']['knocktf']['diff'],
    shell:
        """
        wget --no-verbose '{params.url_m}' -O {output.meta} && \
        wget --no-verbose '{params.url_d}' -O {output.diff}
        """

rule prt_knocktf:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        meta='dbs/hg38/prt/knocktf/meta.csv.gz',
        diff='dbs/hg38/prt/knocktf/diff.csv.gz',
        dir=directory('dbs/hg38/prt/knocktf/')
    params: id=config['zenodo_id']
    shell:
        """
        wget --no-check-certificate --no-verbose \
        'https://zenodo.org/records/{params.id}/files/hg38_prt_knocktf_mat.csv.gz?download=1' \
        -O {output.diff} && \
        wget --no-check-certificate --no-verbose \
        'https://zenodo.org/records/{params.id}/files/hg38_prt_knocktf_meta.csv.gz?download=1' \
        -O {output.meta}
        """

rule prt_knocktf_mm10:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        meta='dbs/mm10/prt/knocktf/meta.csv.gz',
        diff='dbs/mm10/prt/knocktf/diff.csv.gz',
        dir=directory('dbs/mm10/prt/knocktf/')
    params: id=config['zenodo_id']
    shell:
        """
        wget --no-check-certificate --no-verbose \
        'https://zenodo.org/records/{params.id}/files/mm10_prt_knocktf_mat.csv.gz?download=1' \
        -O {output.diff} && \
        wget --no-check-certificate --no-verbose \
        'https://zenodo.org/records/{params.id}/files/mm10_prt_knocktf_meta.csv.gz?download=1' \
        -O {output.meta}
        """
