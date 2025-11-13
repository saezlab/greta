localrules: cre_blacklist, cre_blacklist_mm10, cre_encode, cre_encode_mm10, cre_gwascatalogue, cre_phastcons, cre_promoters,cre_promoters_mm10, cre_zhang21


rule cre_blacklist:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/blacklist/blacklist.bed'
    params:
        url=config['dbs']['hg38']['cre']['blacklist']
    shell:
        """
        wget --no-verbose -O - "{params.url}" | zcat > {output}
        """

rule cre_blacklist_mm10:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output:
        'dbs/mm10/cre/blacklist/blacklist.bed'
    params:
        url=config['dbs']['mm10']['cre']['blacklist']
    shell:
        """
        mkdir -p $(dirname {output})
        wget --no-verbose {params.url} -O - | zcat > {output}
        """

rule cre_encode:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/encode/encode.bed'
    params:
        url=config['dbs']['hg38']['cre']['encode']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        cat {output}.tmp | sort -k 1,1 -k2,2n | bedtools merge -c 6 -o distinct > {output} && \
        rm {output}.tmp
        """

rule cre_encode_mm10:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/mm10/cre/encode/encode.bed'
    params:
        url=config['dbs']['mm10']['cre']['encode']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        cat {output}.tmp | sort -k 1,1 -k2,2n | bedtools merge -c 6 -o distinct > {output} && \
        rm {output}.tmp
        """

rule cre_gwascatalogue:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/gwascatalogue/gwascatalogue.bed'
    params:
        url=config['dbs']['hg38']['cre']['gwascatalogue']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output} && \
        python workflow/scripts/dbs/cre/gwascatalogue.py -i {output} && \
        sort -k 1,1 -k2,2n {output} | bedtools merge -i - -c 4,5 -o distinct,distinct -delim "|" > {output}.tmp && \
        mv {output}.tmp {output}
        """


rule cre_phastcons:
    threads: 1
    singularity: 'workflow/envs/pando.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/phastcons/phastcons.bed'
    params:
        url=config['dbs']['hg38']['cre']['phastcons']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        Rscript -e " \
        df <- get(load('{output}.tmp')); \
        df <- GenomicRanges::reduce(df); \
        df <- as.data.frame(df)[, c('seqnames', 'start', 'end')]; \
        write.table(x=df, file='{output}.tmp', sep = '\t', quote=FALSE, row.names=FALSE, col.names=FALSE)" && \
        sort -k 1,1 -k2,2n {output}.tmp > {output} && \
        rm {output}.tmp
        """


rule cre_promoters:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/promoters/promoters.bed'
    params:
        wsize=config['cre_prom_size']
    shell:
        """
        Rscript workflow/scripts/dbs/cre/promoters.R \
        {params.wsize} \
        {output}
        """

rule cre_promoters_mm10:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/mm10/cre/promoters/promoters.bed'
    params:
        wsize=config['cre_prom_size']
    shell:
        """
        Rscript workflow/scripts/dbs/cre/promoters_mm10.R \
        {params.wsize} \
        {output}
        """


rule cre_zhang21:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/cre/zhang21/zhang21.bed'
    params:
        url=config['dbs']['hg38']['cre']['zhang21']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        zcat {output}.tmp | bedtools merge > {output} && \
        rm {output}.tmp
        """
