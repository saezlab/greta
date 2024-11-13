localrules: c2g_m_eqtlcatalogue, c2g_s_eqtlcatalogue



checkpoint c2g_m_eqtlcatalogue:
    threads: 1
    output: 'dbs/hg38/c2g/eqtlcatalogue/meta.tsv'
    params: url=config['dbs']['hg38']['c2g']['eqtlcatalogue']['meta']
    shell: 
        """
        wget --no-verbose '{params.url}' -O - | \
        awk -F'\t' '$9 == "ge" {{ print $1, $2, $6 }}' OFS='\t' > {output}
        """
        

checkpoint c2g_s_eqtlcatalogue:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_gid_ensmbl.output,
    output: 'dbs/hg38/c2g/eqtlcatalogue/raw/{eqtl_smpl_grp}.{eqtl_tiss}.bed',
    params:
        url=config['dbs']['hg38']['c2g']['eqtlcatalogue']['url'],
        thr_pval=config['dbs']['hg38']['c2g']['eqtlcatalogue']['thr_pval'],
    shell:
        """
        sleep $(shuf -i 1-5 -n 1)
        wget --no-verbose --tries 15 --wait=10 --random-wait --retry-connrefused --max-redirect=3 --waitretry=1 '{params.url}' -O - | \
        zcat | \
        python workflow/scripts/dbs/c2g/eqtlcat_smpl.py \
        {input} \
        '{params.thr_pval}' \
        {output}
        """


def eqtlcat_smpls(wildcards):
    meta_path = checkpoints.c2g_m_eqtlcatalogue.get().output[0]
    grps = []
    tisss = []
    with open(meta_path, 'r') as f:
        for line in f.readlines():
            grp, tiss = line.strip().split('\t')[:2]
            grps.append(grp)
            tisss.append(tiss)
    return expand('dbs/hg38/c2g/eqtlcatalogue/raw/{eqtl_smpl_grp}.{eqtl_tiss}.bed', zip, eqtl_smpl_grp=grps, eqtl_tiss=tisss)


checkpoint c2g_g_eqtlcatalogue:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input:
        smpls=eqtlcat_smpls,
        mta='dbs/hg38/c2g/eqtlcatalogue/meta.tsv'
    output: directory('dbs/hg38/c2g/eqtlcatalogue/raw/genes/'),
    shell:
        """
        mkdir -p {output} && \
        cat {input.smpls} | \
        python workflow/scripts/dbs/c2g/eqtlcat_gene.py \
        {input.mta} \
        {output}
        """


def eqtlcat_genes(wildcards):
    gdir = checkpoints.c2g_g_eqtlcatalogue.get().output[0]
    genes = glob_wildcards(os.path.join(gdir, "{gene}.bed")).gene
    return expand(os.path.join(gdir, "{gene}.bed"), gene=genes)



rule c2g_eqtlcatalogue:
    threads: 32
    singularity: 'workflow/envs/gretabench.sif'
    input: eqtlcat_genes
    output: 'dbs/hg38/c2g/eqtlcatalogue/eqtlcatalogue.bed'
    shell:
        """
        ls dbs/hg38/c2g/eqtlcatalogue/raw/genes/*.bed | xargs -n 1 -P {threads} -I {{}} sh -c \
        'sort -k1,1 -k2,2n "{{}}" | bedtools merge -i - -c 4,5 -o distinct,distinct > "{{}}.tmp"' && \
        cat dbs/hg38/c2g/eqtlcatalogue/raw/genes/*.bed.tmp > {output} && \
        rm dbs/hg38/c2g/eqtlcatalogue/raw/genes/*.bed.tmp
        """
