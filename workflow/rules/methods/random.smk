localrules: lambert, chrom_gene


rule lambert:
    output: 'gdata/tfs/lambert.csv'
    params:
        url='http://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt'
    resources:
        mem_mb=2000
    shell:
        "wget '{params.url}' -O {output.out}"

rule chrom_gene:
    input: rules.download_granges.output.h
    output: 'gdata/chrom_gene/chrom_gene.csv'
    run:
        import pandas as pd
        ref = pd.read_csv(input[0], dtype={'seqnames': str})
        ref = ref.drop_duplicates(['seqnames', 'gene_name']).drop_duplicates('gene_name', keep=False)[['seqnames', 'gene_name']]
        ref = ref[['seqnames', 'gene_name']].sort_values(['seqnames', 'gene_name'])
        ref['chrom'] = 'chr' + ref['seqnames']
        ref = ref[['chrom', 'gene_name']].rename(columns={'gene_name': 'gene'})
        ref.to_csv(output[0], index=False)


rule mdl_random:
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.lambert.output[0],
        cg=rules.chrom_gene.output[0]
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        out='datasets/{dataset}/cases/{case}/runs/random.random.random.random.mdl.csv'
    params:
        g_perc=0.25,
        scale=1,
        tf_g_ratio=0.10,
        seed=42,
    resources:
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'],
    shell:
        """
        python workflow/scripts/methods/random/grn.py \
        -i {input.mdata} \
        -t {input.tf} \
        -c {input.cg} \
        -g {params.g_perc} \
        -n {params.scale} \
        -r {params.tf_g_ratio} \
        -s {params.seed} \
        -o {output.out}
        """
