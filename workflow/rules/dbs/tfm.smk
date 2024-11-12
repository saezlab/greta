localrules: tfm_hpa, tfm_tfmdb


rule tfm_hpa:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_tfs_lambert.output
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


rule tfm_tfmdb:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/tfm/tfmdb/tfmdb.tsv'
    params:
        url=config['dbs']['hg38']['tfm']['tfmdb']
    shell:
        """
        wget --no-verbose '{params.url}' -O {output} && \
        python -c "import pandas as pd; \
        import sys; \
        df = pd.read_csv(sys.argv[1]); \
        df = df[['Gene Name', 'Cell Name', 'Tissue Type']]; \
        df['ctype'] = df['Cell Name'] + ',' + df['Tissue Type']; \
        df = df.groupby('Gene Name', as_index=False)['ctype'].apply(lambda x: ','.join(x)); \
        df['ctype'] = [','.join(sorted(set(s.split(',')))) for s in df['ctype']]; \
        df = df.drop_duplicates(['Gene Name', 'ctype']); \
        df = df.rename(columns={{'Gene Name': 'gene'}}); \
        df = df.sort_values(['gene', 'ctype']); \
        df.to_csv(sys.argv[1], sep='\\t', index=False, header=None)" {output}
        """
