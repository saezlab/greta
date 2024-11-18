localrules: gst_collectri, gst_dorothea, gst_pthw, gst_prog


rule gst_collectri:
    threads: 1
    input: rules.gen_tfs_lambert.output
    output: 'dbs/hg38/gst/collectri.csv'
    params: url=config['dbs']['hg38']['pkn']['collectri']
    run:
        import pandas as pd
        net = pd.read_csv(params.url).drop(columns=['resources', 'sign_decision'])
        net = net.drop(columns='references')
        net = net[['source', 'target', 'weight']]
        tfs = set(pd.read_csv(input[0], header=None).iloc[:, 0].astype('U'))
        net = net[net['source'].isin(tfs)]
        net.to_csv(output[0], index=False)


rule gst_dorothea:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_tfs_lambert.output
    output: 'dbs/hg38/gst/dorothea.csv'
    params: url=config['dbs']['hg38']['pkn']['dorothea']
    shell:
        """
        Rscript -e '
        url <- "{params.url}"; \
        temp_file <- tempfile(); \
        download.file(url, destfile = temp_file, mode = "wb"); \
        load(temp_file); \
        unlink(temp_file); \
        df <- dorothea_hs; \
        df <- df[df$confidence %in% c("A", "B", "C"), ]; \
        df <- df[order(df$tf, df$target), ]; \
        colnames(df) <- c("source", "confidence", "target", "weight"); \
        df <- df[, c("source", "target", "weight")]; \
        tfs <- readLines("{input}"); \
        df <- df[df$source %in% tfs, ]; \
        write.csv(x=df, file="{output}", quote=FALSE, row.names=FALSE);'
        """


rule gst_pthw:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gst/{db}.csv',
    params: url=lambda w: config['dbs']['hg38']['gst'][w.db]
    shell:
        """
        wget --no-verbose '{params.url}' -O {output}.tmp && \
        python -c "import decoupler as dc; \
        gst = dc.read_gmt('{output}.tmp'); \
        gst['source'] = ['_'.join(s.split('_')[1:]) for s in gst['source']]; \
        gst.to_csv('{output}', index=False)" && \
        rm {output}.tmp
        """


rule gst_prog:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    output: 'dbs/hg38/gst/prog.csv',
    params: url=lambda w: config['dbs']['hg38']['gst']['prog']
    shell:
        """
        Rscript -e ' \
        download.file("{params.url}", destfile = "{output}.rda"); \
        load("{output}.rda"); \
        write.csv(model_human_full, "{output}.rda", row.names = FALSE, quote = FALSE);' && \
        python -c "import pandas as pd; \
        prg = pd.read_csv('{output}.rda'); \
        prg = prg.rename(columns={{'gene': 'target', 'pathway': 'source', 'p.value': 'pval'}}); \
        prg = prg[['source', 'target', 'weight', 'pval']]; \
        prg = prg[prg['pval'] < 1e-5]; \
        n = prg.groupby('source').size(); \
        prg = prg[prg['source'].isin(n[n > 5].index)]; \
        prg = prg.sort_values(['source', 'target', 'weight']); \
        prg.to_csv('{output}', index=False)" && \
        rm {output}.rda
        """
