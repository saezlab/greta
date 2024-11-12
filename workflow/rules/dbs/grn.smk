localrules: grn_collectri, grn_dorothea


rule grn_collectri:
    threads: 1
    input: rules.gen_tfs_lambert.output
    output: 'dbs/hg38/grn/collectri.csv'
    params: url=config['dbs']['hg38']['grn']['collectri']
    run:
        import pandas as pd
        net = pd.read_csv(params.url).drop(columns=['resources', 'sign_decision'])
        net = net.drop(columns='references')
        net = net[['source', 'target', 'weight']]
        tfs = set(pd.read_csv(input[0], header=None).iloc[:, 0].astype('U'))
        net = net[net['source'].isin(tfs)]
        net.to_csv(output[0], index=False)


rule grn_dorothea:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.gen_tfs_lambert.output
    output: 'dbs/hg38/grn/dorothea.csv'
    params: url=config['dbs']['hg38']['grn']['dorothea']
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
