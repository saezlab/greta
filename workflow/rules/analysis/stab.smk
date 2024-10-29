def get_stab_paths(config, mthds, baselines, datasets):
    ns = [1024, 2048, 4096, 8192, 16384]
    seeds = [0, 1, 2]
    mthds = ['o_' + m for m in mthds if m != 'scenicplus']  # TODO: remove scenicplus filter when ready
    mthds.extend(baselines)
    d_lst = []
    c_lst = []
    m_lst = []
    for dataset in datasets:
        n_cell = ns[-1]
        for n_gene in ns:
            for seed in seeds:
                n_cre = n_gene * 4
                case = '{n_cell}_{n_gene}_{seed}'.format(n_cell=n_cell, n_gene=n_gene, seed=seed)
                config['datasets'][dataset]['cases'][case] = dict()
                config['datasets'][dataset]['cases'][case]['celltypes'] = 'all'
                config['datasets'][dataset]['cases'][case]['n_sample'] = n_cell
                config['datasets'][dataset]['cases'][case]['seed'] = seed
                config['datasets'][dataset]['cases'][case]['n_hvg'] = n_gene
                config['datasets'][dataset]['cases'][case]['n_hvr'] = n_cre
                for mth in mthds:
                    d_lst.append(dataset)
                    c_lst.append(case)
                    m_lst.append(mth)
        n_gene = ns[-1]
        n_cre = n_gene * 4
        for n_cell in ns:
            for seed in seeds:
                case = '{n_cell}_{n_gene}_{seed}'.format(n_cell=n_cell, n_gene=n_gene, seed=seed)
                config['datasets'][dataset]['cases'][case] = dict()
                config['datasets'][dataset]['cases'][case]['celltypes'] = 'all'
                config['datasets'][dataset]['cases'][case]['n_sample'] = n_cell
                config['datasets'][dataset]['cases'][case]['seed'] = seed
                config['datasets'][dataset]['cases'][case]['n_hvg'] = n_gene
                config['datasets'][dataset]['cases'][case]['n_hvr'] = n_cre
                for mth in mthds:
                    d_lst.append(dataset)
                    c_lst.append(case)
                    m_lst.append(mth)
        return d_lst, c_lst, m_lst


d_lst, c_lst, m_lst = get_stab_paths(config, mthds, baselines, stab_datasets)

rule run_stab:
    input:
        expand(['datasets/{dataset}/cases/{case}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'], zip, dataset=d_lst, case=c_lst, mth=m_lst)
    singularity: 'workflow/envs/gretabench.sif'
    output:
        tmp=temp('analysis/stab/tmp_{dataset}.csv'),
        res='analysis/stab/{dataset}.csv',
    params:
        m=mthds
    shell:
        """
        last_date=$(stat -c %y {input[0]} | cut -d ' ' -f 1)
        last_date=$(date -d '$last_date - 2 days' +%Y-%m-%d)
        sacct -S $last_date -E $(date -d '23:59:59 today' +%Y-%m-%dT%H:%M:%S) --state=COMPLETED --format=Jobname%100,elapsed,MaxRss,State | \
        awk '/^ +mdl_o_/ {{jobname=$1; getline; if ($1 == "batch") print jobname, $2, $3}}' > {output.tmp}
        python workflow/scripts/analysis/stab/run_stab.py \
        -i {output.tmp} \
        -o {output.res}
        """
