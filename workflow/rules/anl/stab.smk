localrules: run_stab, stab_ovsd, stab_cor, stab_unsmthds

def get_stab_paths(config, mthds, baselines, datasets):
    ns = [1024, 2048, 4096, 8192, 16384]
    seeds = [0, 1, 2]
    mthds = ['o_' + m for m in mthds]
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
                config['dts'][dataset]['cases'][case] = dict()
                config['dts'][dataset]['cases'][case]['celltypes'] = 'all'
                config['dts'][dataset]['cases'][case]['n_sample'] = n_cell
                config['dts'][dataset]['cases'][case]['seed'] = seed
                config['dts'][dataset]['cases'][case]['n_hvg'] = n_gene
                config['dts'][dataset]['cases'][case]['n_hvr'] = n_cre
                for mth in mthds:
                    d_lst.append(dataset)
                    c_lst.append(case)
                    m_lst.append(mth)
        n_gene = ns[-1]
        n_cre = n_gene * 4
        for n_cell in ns:
            for seed in seeds:
                case = '{n_cell}_{n_gene}_{seed}'.format(n_cell=n_cell, n_gene=n_gene, seed=seed)
                config['dts'][dataset]['cases'][case] = dict()
                config['dts'][dataset]['cases'][case]['celltypes'] = 'all'
                config['dts'][dataset]['cases'][case]['n_sample'] = n_cell
                config['dts'][dataset]['cases'][case]['seed'] = seed
                config['dts'][dataset]['cases'][case]['n_hvg'] = n_gene
                config['dts'][dataset]['cases'][case]['n_hvr'] = n_cre
                for mth in mthds:
                    d_lst.append(dataset)
                    c_lst.append(case)
                    m_lst.append(mth)
        return d_lst, c_lst, m_lst


d_lst, c_lst, m_lst = get_stab_paths(config, mthds, baselines, stab_datasets)

rule run_stab:
    threads: 1
    container: None
    input:
        expand(['dts/{dat}/cases/{case}/runs/{mth}.{mth}.{mth}.{mth}.grn.csv'], zip, dat=d_lst, case=c_lst, mth=m_lst)
    output:
        res='anl/stab/{dat}.ovc.csv',
        auc='anl/stab/{dat}.auc.csv',
    shell:
        """
        last_date=$(stat -c %y {input[0]} | cut -d ' ' -f 1)
        last_date=$(date -d "$last_date - 2 days" +%Y-%m-%d)
        sacct -S $last_date -E $(date -d '23:59:59 today' +%Y-%m-%dT%H:%M:%S) --state=COMPLETED --format=Jobname%100,elapsed,MaxRss,State | \
        awk '/^ +mdl_/ {{jobname=$1; getline; if ($1 == "batch") print jobname, $2, $3}}' > {output.res}.tmp &&
        python workflow/scripts/anl/stab/run_stab.py \
        -i {output.res}.tmp \
        -r {output.res} \
        -a {output.auc} &&
        rm {output.res}.tmp
        """


rule stab_cor:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.run_stab.output.res
    output:
        wgt='anl/stab/{dat}.wgt.csv',
        cor='anl/stab/{dat}.cor.csv',
    shell:
        """
        python workflow/scripts/anl/stab/seeds.py {input} {output.wgt} {output.cor}
        """


rule stab_ovsd:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.topo_mult.output.sims,
    output: 'anl/stab/{dat}.{case}.ovsd.csv'
    shell:
        """
        python workflow/scripts/anl/stab/ovsd.py {input} {output}
        """


for i in range(4):
    config['dts']['pbmc10k']['cases'][str(i)] = config['dts']['pbmc10k']['cases']['all'].copy()
    config['dts']['pbmc10k']['cases'][str(i)]['n_sample'] = 1000000
    config['dts']['pbmc10k']['cases'][str(i)]['seed'] = str(i)

rule stab_unsmthds:
    threads: 1
    container: None
    input: [[[os.path.join(os.path.dirname(p), '{dat}.' + str(case), f'{mth}.{mth}.{mth}.{mth}.scores.csv') for p in rules.metric_summ.input] for case in list(range(4)) + ['all']] for mth in ['dictys', 'scenicplus']]
    output: 'anl/stab/unsmthds/{dat}.scores.csv'
    shell:
        """
        python workflow/scripts/anl/metrics/aggregate.py \
        -i {input} \
        -o {output} \
        -a
        """
