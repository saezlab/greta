rule extract_case:
    threads: 32
    input:
        lambda w: map_rules('annotate', w.dataset)
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        mdata='datasets/{dataset}/cases/{case}/mdata.h5mu',
    params:
        celltypes=lambda w: config['datasets'][w.dataset]['cases'][w.case]['celltypes'],
        n_sample=lambda w: config['datasets'][w.dataset]['cases'][w.case]['n_sample'] if 'n_sample' in config['datasets'][w.dataset]['cases'][w.case] else '0',
        seed=lambda w: config['datasets'][w.dataset]['cases'][w.case]['seed'] if 'n_sample' in config['datasets'][w.dataset]['cases'][w.case] else '0',
        n_hvg=lambda w: config['datasets'][w.dataset]['cases'][w.case]['n_hvg'],
        n_hvr=lambda w: config['datasets'][w.dataset]['cases'][w.case]['n_hvr'],
        root=lambda w: config['datasets'][w.dataset]['cases'][w.case]['root'] if 'root' in config['datasets'][w.dataset]['cases'][w.case] else 'None',
    shell:
        """
        python workflow/scripts/datasets/extract_case.py \
        -i '{input}' \
        -c '{params.celltypes}' \
        -s '{params.n_sample}' \
        -d '{params.seed}' \
        -g '{params.n_hvg}' \
        -r '{params.n_hvr}' \
        -t '{params.root}' \
        -o '{output.mdata}'
        """
