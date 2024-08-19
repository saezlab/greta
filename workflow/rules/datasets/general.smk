rule extract_case:
    input:
        'datasets/{dataset}/annotated.h5mu'
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

rule download_geneids:
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        hg='gdata/geneids/hg38.csv',
        mm='gdata/geneids/mm10.csv',
        dr=directory('gdata/geneids')
    shell:
        """
        Rscript workflow/scripts/datasets/download_geneids.R \
        {output.hg} \
        {output.mm}
        """
