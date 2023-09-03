rule extract_trajectory:
    input:
        "resources/{dataset}/annotated.h5mu"
    singularity:
        "envs/gretabench.sif"
    output:
        mdata="resources/{dataset}/{trajectory}/mdata.h5mu",
    params:
        celltypes=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celltypes'],
        n_hvg=lambda w: config[w.dataset]['trajectories'][w.trajectory]['n_hvg'],
        n_downsample=lambda w: config[w.dataset]['trajectories'][w.trajectory]['n_downsample']
    shell:
        "python workflow/scripts/extract_trajectory.py -i {input} -c '{params.celltypes}' -g {params.n_hvg} -d {params.n_downsample} -o {output.mdata}"

