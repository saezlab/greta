rule extract_trajectory:
    input:
        "resources/{dataset}/annotated.h5mu"
    conda:
        "../envs/gretabench.yml"
    output:
        mdata="resources/{dataset}/{trajectory}/mdata.h5mu",
        plot="results/{dataset}/{trajectory}/celloracle/{trajectory}.pdf"
    params:
        celltypes=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celltypes']
    shell:
        "python workflow/scripts/extract_trajectory.py -i {input} -c '{params.celltypes}' -p {output.plot} -o {output.mdata}"

