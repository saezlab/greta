rule add_r_env:
    conda:
        "../envs/{env}.yml"
    output:
        "logs/add_r_env/{env}.out"
    shell:
        "Rscript workflow/envs/{wildcards.env}.R > {output}"

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

