rule download_pbmc10k:
    output:
        matrix="resources/pbmc10k/filtered_feature_bc_matrix.h5",
        peak_annot="resources/pbmc10k/atac_peak_annotation.tsv",
        atac_frags="resources/pbmc10k/atac_fragments.tsv.gz",
        atac_index="resources/pbmc10k/atac_fragments.tsv.gz.tbi"
    params:
        matrix=config['data']['pbmc10k']['matrix'],
        peak_annot=config['data']['pbmc10k']['peak_annot'],
        atac_frags=config['data']['pbmc10k']['atac_frags'],
        atac_index=config['data']['pbmc10k']['atac_index']
    shell:
        """
        wget '{params.matrix}' -O {output.matrix}
        wget '{params.peak_annot}' -O {output.peak_annot}
        wget '{params.atac_frags}' -O {output.atac_frags}
        wget '{params.atac_index}' -O {output.atac_index}
        """

rule annotate_pbmc10k:
    input:
        dir="resources/pbmc10k/",
        obj="resources/pbmc10k/filtered_feature_bc_matrix.h5"
    output:
        plot="results/pbmc10k/annotated.pdf",
        mdata="resources/pbmc10k/annotated.h5mu"
    singularity:
        "envs/gretabench.sif"
    shell:
        "python workflow/scripts/annotate/pbmc10k.py -i {input.dir} -p {output.plot} -g {config.use_gpu} -o {output.mdata}"

