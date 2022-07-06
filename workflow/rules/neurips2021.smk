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

rule download_neurips2021:
    output:
        "resources/neurips2021/original.h5ad.gz"
    params:
        url=config['data']['neurips2021']
    shell:
        "wget '{params.url}' -O {output}"

rule decompress_neurips2021:
    input:
        "resources/neurips2021/original.h5ad.gz"
    output:
        "resources/neurips2021/original.h5ad"
    shell:
        "gzip -k -d {input}"
        
rule split_neurips2021:
    input:
        "resources/neurips2021/original.h5ad"
    output:
        gex="resources/neurips2021/gex.h5ad",
        atac="resources/neurips2021/atac.h5ad",
        plot="results/neurips2021/trajectory.pdf"
    conda:
        "../envs/gretabench.yml"
    shell:
        "python workflow/scripts/split_neurips2021.py -i {input} -p {output.plot} -g {output.gex} -a {output.atac}"
        
# snakemake -c8 --use-conda split_neurips2021
# conda env update --file local.yml --prune