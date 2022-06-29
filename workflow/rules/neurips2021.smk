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