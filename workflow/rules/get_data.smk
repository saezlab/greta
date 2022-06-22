# Downloads datasets from urls.
rule download_geo:
    output:
        "resources/{name}/original"
    params:
        url=lambda wildcards: config['data'][wildcards.name]
    shell:
        "wget '{params.url}' -O {output}"
        
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
        "gzip -d {input}"
        
rule split_neurips2021:
    input:
        "resources/neurips2021/original.h5ad"
    output:
        "resources/neurips2021/gex.h5ad"
        "resources/neurips2021/atac.h5ad"
    conda:
        "envs/celloracle.yml"
    shell:
        "python scripts/split_neurips2021.py -i {input} -o resources/neurips2021/"