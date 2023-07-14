
rule download_neurips2021:
    output:
        "resources/neurips2021/original.h5ad.gz"
    params:
        url=config['neurips2021']['data']
    shell:
        "wget '{params.url}' -O {output}"

rule decompress_neurips2021:
    input:
        "resources/neurips2021/original.h5ad.gz"
    output:
        "resources/neurips2021/original.h5ad"
    shell:
        "gzip -d {input}"
        
rule annotate_neurips2021:
    resources:
        mem_mb=48000,
    input: "resources/neurips2021/original.h5ad"
    output:
        mdata="resources/neurips2021/annotated.h5mu"
    conda:
        "../envs/gretabench.yml"
    shell:
        "python workflow/scripts/annotate/neurips2021.py -i {input} -o {output.mdata}"
