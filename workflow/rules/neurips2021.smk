
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
        
rule annotate_neurips2021:
    resources:
        mem_mb=48000,
    input: "resources/neurips2021/original.h5ad"
    output:
        plot="results/neurips2021/annotated.pdf",
        mdata="resources/neurips2021/annotated.h5mu"
    conda:
        "../envs/gretabench.yml"
    shell:
        "python workflow/scripts/annotate/neurips2021.py -i {input} -p {output.plot} -o {output.mdata}"

