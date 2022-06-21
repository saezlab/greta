# Downloads the multiome dataset from the NeurIPS 2021 competition.
rule download_NeurIPS2021:
    shell:
        "wget '{config[data][neurips2021]}' -O resources/neurips2021.h5ad"
