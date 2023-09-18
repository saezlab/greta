
GRaNIE_singularity_path = "envs/GRaNIE.sif"

# TODO check if anything is needed in params

rule preprocess_GRaNIE:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    output:
        counts_rna  = "resources/{dataset}/{trajectory}/GRaNIE/counts_rna.tsv.gz",
        counts_atac = "resources/{dataset}/{trajectory}/GRaNIE/counts_atac.tsv.gz"
        metadata    = "resources/{dataset}/{trajectory}/GRaNIE/metadata.tsv.gz"
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.preprocess_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    params:
        organism=lambda w: config[w.dataset]['organism']
    script: "scripts/GRaNIE/preprocess.R"

rule run_GRaNIE:
    input:
        data=rules.preprocess_GRaNIE.output
    output:
        data = "resources/{dataset}/{trajectory}/GRaNIE/connections_filtered.tsv.gz",
        GRN  = "resources/{dataset}/{trajectory}/GRaNIE/GRN.qs"
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.run_GRaNIE.txt"
    singularity: GRaNIE_singularity_path
    params:
        organism=lambda w: config[w.dataset]['organism'],
    script: "scripts/GRaNIE/runGRaNIE.R"

rule postprocess_GRaNIE:
    input:
        data=rules.run_GRaNIE.output
    output:
        grn="resources/{dataset}/{trajectory}/GRaNIE/grn.csv",
        tri="resources/{dataset}/{trajectory}/GRaNIE/tri.csv"
    singularity: GRaNIE_singularity_path
    benchmark:
        "benchmarks/GRaNIE/{dataset}.{trajectory}.postprocess_GRaNIE.txt"
    params:
        p_thresh=lambda w: config[w.dataset]['trajectories'][w.trajectory]['pando']['p_thresh'],
    script: "scripts/GRaNIE/postprocess.R"
