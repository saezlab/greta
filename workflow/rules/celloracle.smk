rule peak_corr:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.peak_corr.txt"
    output:
        path_all_peaks="resources/{dataset}/{trajectory}/celloracle/all_peaks.csv",
        path_connections="resources/{dataset}/{trajectory}/celloracle/cicero_connections.csv"
    params:
        organism=lambda w: config[w.dataset]['organism']
    shell:
        """
        Rscript scripts/celloracle/peak_corr.R {input.data} {params.organism} {output.path_all_peaks} {output.path_connections}
        """

rule tss_annotation:
    input:
        all_peaks="resources/{dataset}/{trajectory}/celloracle/all_peaks.csv",
        connections="resources/{dataset}/{trajectory}/celloracle/cicero_connections.csv"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.tss_annotation.txt"
    output:
        "resources/{dataset}/{trajectory}/celloracle/processed_peak_file.csv"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        thr_coaccess=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_coaccess']
    shell:
         "python scripts/celloracle/tss_annotation.py -a {input.all_peaks} -c {input.connections} -o {params.organism} -t {params.thr_coaccess} -p {output}"

rule tf_motif_scan:
    input:
        "resources/{dataset}/{trajectory}/celloracle/processed_peak_file.csv"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.tf_motif_scan.txt"
    output:
        "resources/{dataset}/{trajectory}/celloracle/motifs.celloracle.tfinfo"
    resources:
        mem_mb=32000
    params:
        organism=lambda w: config[w.dataset]['organism'],
        fpr=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['fpr']
    shell:
        "python scripts/celloracle/tf_motif_scan.py -p {input} -o {params.organism} -f {params.fpr} -t {output}"

rule build_base_grn:
    input:
        "resources/{dataset}/{trajectory}/celloracle/motifs.celloracle.tfinfo"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.build_base_grn.txt"
    params:
        thr_motif_score=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_motif_score']
    output:
        "resources/{dataset}/{trajectory}/celloracle/base_GRN_dataframe.csv"
    shell:
        "python scripts/celloracle/build_base_grn.py -i {input} -t {params.thr_motif_score} -g {output}"

rule build_grn:
    input:
        mdata="resources/{dataset}/{trajectory}/mdata.h5mu",
        base_grn="resources/{dataset}/{trajectory}/celloracle/base_GRN_dataframe.csv"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.build_grn.txt"
    output:
        "resources/{dataset}/{trajectory}/celloracle/grn.celloracle.links"
    shell:
        "python scripts/celloracle/build_grn.py -m {input.mdata} -b {input.base_grn} -l {output}"

rule filter_grn:
    input:
        grn="resources/{dataset}/{trajectory}/celloracle/grn.celloracle.links",
        base="resources/{dataset}/{trajectory}/celloracle/base_GRN_dataframe.csv"
    singularity:
        "envs/celloracle.sif"
    benchmark:
        "benchmarks/celloracle/{dataset}.{trajectory}.filter_grn.txt"
    params:
        thr_edge_pval=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_edge_pval'],
        thr_top_edges=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_top_edges']
    output:
        grn="resources/{dataset}/{trajectory}/celloracle/grn.csv",
        base="resources/{dataset}/{trajectory}/celloracle/tri.csv"
    shell:
        "python scripts/celloracle/filter_grn.py -l {input.grn} -b {input.base} -p {params.thr_edge_pval} -t {params.thr_top_edges} -g {output.grn} -r {output.base}"
