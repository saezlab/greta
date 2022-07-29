rule peak_corr:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu",
        log="logs/add_r_env/celloracle.out"
    conda:
        "../envs/celloracle.yml"
    output:
        path_all_peaks="resources/{dataset}/{trajectory}/celloracle/all_peaks.csv",
        path_connections="resources/{dataset}/{trajectory}/celloracle/cicero_connections.csv",
        path_plot="results/{dataset}/{trajectory}/celloracle/peak_thr.pdf"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        min_count=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['min_count'],
        max_count=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['max_count']
    envmodules:
        "lib/openssl"
    shell:
        """
        Rscript workflow/scripts/celloracle/peak_corr.R {input.data} {params.organism} {params.min_count} {params.max_count} {output.path_plot} {output.path_all_peaks} {output.path_connections}
        """

rule tss_annotation:
    input:
        all_peaks="resources/{dataset}/{trajectory}/celloracle/all_peaks.csv",
        connections="resources/{dataset}/{trajectory}/celloracle/cicero_connections.csv"
    conda:
        "../envs/celloracle.yml"
    output:
        "resources/{dataset}/{trajectory}/celloracle/processed_peak_file.csv"
    params:
        organism=lambda w: config[w.dataset]['organism'],
        thr_coaccess=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_coaccess']
    shell:
         "python workflow/scripts/celloracle/tss_annotation.py -a {input.all_peaks} -c {input.connections} -o {params.organism} -t {params.thr_coaccess} -p {output}"

rule tf_motif_scan:
    input:
        "resources/{dataset}/{trajectory}/celloracle/processed_peak_file.csv"
    conda:
        "../envs/celloracle.yml"
    output:
        "resources/{dataset}/{trajectory}/celloracle/motifs.celloracle.tfinfo"
    resources:
        mem_mb=32000
    params:
        organism=lambda w: config[w.dataset]['organism'],
        fpr=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['fpr']
    shell:
        "python workflow/scripts/celloracle/tf_motif_scan.py -p {input} -o {params.organism} -f {params.fpr} -t {output}"

rule build_base_grn:
    input:
        "resources/{dataset}/{trajectory}/celloracle/motifs.celloracle.tfinfo"
    conda:
        "../envs/celloracle.yml"
    params:
        thr_motif_score=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_motif_score']
    output:
        "resources/{dataset}/{trajectory}/celloracle/base_GRN_dataframe.csv"
    shell:
        "python workflow/scripts/celloracle/build_base_grn.py -i {input} -t {params.thr_motif_score} -g {output}"

rule build_grn:
    input:
        mdata="resources/{dataset}/{trajectory}/mdata.h5mu",
        base_grn="resources/{dataset}/{trajectory}/celloracle/base_GRN_dataframe.csv"
    conda:
        "../envs/celloracle.yml"
    output:
        "resources/{dataset}/{trajectory}/celloracle/links.hdf5"
    shell:
        "python workflow/scripts/celloracle/build_grn.py -m {input.mdata} -b {input.base_grn} -l {output}"

rule filter_grn:
    input:
        "resources/{dataset}/{trajectory}/celloracle/links.hdf5"
    conda:
        "../envs/celloracle.yml"
    params:
        thr_edge_pval=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_edge_pval'],
        thr_top_edges=lambda w: config[w.dataset]['trajectories'][w.trajectory]['celloracle']['thr_top_edges']
    output:
        "resources/{dataset}/{trajectory}/celloracle/grns/celloracle_{celltype}.csv"
    shell:
        "python workflow/scripts/celloracle/filter_grn.py -l {input} -p {params.thr_edge_pval} -t {params.thr_top_edges} -o {output}"

