rule run_pando:
    input:
        data="resources/{dataset}/{trajectory}/mdata.h5mu",
        log="logs/add_r_env/pando.out"
    conda:
        "../envs/pando.yml"
    output:
        grn="resources/{dataset}/{trajectory}/pando/grn.csv",
        gof="results/{dataset}/{trajectory}/pando/goodness_of_fit.pdf",
        modules="results/{dataset}/{trajectory}/pando/mod_distrs.pdf",
        topo="results/{dataset}/{trajectory}/pando/topology.pdf"
    params:
        organism=lambda w: config[w.dataset]['organism']
    envmodules:
        "lib/openssl"
    shell:
        """
        Rscript workflow/scripts/pando/run_pando.R {input.data} {params.organism} {output.gof} {output.modules} {output.topo} {output.grn}
        """
