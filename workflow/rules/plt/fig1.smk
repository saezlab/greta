localrules: plt_stab, plt_sims, plt_AREG, plt_fig1


rule plt_stab:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        ovc='anl/stab/pitupair.ovc.csv',
        auc='anl/stab/pitupair.auc.csv'
    output:
        stab='plt/fig1/stab.pdf',
        cors='plt/fig1/cors.pdf',
    shell:
        """
        python workflow/scripts/plt/fig1/stab.py {input.ovc} {input.auc} {output.stab}
        python workflow/scripts/plt/fig1/cors.py {input.ovc} {output.cors}
        """


rule plt_sims:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sims='anl/topo/pitupair.all.sims_mult.csv',
        stats='anl/topo/pitupair.all.stats_mult.csv',
        tss=rules.tss_aggr.output,
        dst='anl/tss/pitupair.all.dist.csv',
        net='anl/topo/pitupair.all.inter.csv',
    output: 'plt/fig1/sims.pdf'
    shell:
        """
        python workflow/scripts/plt/fig1/sims.py \
        {input.sims} {input.stats} {input.tss} {input.dst} {input.net} {output}
        """


rule plt_AREG:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sims='anl/topo/pitupair.all.sims_mult.csv',
        gann='dbs/hg38/gen/ann/dictys/ann.bed',
    output: 'plt/fig1/links_AREG.pdf'
    params:
        gene='AREG',
        tfs=['FOSL1', 'FOSL2', 'JUNB', 'JUND'],
        wsize=250000
    shell:
        """
        python workflow/scripts/plt/fig1/links.py \
        -s {input.sims} \
        -g {params.gene} \
        -t {params.tfs} \
        -a {input.gann} \
        -w {params.wsize} \
        -o {output}
        """

rule plt_fig1:
    threads: 1
    input:
        stab='plt/fig1/stab.pdf',
        cors='plt/fig1/cors.pdf',
        sims='plt/fig1/sims.pdf',
        areg='plt/fig1/links_AREG.pdf'
    output: 'plt/fig1/fig1.pdf'
    shell:
        """
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input.stab} {input.cors} {input.sims} {input.areg}
        """
