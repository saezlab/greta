localrules: plt_dwns, plt_sims, plt_AREG, fig_stability


rule plt_dwns:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        ovc='anl/stab/pitupair.ovc.csv',
        auc='anl/stab/pitupair.auc.csv',
        wgt='anl/stab/pitupair.wgt.csv',
        cor='anl/stab/pitupair.cor.csv',
    output:
        stab='plt/stab/dwns.pdf',
        cors='plt/stab/cors.pdf',
    shell:
        """
        python workflow/scripts/plt/stab/stab.py {input.ovc} {input.auc} {output.stab}
        python workflow/scripts/plt/stab/cors.py {input.wgt} {input.cor} {output.cors}
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
    output: 'plt/stab/sims.pdf'
    shell:
        """
        python workflow/scripts/plt/stab/sims.py \
        {input.sims} {input.stats} {input.tss} {input.dst} {input.net} {output}
        """


rule plt_AREG:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        sims='anl/topo/pitupair.all.sims_mult.csv',
        gann='dbs/hg38/gen/ann/dictys/ann.bed',
    output: 'plt/stab/links_AREG.pdf'
    params:
        gene='AREG',
        tfs=['FOSL1', 'FOSL2', 'JUNB'],
        wsize=250000
    shell:
        """
        python workflow/scripts/plt/stab/links.py \
        -s {input.sims} \
        -g {params.gene} \
        -t {params.tfs} \
        -a {input.gann} \
        -w {params.wsize} \
        -o {output}
        """


rule fig_stability:
    threads: 1
    input:
        stab='plt/stab/dwns.pdf',
        cors='plt/stab/cors.pdf',
        sims='plt/stab/sims.pdf',
        areg='plt/stab/links_AREG.pdf'
    output: 'plt/stab/fig.pdf'
    shell:
        """
        gs -dBATCH -dNOPAUSE -q -sDEVICE=pdfwrite -sOutputFile={output} {input.stab} {input.cors} {input.sims} {input.areg}
        """
