rule pair_corr:
    threads: 1
    input:
        pair='datasets/{dname}pair/cases/{case}/mdata.h5mu',
        npair='datasets/{dname}npair/cases/{case}/mdata.h5mu',
    output:
        cors='analysis/pair/{dname}.{case}.cors_vals.csv',
        stat='analysis/pair/{dname}.{case}.cors_stat.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/analysis/pair/cors_and_stat.py \
        -a {input.pair} \
        -b {input.npair} \
        -c {output.cors} \
        -d {output.stat}
        """


rule pair_closek:
    threads: 1
    input:
        mdata='datasets/{dname}pair/cases/{case}/mdata.h5mu',
        barmap='datasets/fake{dname}pair/barmap.csv',
    output:
        ks='analysis/pair/{dname}.{case}.closek.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/analysis/pair/closek.py \
        -a {input.mdata} \
        -b {input.barmap} \
        -c {output.ks} \
        """


rule pair_sim:
    threads: 1
    input:
        p='analysis/topo/{dname}pair.{case}.sims_mult.csv',
        n='analysis/topo/{dname}npair.{case}.sims_mult.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    output:
        make_combs(
            path='analysis/pair/{dname}.{case}/',
            mthds=mthds,
            name='scores',
        )
    shell:
        """
        python workflow/scripts/analysis/pair/compute_pairsim.py -i {input.p}
        """