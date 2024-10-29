rule pair_real_cor:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        pair='datasets/{dname}pair/cases/{case}/mdata.h5mu',
        npair='datasets/{dname}npair/cases/{case}/mdata.h5mu',
    output:
        cors='analysis/pair/{dname}.{case}.real_corvals.csv',
        stat='analysis/pair/{dname}.{case}.real_corsstat.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/analysis/pair/real_cors.py \
        -a {input.pair} \
        -b {input.npair} \
        -c {output.cors} \
        -d {output.stat}
        """


rule pair_fake_stats:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata='datasets/{dname}pair/cases/{case}/mdata.h5mu',
        barmap='datasets/fake{dname}pair/barmap.csv',
    output:
        knn='analysis/pair/{dname}.{case}.fake_knn.csv',
        cor='analysis/pair/{dname}.{case}.fake_cor.csv',
        prp='analysis/pair/{dname}.{case}.fake_prp.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/analysis/pair/fake_stats.py \
        -a {input.mdata} \
        -b {input.barmap} \
        -c {output.knn} \
        -d {output.cor} \
        -e {output.prp}
        """


rule pair_sim:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
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
