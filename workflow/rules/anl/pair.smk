localrules: pair_realsim, pair_fakesim


rule pair_real_cor:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        pair='dts/{org}/{dname}pair/cases/{case}/mdata.h5mu',
        npair='dts/{org}/{dname}npair/cases/{case}/mdata.h5mu',
    output:
        cors='anl/pair/{org}.{dname}.{case}.real_corvals.csv',
        stat='anl/pair/{org}.{dname}.{case}.real_corsstat.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/anl/pair/real_cors.py \
        -a {input.pair} \
        -b {input.npair} \
        -c {output.cors} \
        -d {output.stat}
        """


rule pair_fake_stats:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        mdata='dts/{org}/{dname}pair/cases/{case}/mdata.h5mu',
        barmap='dts/{org}/fake{dname}pair/barmap.csv',
    output:
        knn='anl/pair/{org}.{dname}.{case}.fake_knn.csv',
        cor='anl/pair/{org}.{dname}.{case}.fake_cor.csv',
        prp='anl/pair/{org}.{dname}.{case}.fake_prp.csv',
    singularity:
        'workflow/envs/gretabench.sif'
    shell:
        """
        python workflow/scripts/anl/pair/fake_stats.py \
        -a {input.mdata} \
        -b {input.barmap} \
        -c {output.knn} \
        -d {output.cor} \
        -e {output.prp}
        """


rule pair_realsim:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        p='anl/topo/{org}.{dname}pair.{case}.sims_mult.csv',
        n='anl/topo/{org}.{dname}npair.{case}.sims_mult.csv',
    output: 'anl/pair/{org}.{dname}.{case}.pvsn.csv'
    shell:
        """
        python workflow/scripts/anl/pair/pairsim.py \
        -a {input.p} \
        -b {input.n} \
        -o {output}
        """


rule pair_fakesim:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        p='anl/topo/{org}.{dname}pair.{case}.sims_mult.csv',
        f='anl/topo/{org}.fake{dname}pair.{case}.sims_mult.csv',
    output: 'anl/pair/{org}.{dname}.{case}.pvsf.csv'
    shell:
        """
        python workflow/scripts/anl/pair/pairsim.py \
        -a {input.p} \
        -b {input.f} \
        -o {output}
        """


rule pair_real_qc:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        pair='dts/{org}/{dname}pair/cases/{case}/mdata.h5mu',
        npair='dts/{org}/{dname}npair/cases/{case}/mdata.h5mu',
    output:
        qc='anl/pair/{org}.{dname}.{case}.qc.csv',
        nc='anl/pair/{org}.{dname}.{case}.ncells.csv'
    shell:
        """
        python workflow/scripts/anl/pair/realqc.py {input.pair} {input.npair} {output.qc} {output.nc}
        """
