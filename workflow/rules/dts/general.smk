localrules: dts_gzip, compress_all

rule extract_case:
    threads: 16
    singularity: 'workflow/envs/gretabench.sif'
    input: lambda w: map_rules('annotate', w.dat)
    output:
        mdata='dts/{org}/{dat}/cases/{case}/mdata.h5mu',
    params:
        celltypes=lambda w: config['dts'][w.dat]['cases'][w.case]['celltypes'],
        n_sample=lambda w: config['dts'][w.dat]['cases'][w.case]['n_sample'] if 'n_sample' in config['dts'][w.dat]['cases'][w.case] else '0',
        seed=lambda w: config['dts'][w.dat]['cases'][w.case]['seed'] if 'n_sample' in config['dts'][w.dat]['cases'][w.case] else '0',
        n_hvg=lambda w: config['dts'][w.dat]['cases'][w.case]['n_hvg'],
        n_hvr=lambda w: config['dts'][w.dat]['cases'][w.case]['n_hvr'],
        root=lambda w: config['dts'][w.dat]['cases'][w.case]['root'] if 'root' in config['dts'][w.dat]['cases'][w.case] else 'None',
    shell:
        """
        python workflow/scripts/dts/extract_case.py \
        -i '{input}' \
        -c '{params.celltypes}' \
        -s '{params.n_sample}' \
        -d '{params.seed}' \
        -g '{params.n_hvg}' \
        -r '{params.n_hvr}' \
        -t '{params.root}' \
        -o '{output.mdata}'
        """

rule dts_gzip:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input: rules.extract_case.output.mdata,
    output: 'dts/compressed/{org}_dts_{dat}_{case}.h5mu.gz'
    shell:
        """
        python workflow/scripts/dts/compress.py {input} {output}
        """

rule compress_all:
    threads: 1
    singularity: 'workflow/envs/gretabench.sif'
    input:
        [
            'dts/compressed/hg38_dts_brain_all.h5mu.gz',
            'dts/compressed/hg38_dts_breast_all.h5mu.gz',
            'dts/compressed/hg38_dts_embryo_all.h5mu.gz',
            'dts/compressed/hg38_dts_eye_all.h5mu.gz',
            'dts/compressed/hg38_dts_heart_all.h5mu.gz',
            'dts/compressed/hg38_dts_kidney_all.h5mu.gz',
            'dts/compressed/hg38_dts_pbmc10k_all.h5mu.gz',
            'dts/compressed/hg38_dts_pitupair_all.h5mu.gz',
            'dts/compressed/hg38_dts_reprofibro_all.h5mu.gz',
            'dts/compressed/hg38_dts_skin_all.h5mu.gz',
            'dts/compressed/mm10_dts_epalate_all.h5mu.gz',
        ]
    output:
        'dts/compressed/done.txt'
    shell:
        """
        touch {output}
        """
