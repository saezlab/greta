rule mdl_o_scgpt:
    threads: 4
    conda: '../../envs/scgpt.yaml'
    input:
        mdata=rules.extract_case.output.mdata,
        tf=rules.gen_tfs_lambert.output,
        cg=rules.cre_promoters.output,
    output:
        out='dts/{org}/{dat}/cases/{case}/runs/o_scgpt.o_scgpt.o_scgpt.o_scgpt.mdl.csv'
    params:
        batch_size=config['methods']['scgpt']['batch_size'],
        device=config['methods']['scgpt']['device'],
        n_hvg=config['methods']['scgpt']['n_hvg'],
        min_score=config['methods']['scgpt']['min_score'],
        model_dir='dbs/hg38/gen/genome/scgpt',
    resources:
        partition='gpu-single',
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
        slurm="gres=gpu:1",
    shell:
        """
        export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:${{LD_LIBRARY_PATH:-}}"
        python workflow/scripts/mth/scgpt/src.py \
        --data {input.mdata} \
        --promoters {input.cg} \
        --batch_size {params.batch_size} \
        --tfs {input.tf} \
        --device {params.device} \
        --n_hvg {params.n_hvg} \
        --min_score {params.min_score} \
        --min_top_q 5 \
        --verbose \
        --model_dir {params.model_dir} \
        --out {output.out}
        """
