rule mdl_o_scdori:
    threads: 32
    conda: '../../envs/scdori.yaml'
    container: None
    input:
        mdata=rules.extract_case.output.mdata,
    output:
        out='dts/{dat}/cases/{case}/runs/o_scdori.o_scdori.o_scdori.o_scdori.mdl.csv'
    params:
        ext=config['methods']['scdori']['ext'] // 2,
    resources:
        partition='gpu-single',
        mem_mb=restart_mem,
        runtime=config['max_mins_per_step'] * 2,
        slurm="gres=gpu:1",
    shell:
        """
        path_out=$(dirname {output.out})
        path_out=$path_out/tmp_o_scdori
        mkdir -p $path_out
        python  workflow/scripts/mth/scdori/src.py \
        --config workflow/scripts/mth/scdori/config.yaml \
        window_size={params.ext} \
        base_dir=$path_out \
        data_dir=$path_out/data \
        weight_dir=$path_out/weights \
        mudata_file_name={input.mdata} \
        grn_file_out={output.out}
        """
