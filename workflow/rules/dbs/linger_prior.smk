import os

rule linger_prior:
    threads: 1
    singularity: 'workflow/envs/linger.sif'
    input: 'workflow/envs/linger.sif'
    output:
        archive='dbs/lingerGRN/data_bulk.tar.gz',
        dir=directory('dbs/lingerGRN/data_bulk')
    params:
        dir_name=lambda wildcards, output: os.path.dirname(output.archive)
    shell:
        """
        mkdir -p {params.dir_name}
        wget -nv --load-cookies /tmp/cookies.txt \
        "https://drive.usercontent.google.com/download?export=download&confirm=$(wget --quiet --save-cookies /tmp/cookies.txt --keep-session-cookies --no-check-certificate 'https://drive.usercontent.google.com/download?id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b'  -O- | sed -rn 's/.*confirm=([0-9A-Za-z_]+).*/\1/p')&id=1jwRgRHPJrKABOk7wImKONTtUupV7yJ9b" \
        -O '{output.archive}'
        rm -f /tmp/cookies.txt
        tar -xzf '{output.archive}' -C '{params.dir_name}'
        """