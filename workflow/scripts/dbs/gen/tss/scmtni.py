import os
import subprocess
import sys
import tempfile

import pandas as pd

path_inp = sys.argv[1]
path_chain = sys.argv[2]
path_out = sys.argv[3]

bed = pd.read_csv(path_inp, sep='\t', header=None)
bed = bed[[0, 1, 2, 6]].dropna().drop_duplicates()

with tempfile.TemporaryDirectory() as tmpdir:
    path_tmp_in = os.path.join(tmpdir, 'in.bed')
    path_tmp_out = os.path.join(tmpdir, 'out.bed')
    path_tmp_unmap = os.path.join(tmpdir, 'unmap.bed')

    bed.to_csv(path_tmp_in, sep='\t', index=False, header=None)

    subprocess.run(
        ['liftOver', path_tmp_in, path_chain, path_tmp_out, path_tmp_unmap],
        check=True
    )

    bed_hg38 = pd.read_csv(path_tmp_out, sep='\t', header=None)

bed_hg38.to_csv(path_out, sep='\t', index=False, header=None, compression='gzip')
