import mudata as mu
import gzip
import shutil
import os
import sys
path_inp = sys.argv[1]
path_out_gz = sys.argv[2] 
mdata = mu.read(path_inp)
rna = mdata.mod['rna']
rna.X = rna.layers['counts'].astype(float)
del rna.layers
atac = mdata.mod['atac']
atac.X = atac.layers['counts'].astype(float)
del atac.layers
path_out = path_out_gz.replace('.gz', '')
mdata.write(path_out)
with open(path_out, 'rb') as f_in:
    with gzip.open(path_out_gz, 'wb', compresslevel=6) as f_out:
        shutil.copyfileobj(f_in, f_out)
os.remove(path_out)
