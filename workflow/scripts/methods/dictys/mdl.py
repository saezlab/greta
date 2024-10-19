import argparse, os, sys
import pandas as pd
import numpy as np
import mudata as md
import dictys


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--path_data', required=True)
parser.add_argument('-g', '--p2g_data', required=True)
parser.add_argument('-b', '--tfb_data', required=True)
parser.add_argument('-p', '--path_out', required=True)
parser.add_argument('-v', '--device', required=True)
parser.add_argument('-t', '--threads', required=True)
parser.add_argument('-e', '--ext', required=True)
parser.add_argument('-a', '--gene_annotation', required=True)

args = vars(parser.parse_args())
mudata = args['path_data']
p2g_data = args['p2g_data']
tfb_data = args['tfb_data']
path_out = args['path_out']
device = f'cuda:{args["device"]}'
threads = int(args['threads'])
distance = int(args['ext'])
annot = args['gene_annotation']

working_dir = os.path.dirname(tfb_data)
pre = os.path.basename(tfb_data).split('.')[0]
p2g = os.path.basename(tfb_data).split('.')[1]
tfb = os.path.basename(tfb_data).split('.')[2]
tmp_path = f"{pre}.{p2g}.{tfb}.dictys_tmp"
os.makedir(tmp_path)

# QC Expression matrix
rna_filename = os.path.join(working_dir, tmp_path, "expression.tsv.gz")
proc_rna_filename = os.path.join(working_dir, tmp_path, "expression_proc.tsv.gz")

data = md.read(mudata)
rna_X = pd.DataFrame(np.array(data['rna'].layers['counts'].todense()).T, columns=data['rna'].obs.index, index=data['rna'].var.index)
if not os.path.exists(rna_filename):
    rna_X.to_csv(rna_filename, sep="\t", compression="gzip")
dictys.preproc.qc_reads(rna_filename, proc_rna_filename, 50, 10, 0, 0, 0, 0)  # all cell filters are turned off
rna_X_proc = pd.read_csv(proc_rna_filename, sep='\t', index_col=0)
genes = rna_X_proc.index.to_numpy()


# Write ATAC to working directory
atac_filename = os.path.join(working_dir, tmp_path, "atac_peak.tsv.gz")
atac_peak_names = [n.replace('-', ':') for n in data['atac'].var.index]
atac_X = pd.DataFrame(np.zeros((data['atac'].var.index.shape[0], 1)), index=atac_peak_names, columns=['placeholder'])
if not os.path.exists(atac_filename):
    atac_X.to_csv(atac_filename, sep="\t", compression="gzip")

# Identify all peaks that are within Xbp of annotated TSS
dist_filename = os.path.join(working_dir, tmp_path, "tssdist.tsv.gz")
os.system(f'python3 -m dictys chromatin tssdist --cut {distance//2} {proc_rna_filename} {atac_filename} {annot} {dist_filename}')


# Subset TF-Peak linking to preserve non-spurious TFs
tfb = pd.read_csv(tfb_data)
tfb.columns = ['loc', 'TF', 'score']  # What if the method output scores where sign matters, unlike distance?
tfb = tfb[tfb.TF.apply(lambda x: x in genes)]
tfb_name = os.path.join(working_dir, tmp_path, 'tfb_reformat.tsv.gz')
tfb[['TF', 'loc', 'score']].to_csv(tfb_name, sep='\t', index=False)


# Compute TF-Gene linking from TF-Peak and Peak-Gene linkings
tfg_name = os.path.join(working_dir, tmp_path, 'linking.tsv.gz')
mask_name = os.path.join(working_dir, tmp_path, 'binlinking.tsv.gz')
os.system(f'python3 -m dictys chromatin linking {tfb_name} {dist_filename} {tfg_name}')
os.system(f'python3 -m dictys chromatin binlinking {tfg_name} {mask_name} 20')


# Network inference
weight_name = os.path.join(working_dir, tmp_path, 'net_weight.tsv.gz')
meanvar_name = os.path.join(working_dir, tmp_path, 'net_meanvar.tsv.gz')
cov_name = os.path.join(working_dir, tmp_path, 'net_covfactor.tsv.gz')
loss_name = os.path.join(working_dir, tmp_path, 'net_loss.tsv.gz')
stats_name = os.path.join(working_dir, tmp_path, 'net_stats.tsv.gz')
normalized_weight_name = os.path.join(working_dir, tmp_path, 'net_nweight.tsv.gz')

os.system(f"python3 -m dictys network reconstruct --device {device} --nth {threads} {proc_rna_filename} {mask_name} {weight_name} {meanvar_name} {cov_name} {loss_name} {stats_name}")
os.system(f"python3 -m dictys network normalize --nth {threads} {weight_name} {meanvar_name} {cov_name} {normalized_weight_name}")


# Output TF-Gene scores. Only keep pairs that deemed to have TF-Gene links
weights = pd.read_csv(normalized_weight_name, sep="\t", index_col=0)
mask = pd.read_csv(mask_name, sep="\t", index_col=0)
mask = mask[weights.columns].loc[weights.index]
df = list()
for i in np.arange(weights.shape[0]):
    tmpdf = [(weights.index[i], weights.columns[j], weights.iloc[i,j])
             for j in np.arange(weights.shape[1]) if mask.iloc[i,j]]
    df.extend(tmpdf)
df = np.array(df)
df = pd.DataFrame(df, columns=['source', 'target', 'score'])
df['pval'] = 1
df.to_csv(path_out, index=False)
