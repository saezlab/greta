import argparse, os, sys
import pandas as pd
import numpy as np
import mudata as md


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-d', '--path_data', required=True)
parser.add_argument('-f', '--path_frag', required=True)
parser.add_argument('-w', '--working_dir', required=True)
parser.add_argument('-p', '--path_out', required=True)
parser.add_argument('-k', '--peak_file', required=True)
parser.add_argument('-m', '--motif_file', required=True)
parser.add_argument('-g', '--genome', required=True)
parser.add_argument('-t', '--threads', required=True)

args = vars(parser.parse_args())
mudata = args['path_data']
frag = args['path_frag']
working_dir = args['working_dir']
path_out = args['path_out']
peak_file = args['peak_file']
motif_file = args['motif_file']
genome_file = args['genome']
threads = int(args['threads'])


# Read in p2g and keep only peaks that are wide enough for footprinting
all_atac_peak = np.unique(pd.read_csv(peak_file)['cre'])
all_atac_peak = pd.DataFrame([n.split(':') for n in all_atac_peak])
all_atac_peak.columns = ['chr', 'srt', 'end']
all_atac_peak['srt'] = all_atac_peak['srt'].astype(int)
all_atac_peak['end'] = all_atac_peak['end'].astype(int)
all_atac_peak = all_atac_peak[(all_atac_peak.end - all_atac_peak.srt) >= 100]
all_atac_peak = all_atac_peak.sort_values(by=['chr', 'srt'])
valid_peak_fname = os.path.join(working_dir, "valid_peak.bed")
all_atac_peak.to_csv(valid_peak_fname, sep='\t', header=False, index=False)


# Write the RNA matrix to working directory
data = md.read(mudata)
rna_X = pd.DataFrame(np.array(data['rna'].layers['counts'].todense()).T, columns=data['rna'].obs.index, index=data['rna'].var.index)
rna_filename = os.path.join(working_dir, "expression.tsv.gz")
rna_X.to_csv(rna_filename, sep="\t", compression="gzip") 


# For each cluster of cells, compute TF to peak linking
clus = sorted(data.obs['celltype'].unique())
for count, c in enumerate(clus):
    cell_name = data[data.obs['celltype']==c].obs.index
    atac_fname = os.path.join(working_dir, f"name_clus_{count}.txt")
    with open(atac_fname, "w") as f:
        for n in cell_name:
            f.write(f"{n}\n")

    
    # Generate pseudo-bam file from fragments belong to the cell cluster
    bam_name = os.path.join(working_dir, f"clus_{count}.bam")
    bai_name = os.path.join(working_dir, f"clus_{count}.bai")
    command = (f'zcat {frag} | python {os.path.dirname(os.path.abspath(__file__))}/frag_to_bam.py --fname {atac_fname} | ' + 
               f'samtools view -b | samtools sort -o {bam_name} && samtools index {bam_name} {bai_name}')
    os.system(command)

    
    # Perform footprinting analysis and motif enrichment analysis
    footprint_name = os.path.join(working_dir, f"clus_{count}_footprints.tsv.gz")
    homer_name = os.path.join(working_dir, f"clus_{count}_homer.tsv.gz")
    wellington_name = os.path.join(working_dir, f"clus_{count}_wellington.tsv.gz")
    motif_name = os.path.join(working_dir, f"clus_{count}_motifs.tsv.gz")
    binding_name = os.path.join(working_dir, f"clus_{count}_binding.tsv.gz")
    os.system(f'python3 -m dictys chromatin wellington {bam_name} {bai_name} {valid_peak_fname} {footprint_name} --nth {threads}')
    os.system(f'python3 -m dictys chromatin homer {footprint_name} {motif_file} {genome_file} {rna_filename} {motif_name} {wellington_name} {homer_name}')
    os.system(f'python3 -m dictys chromatin binding {wellington_name} {homer_name} {binding_name}')


# Convert TF-TFBS files to bed format
for count, c in enumerate(clus):
    binding_name = os.path.join(working_dir, f"clus_{count}_binding.tsv.gz")
    df = pd.read_csv(binding_name, sep='\t')
    df['chr'] = df['loc'].apply(lambda x: x.split(':')[0])
    df['srt'] = df['loc'].apply(lambda x: int(x.split(':')[1]))
    df['end'] = df['loc'].apply(lambda x: int(x.split(':')[2]))
    footprint_tf_name = os.path.join(working_dir, f"clus_{count}_footprint_to_TF.bed")
    df[['chr', 'srt', 'end', 'TF', 'score']].to_csv(footprint_tf_name, sep='\t', header=False, index=False)
    
    valid_peak_to_tf_name = os.path.join(working_dir, f"valid_peak_to_clus_{count}_TF.bed")
    os.system(f'bedtools intersect -a {valid_peak_fname} -b {footprint_tf_name} -wa -wb > {valid_peak_to_tf_name}')


# Convert TF-TFBS to TF-Peak linking by mapping TFBS to Peak and averaging all TF activity within the Peak
overall_peak_to_TF = pd.concat([
    pd.read_csv(os.path.join(working_dir, f"valid_peak_to_clus_{count}_TF.bed"), sep='\t', header=None) 
    for count, c in enumerate(clus)
], axis=0)
overall_peak_to_TF.columns = ['peakchr', 'peaksrt', 'peakend', 'footprintchr', 'footprintsrt', 'footprintend', 'TF', 'score']
overall_peak_to_TF['peak'] = overall_peak_to_TF.peakchr + ':' + overall_peak_to_TF.peaksrt.astype(str) + ':' + overall_peak_to_TF.peakend.astype(str)
overall_peak_to_TF = overall_peak_to_TF.groupby(['peak', 'TF'], as_index=False)['score'].mean()
overall_peak_to_TF.columns = ['cre', 'tf', 'score']
overall_peak_to_TF.to_csv(path_out, index=False)
