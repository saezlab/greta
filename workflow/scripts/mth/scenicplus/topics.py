from pycisTopic.pseudobulk_peak_calling import export_pseudobulk, peak_calling
from pycisTopic.cistopic_class import create_cistopic_object
from pycisTopic.iterative_peak_calling import get_consensus_peaks
from pycisTopic.qc import get_barcodes_passing_qc_for_sample
import pycisTopic.fragments
import pandas as pd
import polars as pl
import mudata as mu
import tempfile
import scipy
import os
import yaml
import pickle
import argparse


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-m','--path_mdata', required=True)
parser.add_argument('-s','--lst_smpls', required=True, nargs='+')
parser.add_argument('-b','--path_bl', required=True)
parser.add_argument('-t','--path_tss', required=True)
parser.add_argument('-o','--organism', required=True)
parser.add_argument('-n','--n_topics', required=True)
parser.add_argument('-c','--n_cores', required=True)
parser.add_argument('-d','--out_dir', required=True)
parser.add_argument('-p','--path_pre', required=True)
args = vars(parser.parse_args())

path_mdata = args['path_mdata']
lst_smpls = args['lst_smpls']
path_bl = args['path_bl']
path_tss = args['path_tss']
organism = args['organism']
n_topics = int(args['n_topics'])
n_cores = int(args['n_cores'])
out_dir = args['out_dir']
path_pre = args['path_pre']

# Generate frags dict
frags_dict = {os.path.basename(path).split('.')[0]: path for path in lst_smpls}

# Read metadata
mdata = mu.read(path_mdata)
obs = mdata.obs.copy()
obs['celltype'] = obs['celltype'].str.replace(r'\s+', '_', regex=True)
# Write gex
os.makedirs(out_dir, exist_ok=True)
path_rna = os.path.join(out_dir, 'rna.h5ad')
mdata.mod['rna'].write(path_rna)
del mdata

# Get chromsizes
chromsizes = pd.read_table(
    f"http://hgdownload.cse.ucsc.edu/goldenPath/{organism}/bigZips/{organism}.chrom.sizes",
    header = None,
    names = ["Chromosome", "End"]
)
chromsizes.insert(1, "Start", 0)

# Compute psbulk per celltype
os.makedirs(os.path.join(out_dir, "consensus_peak_calling"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"), exist_ok=True)
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"), exist_ok=True)
bw_paths, bed_paths = export_pseudobulk(
    input_data = obs,
    variable = "celltype",
    sample_id_col = "batch",
    chromsizes = chromsizes,
    bed_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bed_files"),
    bigwig_path = os.path.join(out_dir, "consensus_peak_calling/pseudobulk_bw_files"),
    path_to_fragments = frags_dict,
    n_cpu = n_cores,
    normalize_bigwig = True,
    temp_dir = tempfile.gettempdir(),
)

# Call peaks
os.makedirs(os.path.join(out_dir, "consensus_peak_calling/MACS"), exist_ok = True)
organism_to_genome_size = {
    'hg38': 'hs',
    'mm10': 'mm',
    'dm6': 'dm',
}
narrow_peak_dict = peak_calling(
    macs_path = "macs2",
    bed_paths = bed_paths,
    outdir = os.path.join(os.path.join(out_dir, "consensus_peak_calling/MACS")),
    genome_size = organism_to_genome_size[organism],
    n_cpu = n_cores,
    input_format = 'BEDPE',
    shift = 73,
    ext_size = 146,
    keep_dup = 'all',
    q_value = 0.05,
    _temp_dir = tempfile.gettempdir()
)

# Get consensus peaks
consensus_peaks = get_consensus_peaks(
    narrow_peaks_dict = narrow_peak_dict,
    peak_half_width = 250,
    chromsizes = chromsizes,
    path_to_blacklist = path_bl
)

# Save to bed
path_cpeaks = os.path.join(out_dir, "consensus_peak_calling/consensus_regions.bed")
consensus_peaks.to_bed(
    path = path_cpeaks,
    keep = True,
    compression = 'infer',
    chain = False
)

# Run QC
path_qc = os.path.join(out_dir, 'qc')
os.makedirs(path_qc, exist_ok = True)
command = ''
for sample, path_frag in frags_dict.items():
    command += f"""printf "%s" "pycistopic qc run --fragments {path_frag} --regions {path_cpeaks} --tss {path_tss} --output {path_qc}/{sample}" """
command += f" | xargs -0 -n 1 -P {min([n_cores, len(frags_dict.keys())])}" + ' -I {} sh -c "{}"'
print(command, flush=True)
os.system(command)

# Generate binary mats
cistopic_obj_list = []
for sample in frags_dict:
    binary_matrix, cbs, region_ids = pycisTopic.fragments.create_fragment_matrix_from_fragments(
       frags_dict[sample],
       path_cpeaks,
       os.path.join(path_qc, f"{sample}.cbs_for_otsu_thresholds.tsv"),
    )
    binary_matrix.data.fill(1)
    cistopic_obj = create_cistopic_object(
        fragment_matrix=binary_matrix,
        cell_names=cbs,
        region_names=region_ids,
        path_to_blacklist=path_bl,
        project=sample,
        tag_cells=False,
        split_pattern="_",
    )
    cistopic_obj_list.append(cistopic_obj)
if len(cistopic_obj_list) > 1:
    cistopic_obj = cistopic_obj_list[0].merge(cistopic_obj_list[1:])
else:
    cistopic_obj = cistopic_obj_list[0]

# Write raw data
path_bmat = os.path.join(out_dir, "bool_mat.mtx")
path_cres = os.path.join(out_dir, "cres.txt")
path_brcs = os.path.join(out_dir, "brcs.txt")
scipy.io.mmwrite(path_bmat, scipy.sparse.csr_matrix(cistopic_obj.binary_matrix != 0.))
with open(path_cres, 'w') as f:
    for name in cistopic_obj.region_names:
        f.write(f"{name}\n")

with open(path_brcs, 'w') as f:
    for name in cistopic_obj.cell_names:
        f.write(f"{name}\n")

# Create corpus
command = f"pycistopic topic_modeling mallet create_corpus -i {path_bmat} -o {path_bmat}.mallet"
print(command, flush=True)
os.system(command)

# Compute topics
path_topics = os.path.join(out_dir, "topics")
os.makedirs(path_topics, exist_ok=True)
command = f"pycistopic topic_modeling mallet run \
-i {path_bmat}.mallet \
-o {path_topics}/topic_ \
-t {n_topics} \
-p {n_cores}"
print(command, flush=True)
os.system(command)

# Get topic model
command = f"pycistopic topic_modeling mallet stats \
-i {path_bmat} \
-c {path_brcs} \
-r {path_cres} \
-o {path_topics}/topic_ \
-t {n_topics}"
print(command, flush=True)
os.system(command)

# Load model
path_model = os.path.join(path_topics, f'topic_.{n_topics}_topics.model.pkl')
with open(path_model, 'rb') as pickle_file:
    model = pickle.load(pickle_file)
cistopic_obj.add_LDA_model(model)

# Write object
path_obj = os.path.join(out_dir, "cistopic_obj.pkl")
with open(path_obj, 'wb') as f:
    pickle.dump(cistopic_obj, f)

# Run scenicplus preprocessing
command = f"scenicplus prepare_data prepare_GEX_ACC \
--cisTopic_obj_fname {path_obj} \
--GEX_anndata_fname {path_rna} \
--out_file {path_pre} \
--do_not_use_raw_for_GEX_anndata"
print(command, flush=True)
os.system(command)

