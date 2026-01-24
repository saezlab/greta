import matplotlib.pyplot as plt  # Import else compiling error
import pandas as pd
import numpy as np
import scanpy as sc
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
import mudata as md
import os
import re
import argparse
import subprocess
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
import shutil
import glob
import itertools


# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_mudata', required=True)
parser.add_argument('-b','--path_output', required=True)
parser.add_argument('-c','--path_outfile', required=False, default=None)
parser.add_argument('-d','--path_chainFiles', required=True)
parser.add_argument('-e','--motif_dir', required=True)
parser.add_argument('-f','--promoter_dir', required=True)
parser.add_argument('-g','--window_size', required=True)
parser.add_argument('-i','--threads', required=True)
parser.add_argument('-j','--path_chainFiles_rev', required=True)


args = vars(parser.parse_args())

path_mudata = args['path_mudata']
path_output = args['path_output']
path_outfile = args['path_outfile']
path_chainFiles = args['path_chainFiles']
motif_dir = args['motif_dir']
motif_file = os.path.join(motif_dir, "human_all_motifs_sorted_clean.txt")
promoter_dir = args['promoter_dir']
promoter_scmtni = os.path.join(promoter_dir, "Homo_sapiens.GRCh37.74.TSS.5000.bed")
promoter_file = os.path.join(promoter_dir, "Homo_sapiens.GRCh37.74.TSS.customWindow.bed")
window_size = int(args['window_size'])
threads = int(args['threads'])
path_chainFiles_rev = args['path_chainFiles_rev']


# liftover files
liftover_path = "liftOver"  # path to the liftOver binary in sif image
input_bed = os.path.join(path_output, "peaks_hg38.bed")
output_bed = os.path.join(path_output, "peaks_hg19.bed")
unmapped_file = os.path.join(path_output, "unmapped.bed")

# narrow peak Paths
narrowPeak_paths = "peaks_by_cluster"

#-------------------------------------------------------------------------------------

# Load data and make sure celltypes don't have spaces
mdata = md.read_h5mu(path_mudata)
mdata.obs['celltype'] = mdata.obs['celltype'].str.replace(" ", "_")


#--------------------------------------------------------------------------------------
# Compute cell lineage tree

lineage_tree_file = os.path.join(path_output, 'lineage_tree.txt')
if not os.path.exists(lineage_tree_file):

    print("Computing cell lineage tree...")
    ## Get clusters
    clusters = mdata.obs['celltype']
    cluster_names = sorted(clusters.unique())

    rna = mdata.mod['rna']
    expr_df = pd.DataFrame(
        rna.X.toarray() if not isinstance(rna.X, np.ndarray) else rna.X,
        index=rna.obs_names,
        columns=rna.var_names
    )

    ## Save all genes to file
    # Extract gene names and save to a text file
    pd.Series(expr_df.columns).to_csv(os.path.join(path_output, "allGenes.txt"), index=False, header=False)

    ## Extract expression data and get cluster mean
    expr_df['cluster'] = clusters.loc[expr_df.index].values
    cluster_means = expr_df.groupby('cluster').mean()

    ## Compute pairwise distances between cluster means
    dist_matrix = squareform(pdist(cluster_means.values, metric="euclidean"))

    ## Compute minimum spanning tree
    mst_sparse = minimum_spanning_tree(dist_matrix)
    mst_arr = mst_sparse.toarray()

    ## Create graph from MST
    G = nx.Graph()
    for i, a in enumerate(cluster_names):
        for j, b in enumerate(cluster_names):
            w = mst_arr[i, j]
            if w != 0 and not np.isnan(w):
                G.add_edge(a, b, weight=float(w))

    ## Determine root cluster: Find node with the highest degree in MST
    degrees = dict(G.degree())
    root_cluster = max(degrees, key=degrees.get)
    print(f"Automatically selected root cluster: {root_cluster}")

    ## Generate lineage tree and celltype order
    ## DFS traversal from root cluster
    visited = set()
    edges = []
    def dfs(node):
        visited.add(node)
        for neighbor in sorted(G.neighbors(node)):
            if neighbor not in visited:
                dist = 0.2  # fixed time
                edges.append((node, neighbor, dist, dist))
                dfs(neighbor)

    dfs(root_cluster)

    edges_df = pd.DataFrame(edges, columns=["parent", "child", "time1", "time2"])
    edges_df = edges_df[["child", "parent", "time1", "time2"]]

    order = list(nx.dfs_tree(G, root_cluster).nodes())
    order_df = pd.DataFrame(order, columns=["celltype"])

    ## Save lineage tree and order
    edges_df.to_csv(os.path.join(path_output, "lineage_tree.txt"), sep="\t", index=False, header=False)
    order_df.to_csv(os.path.join(path_output, "celltype_order.txt"), sep="\t", index=False, header=False)

    print(f"Saved cell lineage tree to {os.path.join(path_output, 'lineage_tree.txt')}")
    print(f"Saved celltype order to {os.path.join(path_output, 'celltype_order.txt')}")

else:
    print(f"Skipping computing lineage tree, file already exists.")

#-------------------------------------------------------------------------------------
# LiftOver Peak Matrix to hg19 

if not os.path.exists(output_bed):

    print("Running liftOver to hg19...")
    ## Extract peak matrix
    peaks = pd.Series(mdata.mod['atac'].var_names)

    ## Convert peaks to a DataFrame in BED format
    bed_df = peaks.str.split('[:-]', expand=True)
    bed_df.columns = ["chrom", "start", "end"]
    bed_df["start"] = bed_df["start"].astype(int)
    bed_df["end"] = bed_df["end"].astype(int)

    bed_df.to_csv(input_bed, sep="\t", header=False, index=False)

    ## Run liftOver from Python
    subprocess.run([
        liftover_path,
        input_bed,
        path_chainFiles,
        output_bed,
        unmapped_file
    ], check=True)

else:
    print(f"Skipping liftOver, {output_bed} already exists.")


#---------------------------------------------------------------------
# Annotate mudata

annotated_file = os.path.join(path_output, "annotated_hg19.h5mu")
if not os.path.exists(annotated_file):

    ## Load lifted peaks
    lifted_df = pd.read_csv(output_bed, sep="\t", header=None, names=["chrom", "start", "end"])
    lifted_peaks = lifted_df["chrom"] + "-" + lifted_df["start"].astype(str) + "-" + lifted_df["end"].astype(str)

    ## Load unmapped peaks for filtering
    unmapped = pd.read_csv(unmapped_file, sep="\t", header=None, names=["chrom", "start", "end"])

    ## Filter out unmapped peaks
    ## Drop rows where 'chrom' starts with '#' or nan
    unmapped_clean = unmapped[~unmapped["chrom"].astype(str).str.startswith("#")]
    unmapped_clean = unmapped_clean.dropna(subset=["start", "end"])
    unmapped_clean['start'] = unmapped['start'].astype(str)
    unmapped_clean['end'] = unmapped['end'].astype(str)
    unmapped_clean = unmapped_clean["chrom"] + "-" + unmapped_clean["start"].astype(str) + "-" + unmapped_clean["end"].astype(str)

    # Split each string by '-', convert start and end to int, then rebuild string without .0
    def clean_peak_name(peak_str):
        chrom, start, end = peak_str.split('-')
        start = int(float(start))
        end = int(float(end))
        return f"{chrom}-{start}-{end}"

    cleaned = unmapped_clean.apply(clean_peak_name)

    ## Do filtering
    print("Filtering peaks...")
    peaks_to_remove = set(cleaned)
    current_peaks = mdata.mod['atac'].var_names
    peaks_to_keep = [peak for peak in current_peaks if peak not in peaks_to_remove]
    mdata.mod['atac'] = mdata.mod['atac'][:, peaks_to_keep]

    print(f"Peaks before: {len(current_peaks)}")
    print(f"Peaks after: {mdata.mod['atac'].n_vars}")

    ## Replace var_names with lifted peaks
    mdata.mod['atac'].var_names = lifted_peaks.values

    # Re-create MuData object from filtered modalities and save
    mdata_filtered = md.MuData(mdata.mod)
    mdata_filtered.obs = mdata.obs.copy()
    mdata_filtered.write_h5mu(os.path.join(path_output, "annotated_hg19.h5mu"))
    print("Saved filtered MuData object to annotated_hg19.h5mu")

else:
    print(f"Skipping annotated_hg19.h5mu, file already exists.")

#-------------------------------------------------------------------------------------
# Split peak matrix by cell type and save as narrowPeak files

for celltype in mdata.obs['celltype'].unique():
    out_path = os.path.join(path_output, f"{celltype}.narrowPeak")
    if os.path.exists(out_path):
        print(f"Skipping {out_path}, already exists.")
        continue

    print("Saving peaks by cell type as narrowPeak files...")

    # Get ATAC modality
    mdata = md.read_h5mu(os.path.join(path_output, "annotated_hg19.h5mu"))
    atac = mdata.mod['atac']

    ## narrowPeak columns: chrom, start, end, name, score, strand, signalValue, pValue, qValue, peak
    template_cols = ["chrom", "start", "end", "name", "score", "strand", "signalValue", "pValue", "qValue", "peak"]

    # Iterate over clusters
    for celltype, cells in mdata.obs.groupby("celltype").groups.items():
        # Subset ATAC matrix to cluster cells
        submatrix = atac[cells, :].X

        # Convert sparse to dense if needed
        if not isinstance(submatrix, np.ndarray):
            submatrix = submatrix.toarray()

        # Find peaks present in this cluster
        present_peaks = (submatrix > 0).sum(axis=0) > 0

        # Mean accessibility across cells in cluster
        mean_signal = submatrix.mean(axis=0)

        # Build DataFrame for narrowPeak
        df = lifted_df.loc[present_peaks].copy()
        df["name"] = [f"peak{i+1}_{celltype}" for i in range(len(df))]
        df["score"] = 0
        df["strand"] = "."
        df["signalValue"] = mean_signal[present_peaks]
        df["pValue"] = -1
        df["qValue"] = -1
        df["peak"] = -1

        # Save in narrowPeak format
        df = df[template_cols]
        out_path = os.path.join(path_output, f"{celltype}.narrowPeak")
        df.to_csv(out_path, sep="\t", header=False, index=False)


    print(f"Saved peaks for {len(mdata.obs['celltype'].unique())} clusters to {path_output}/")



#-------------------------------------------------------------------------------------
# Split expression matrix by cell type and save as txt files

for celltype in mdata.obs['celltype'].unique():
    out_path = os.path.join(path_output, f"{celltype}.table")
    if os.path.exists(out_path):
        print(f"Skipping {out_path}, already exists.")
        continue

    print("Saving expression matrices by cell type...")

    # Get expression matrix and celltype annotations
    mdata = md.read_h5mu(os.path.join(path_output, "annotated_hg19.h5mu"))
    adata = mdata.mod["rna"]
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    expr = adata.to_df()   # cells x genes
    celltypes = mdata.obs["celltype"]


    # Loop over celltypes and save one table per cluster
    for cluster in celltypes.unique():
        cells_in_cluster = celltypes[celltypes == cluster].index
        expr_cluster = expr.loc[cells_in_cluster]  # subset matrix

        # transpose → rows = genes, cols = cells
        expr_cluster = expr_cluster.T

        # rename genes to include cluster suffix
        expr_cluster.index = [f"{g}_{cluster}" for g in expr_cluster.index]

        # save to file
        outfile = os.path.join(path_output, f"{cluster}.table")
        expr_cluster.to_csv(outfile, sep="\t")

        print(f"Saved: {outfile}")

#--------------------------------------------------------------------------------------
# Create filelist
print("Creating filelist...")

# directory where you stored the per-cluster tables
tables_dir = Path(path_output)
config_file = Path(os.path.join(path_output, "filelist.txt"))

# find all cluster tables (assuming naming: clusterX.table)
cluster_tables = sorted(tables_dir.glob("*.table"))

with open(config_file, "w") as f:
    for table in cluster_tables:
        cluster_name = table.stem  # "cluster3" from "cluster3.table"
        f.write(f"{cluster_name}\t{table}\n")

print(f"Config file written to {config_file}")


#--------------------------------------------------------------------------------------
# Download TFs from Lambert
print("Downloading TF list from Lambert et al. 2018...")
path = os.path.join(path_output, "regulators.txt")

subprocess.run(["wget", "-O", path, "https://humantfs.ccbr.utoronto.ca/download/v_1.01/TF_names_v_1.01.txt"])


#-------------------------------------------------------------------------------------
# Prepare input files for prior network generation
print("Preparing input files for prior network generation...")

# define paths
indir = path_output
filelist = os.path.join(path_output, "filelist.txt")
regfile = os.path.join(path_output, "regulators.txt")
outdir = Path(os.path.join(path_output, "Results/"))

# make sure output dir exists
outdir.mkdir(parents=True, exist_ok=True)

# build the command
cmd = [
    "python", "/opt/scMTNI/Scripts/PreparescMTNIinputfiles.py",
    "--filelist", str(filelist),
    "--regfile", str(regfile),
    "--indir", str(indir),
    "--outdir", str(outdir),
    "--splitgene", "25",
    "--motifs", "1"
]

# run it
subprocess.run(cmd, check=True)

print(f"PreparescMTNIinputfiles.py finished, results in {outdir}")


#------------------------------------------------
#Extend promoter regions
#------------------------------------------------

# Correct window: subtract 10 kb (because ±5 kb already included)
corrected_window = max(window_size - 10000, 0)
half_extension = corrected_window // 2

cmd = f"""
awk -v OFS='\\t' -v ext={half_extension} '
BEGIN {{ FS=OFS="\t" }}
{{ 
    start = $2 - ext; 
    end = $3 + ext;
    if (start < 0) start = 0;
    print $1, start, end, $4, $5, $6, $7;
}}' "{promoter_scmtni}" > "{promoter_file}"
"""

print(f"Running:\n{cmd}")
subprocess.run(cmd, shell=True, check=True)
print(f"Extended promoter regions saved to {promoter_file}")


#----------------------------------------------
# Define preprocessing functions
#----------------------------------------------
# 1. Sort BED files
#------------------------
def sort_bed(in_file, out_file):
    in_file = Path(in_file)
    out_file = Path(out_file)
    if out_file.exists():
        print(f"Skipping sorting: {out_file} already exists")
        return
    subprocess.run(
        ["sort", "-k1,1", "-k2,2n", str(in_file)],
        stdout=open(out_file, "w"),
        check=True
    )
    print(f"Sorted {in_file} -> {out_file}")


#------------------------
# 2. Filter to standard chromosomes
#------------------------
def filter_standard_chroms(in_file):
    """Filter BED file to keep only standard chromosomes (chr1-22, X, Y, M)."""
    in_file = Path(in_file)
    tmp_file = in_file.with_suffix(".tmp")
    subprocess.run(
        f"grep -E '^chr([0-9]+|X|Y|M)\t' {in_file} > {tmp_file} && mv {tmp_file} {in_file}",
        shell=True,
        check=True
    )
    print(f"Filtered {in_file} to standard chromosomes")


#------------------------
# 3. Run bedtools intersects
#------------------------
def run_bedtools_intersects(cell_types, outdir, motif_file, promoter_file, bedtools_path="bedtools"):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for cluster in cell_types:
        peaks_file = Path(outdir) / f"{cluster}.sorted.narrowPeak"
        motif_out = Path(outdir) / f"motifs_in_{cluster}"
        tss_out = Path(outdir) / f"TSS_in_{cluster}"

        if not motif_out.exists():
            subprocess.run([
                bedtools_path, "intersect", "-a", str(peaks_file),
                "-b", motif_file, "-wb", "-sorted"
            ], stdout=open(motif_out, "w"), check=True)
            print(f"Generated {motif_out}")
        else:
            print(f"Skipping {motif_out}, already exists")

        if not tss_out.exists():
            subprocess.run([
                bedtools_path, "intersect", "-a", str(peaks_file),
                "-b", promoter_file, "-wb", "-sorted"
            ], stdout=open(tss_out, "w"), check=True)
            print(f"Generated {tss_out}")
        else:
            print(f"Skipping {tss_out}, already exists")


#------------------------
# 3. Map motifs to genes
#------------------------

def run_map_motifs_to_genes(cell_types, outdir, motif2tf_file):
    Path(outdir).mkdir(parents=True, exist_ok=True)
    for cluster in cell_types:
        outfile = Path(outdir) / f"{cluster}_network.txt"
        if outfile.exists():
            print(f"Skipping {outfile}, already exists")
            continue
        subprocess.run([
            "python", "/opt/scMTNI/Scripts/genPriorNetwork/mapMot2Gene.py",
            "--mot2tf", motif2tf_file,
            "--mot2peak", str(Path(outdir) / f"motifs_in_{cluster}"),
            "--peak2gene", str(Path(outdir) / f"TSS_in_{cluster}"),
            "--outfile", str(outfile)
        ], check=True)
        print(f"Mapped motifs -> genes for {cluster}")


#------------------------
# 4. Filter prior network
#------------------------
def run_filter_prior_network(cell_types, outdir, outdir_prior):
    Path(outdir_prior).mkdir(parents=True, exist_ok=True)
    for sample in cell_types:
        regfile = Path(outdir) / f"{sample}_allregulators.txt"
        genefile = Path(outdir) / f"{sample}_allGenes.txt"
        netfile0 = Path(outdir) / f"{sample}_network.txt"
        netfile = Path(outdir) / f"{sample}.txt"
        outfile = Path(outdir_prior) / f"{sample}_network.txt"

        if outfile.exists():
            print(f"Skipping {outfile}, already exists")
            continue

        # Rename nodes
        with open(netfile0) as fin, open(netfile, "w") as fout:
            for line in fin:
                parts = line.strip().split("\t")
                parts[0] = f"{parts[0]}_{sample}"
                parts[1] = f"{parts[1]}_{sample}"
                fout.write("\t".join(parts) + "\n")

        # Run filtering script
        subprocess.run([
            "python", "/opt/scMTNI/Scripts/genPriorNetwork/filterpriornetwork.py",
            "--regfile", str(regfile),
            "--genefile", str(genefile),
            "--netfile", str(netfile),
            "--outfile", str(outfile)
        ], check=True)
        print(f"Filtered prior network for {sample}")


#------------------------
# 5. Filter top edges
#------------------------
def run_filter_top_edges(outdir, outdir_prior):
    # Here we assume the Rscript produces files in outdir4
    subprocess.run([
        "Rscript", "--vanilla", "/opt/scMTNI/Scripts/genPriorNetwork/filtertop20Pedges.R",
        str(outdir), 
        str(outdir_prior)
    ], check=True)
    print(f"Filtered top edges")


#------------------------
# 6. Percentile ranking
#------------------------
def run_percentile_ranking(cell_types, outdir_top02, outdir_ranked):
    Path(outdir_ranked).mkdir(parents=True, exist_ok=True)
    for sample in cell_types:
        networkfile = Path(outdir_top02) / f"{sample}_network.txt"
        outfile = Path(outdir_ranked) / f"{sample}_network.txt"
        if outfile.exists():
            print(f"Skipping {outfile}, already exists")
            continue
        subprocess.run([
            "/opt/scMTNI/Scripts/genPriorNetwork/rankEdges",
            str(networkfile),
            str(outfile),
            "incr"
        ], check=True)
        print(f"Percentile ranked network for {sample}")
        
        
#========================================================
# Run preprocessing
#========================================================

# Define celltypes
celltypes = list(mdata.obs["celltype"].unique())


# Sorted peak files 
print("Sorting peak files...")


for ct in celltypes:
    inpath = Path(path_output) / f"{ct}.narrowPeak"
    sorted_peaks = inpath.with_name(inpath.stem + ".sorted" + inpath.suffix)

    sort_bed(inpath, sorted_peaks)
    filter_standard_chroms(sorted_peaks)


print("Running bedtools intersects...")
run_bedtools_intersects(cell_types=celltypes, 
                        outdir = path_output,
                        motif_file = motif_file,
                        promoter_file = promoter_file)


print("Mapping motifs to genes...")
run_map_motifs_to_genes(cell_types=celltypes,
                        outdir = path_output,
                        motif2tf_file= "/opt/scMTNI/ExampleData/motifs/cisbp_motif2tf.txt")

print("Filtering prior network...")
run_filter_prior_network(cell_types=celltypes, 
                        outdir= path_output, 
                        outdir_prior= os.path.join(path_output, "prior_networks/"))

print("Filtering top 20% edges...")
run_filter_top_edges(outdir = path_output, 
                    outdir_prior = os.path.join(path_output, "prior_networks/"))


print("Running percentile ranking...")
run_percentile_ranking(cell_types=celltypes, 
                    outdir_top02 = os.path.join(path_output, "prior_networks_top0.2/"), 
                    outdir_ranked = os.path.join(path_output, "prior_networks_ranked/"))

#---------------------------------------------------------------------------------
# Move filtered and ranked prior networks to main output dir

def move_ranked_networks(cell_types, outdir_ranked, outdir):
    outdir_ranked = Path(outdir_ranked)
    outdir = Path(outdir)

    for sample in cell_types:
        ranked_file = outdir_ranked / f"{sample}_network.txt"
        dest_file = outdir / f"{sample}_network.txt"

        if not ranked_file.exists():
            print(f"Skipping {sample}, ranked file does not exist: {ranked_file}")
            continue

        if dest_file.exists():
            dest_file.unlink()  # remove existing file

        shutil.move(str(ranked_file), str(dest_file))
        print(f"Moved {ranked_file} -> {dest_file} (overwriting if existed)")


move_ranked_networks(cell_types=celltypes, 
                     outdir_ranked = os.path.join(path_output, "prior_networks_ranked/"),
                     outdir = path_output)


#---------------------------------------------------------------------------------
# -------- Helper to get AllGenes split files --------


os.environ["OMP_NUM_THREADS"] = "1"
os.environ["MKL_NUM_THREADS"] = "1"


def get_allgenes_files(ogids_dir):
    ogids_dir = Path(ogids_dir)
    pattern = re.compile(r"AllGenes(\d+)\.txt")
    files = []
    for f in ogids_dir.glob("AllGenes*.txt"):
        if pattern.match(f.name):
            files.append(f)
    # sort by index
    files.sort(key=lambda x: int(pattern.match(x.name).group(1)))
    return files

# -------- Function to run scMTNI for a single batch --------
def run_scmtni_for_batch(datadir, gene_file):
    data = Path(datadir)
    cmd = [
        "/opt/scMTNI/Code/scMTNI",
        "-f", str(data / "testdata_config.txt"),
        "-x5",
        "-l", str(data / "TFs_OGs.txt"),
        "-n", str(gene_file),
        "-d", str(data / "lineage_tree.txt"),
        "-m", str(data / "testdata_ogids.txt"),
        "-s", str(data / "celltype_order.txt"),
        "-p", "0.5",
        "-c", "yes",
        "-b", "-0.9",
        "-q", "2"
    ]
    print("Running:", " ".join(cmd))
    subprocess.run(cmd, check=True)
    return gene_file

# -------- Main parallel runner --------
batch_size = threads  # number of parallel scMTNI runs at a time

def parallel_scmtni(datadir, max_workers=batch_size):
    ogids_dir = Path(datadir) / "ogids"
    gene_files = get_allgenes_files(ogids_dir)

    
    print(f"Found {len(gene_files)} split files. Launching parallel execution...")


    for batch in [gene_files[i:i+batch_size] for i in range(0, len(gene_files), batch_size)]:
        with ThreadPoolExecutor(max_workers=batch_size) as executor:
            futures = {executor.submit(run_scmtni_for_batch, path_output, gf): gf for gf in batch}
            for future in as_completed(futures):
                gf = futures[future]
                try:
                    future.result()
                    print(f"Finished {gf}")
                except Exception as e:
                    print(f"Error in {gf}: {e}")


# Run parallel scMTNI
parallel_scmtni(
    datadir=path_output,
    max_workers=batch_size
)

#----------------------------------
# Remap TF-gene edge back to CRE
#----------------------------------
# Helper functions
#------------------------

def extract_cres_by_peakname(network_file, mot2tf_file, mot2peak_file, peak2gene_file, topn=None):
    # --- Load motif2TF ---
    mot2tf = {}
    with open(mot2tf_file) as f:
        for l in f:
            parts = l.strip().split('\t')
            if len(parts) < 2:
                continue
            mot2tf[parts[0]] = parts[1].split('::')

    # --- Load motif2peak (peak names only) ---
    mot2peak = {}
    with open(mot2peak_file) as f:
        for l in f:
            parts = l.strip().split('\t')
            if len(parts) < 15:
                continue
            peak_name = parts[3]
            motif = parts[-3]
            mot2peak.setdefault(motif, set()).add(peak_name)

    # --- Load peak2gene (peak name to genes, base names only) ---
    peak2gene = {}
    with open(peak2gene_file) as f:
        for l in f:
            parts = l.strip().split('\t')
            if len(parts) < 17:
                continue
            peak_name = parts[3]
            gene = parts[-1]
            # remove potential cluster suffix from gene
            gene_base = gene.split("_")[0]
            peak2gene.setdefault(peak_name, set()).add(gene_base)

    # --- Load TF–gene network ---
    net_df = pd.read_csv(network_file, sep="\t", header=None, names=["tf", "gene", "score"])
    rows = []

    # --- Map TF -> motif -> peak name -> gene ---
    for _, row in net_df.iterrows():
        tf_full, gene_full, score = row["tf"], row["gene"], row["score"]
        tf_name = tf_full.split("_")[0]      # removing celltype suffix
        gene_base = gene_full.split("_")[0]  # removing celltype suffix

        motifs = [mot for mot, tfs in mot2tf.items() if tf_name in tfs]
        if not motifs:
            continue

        for mot in motifs:
            if mot not in mot2peak:
                continue
            for peak_name in mot2peak[mot]:
                if peak_name not in peak2gene:
                    continue
                if gene_base not in peak2gene[peak_name]:
                    continue
                rows.append((tf_full, peak_name, gene_full, score))

    # --- Create DataFrame ---
    df = pd.DataFrame(rows, columns=["source", "CRE", "target", "score"])

    # --- Collapse duplicates by max score ---
    if not df.empty:
        df = df.groupby(["source", "CRE", "target"], as_index=False).score.max()

        # --- Optional: keep top-n CREs per TF–gene pair ---
        if topn is not None:
            df = df.sort_values("score", ascending=False)
            df = df.groupby(["source", "target"]).head(topn)

    return df



def replace_peak_names_with_coords(df, peak2gene_file, cre_net):

    # --- Build peak name -> coordinate mapping ---
    peak2coord = {}
    with open(peak2gene_file) as f:
        for line in f:
            if not line.strip():
                continue
            parts = line.strip().split('\t')
            if len(parts) < 17:
                continue
            chrom, start, end = parts[0], parts[1], parts[2]
            peak_name = parts[3]
            coord = f"{chrom}-{start}-{end}"
            peak2coord[peak_name] = coord

    # --- Replace peak names in dataframe ---
    df["CRE"] = df["CRE"].map(lambda x: peak2coord.get(x, x)) 

    df.to_csv(cre_net, index=None)
    return df



#------------------------
# Map all cell types
#------------------------
def remap_all_celltypes(results_dir, outdir, cell_types, mot2tf_file):
    results_dir = Path(results_dir)
    outdir = Path(outdir)

    for sample in cell_types:

        fold0_dir = results_dir / f"{sample}/fold0/"
        if not fold0_dir.exists():
            print(f"No fold0 directory found for {sample}, skipping")
            continue

        # assume existing network is named "network.txt"
        network_file = fold0_dir / "var_mb_pw_k5.txt"
        if not network_file.exists():
            print(f"No var_mb_pw_k5.txt found in {fold0_dir}, skipping")
            continue
        
        cre_net = fold0_dir / "cre_net.txt"

        # peak2gene and motif2peak files for this cell type
        peak2gene_file = outdir / f"TSS_in_{sample}"
        mot2peak_file  = outdir / f"motifs_in_{sample}"

        if not peak2gene_file.exists() or not mot2peak_file.exists():
            print(f"Missing peak2gene or motif2peak for {sample}, skipping")
            continue


        print(f"Processing cell type {sample}...")

        final_net = extract_cres_by_peakname(
        network_file= network_file,                        
        mot2tf_file=mot2tf_file,                            
        mot2peak_file=mot2peak_file,                        
        peak2gene_file=peak2gene_file                       
        )
        
        
        final_net_coords = replace_peak_names_with_coords(
        final_net, 
        peak2gene_file=peak2gene_file,                      
        cre_net=cre_net
        )    

        print(f"Done {sample}\n")


#------------------------
#  Run remapping
#------------------------

print(celltypes)
resultsdir=os.path.join(path_output, "Results/")
print(resultsdir)

remap_all_celltypes(results_dir=resultsdir,
                    outdir=path_output, 
                    cell_types=celltypes,
                    mot2tf_file="/opt/scMTNI/ExampleData/motifs/cisbp_motif2tf.txt")


#-------------------------------------------------------------------------------------
# Merge networks across cell types

def merge_networks(outdir):
    edges = {}
    files = glob.glob(outdir + "/Results/*/fold0/cre_net.txt")

    for file in files:
        df = pd.read_csv(file, sep=",")
        # strip cluster suffix to merge across cell typess
        df["source"] = df["source"].str.replace("_.*", "", regex=True)
        df["target"] = df["target"].str.replace("_.*", "", regex=True)
        df["score"] = pd.to_numeric(df["score"], errors="coerce")

        for _, row in df.iterrows():
            pair = (row["source"], row["CRE"], row["target"])
            edges.setdefault(pair, []).append(row["score"])

    consensus = []
    for (tf, cre, target), score in edges.items():
        arr = np.array(score)
        # maximum absolute value
        max_abs = np.max(np.abs(arr))
        # sign product (ignoring zeros)
        signs = np.sign(arr[arr != 0])
        if len(signs) == 0:
            sign = 0
        else:
            sign = np.prod(signs)
        consensus_val = max_abs * sign
        consensus.append((tf, cre, target, consensus_val))

    consensus_df = pd.DataFrame(consensus, columns=["source", "cre", "target", "score"])

    return consensus_df


#-------------------------------------------------------------------------------------
# LiftOver CRE coordinates from hg19 back to hg38

def liftover_cre_to_hg38(consensus_df, chain_file, output_dir):
    """
    Lift CRE coordinates from hg19 back to hg38.
    """
    # Extract unique CREs
    cres = consensus_df['cre'].unique()

    # Write to BED file with name column for tracking
    input_bed = os.path.join(output_dir, "cre_hg19.bed")
    output_bed = os.path.join(output_dir, "cre_hg38.bed")
    unmapped = os.path.join(output_dir, "cre_unmapped.bed")

    with open(input_bed, 'w') as f:
        for cre in cres:
            parts = cre.split('-')
            if len(parts) == 3:
                # BED format: chrom, start, end, name
                f.write(f"{parts[0]}\t{parts[1]}\t{parts[2]}\t{cre}\n")

    # Run liftOver
    subprocess.run([
        "liftOver",
        input_bed,
        chain_file,
        output_bed,
        unmapped
    ], check=True)

    # Build mapping using name column
    hg19_to_hg38 = {}
    with open(output_bed) as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 4:
                hg38_coord = f"{parts[0]}-{parts[1]}-{parts[2]}"
                hg19_coord = parts[3]  # original name
                hg19_to_hg38[hg19_coord] = hg38_coord

    # Update consensus dataframe
    n_before = len(consensus_df)
    consensus_df['cre'] = consensus_df['cre'].map(lambda x: hg19_to_hg38.get(x, None))

    # Remove rows where CRE couldn't be lifted
    consensus_df = consensus_df.dropna(subset=['cre'])
    n_after = len(consensus_df)

    if n_before != n_after:
        print(f"Removed {n_before - n_after} edges with unmappable CREs")

    print(f"Successfully lifted {len(hg19_to_hg38)} CREs from hg19 to hg38")

    return consensus_df


# Merge networks across cell types
print("Merging networks across cell types...")
consensus = merge_networks(outdir=path_output)

# Lift CRE coordinates back to hg38
print("Lifting CRE coordinates from hg19 to hg38...")
consensus = liftover_cre_to_hg38(consensus, path_chainFiles_rev, path_output)

# Save final output
consensus.to_csv(path_outfile, index=False)
print(f"Saved final network with hg38 coordinates to {path_outfile}")

# Safe cleanup
shutil.rmtree(path_output, ignore_errors=True)
print(f"Cleaned up temporary directory {path_output}")
