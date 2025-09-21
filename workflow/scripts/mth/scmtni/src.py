import matplotlib.pyplot as plt  # Import else compiling error
import pandas as pd
import numpy as np
import scanpy as sc
import networkx as nx
from scipy.spatial.distance import pdist, squareform
from scipy.sparse.csgraph import minimum_spanning_tree
import mudata as md
import os
import argparse
import subprocess
from pathlib import Path
import shutil



# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a','--path_mudata', required=True)
parser.add_argument('-b','--path_output', required=True)
parser.add_argument('-c','--path_outfile', required=False, default=None)
parser.add_argument('-d','--path_chainFiles', required=True)
parser.add_argument('-e','--motif_file', required=True)
parser.add_argument('-f','--promoter_file', required=True)

args = vars(parser.parse_args())

path_mudata = args['path_mudata']
path_output = args['path_output']
path_outfile = args['path_outfile']
path_chainFiles = args['path_chainFiles']
motif_file = args['motif_file']
promoter_file = args['promoter_file']


# liftover files
liftover_path = "liftOver"  # path to the liftOver binary in sif image
input_bed = os.path.join(path_output, "peaks_hg38.bed")
output_bed = os.path.join(path_output, "peaks_hg19.bed")
unmapped_file = os.path.join(path_output, "unmapped.bed")

# narrow peak Paths
narrowPeak_paths = "peaks_by_cluster"



n_jobs = 32

#-------------------------------------------------------------------------------------

# Load data and make sure celltypes don't have spaces
mdata = md.read_h5mu(path_mudata)
mdata.obs['celltype'] = mdata.obs['celltype'].str.replace(" ", "_")


#--------------------------------------------------------------------------------------
# Compute cell lineage tree
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
            dist = 0.2  # fixed time, or use: G[node][neighbor]["weight"]
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


#-------------------------------------------------------------------------------------
# LiftOver Peak Matrix to hg19 
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

#-------------------------------------------------------------------------------------
# Split peak matrix by cell type and save as narrowPeak files
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
    "--splitgene", "50",
    "--motifs", "1"
]

# run it
subprocess.run(cmd, check=True)

print(f"PreparescMTNIinputfiles.py finished, results in {outdir}")


#-------------------------------------------------------------------------------------
# Generate prior network
print("Generating prior network...")


#------------------------
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
# 2. Run bedtools intersects
#------------------------
def run_bedtools_intersects(cell_types, outdir, outdir2, motifs_file, promoter_file, bedtools_path="/opt/scMTNI/Scripts/genPriorNetwork/bedtools"):
    Path(outdir2).mkdir(parents=True, exist_ok=True)
    for cluster in cell_types:
        peaks_file = Path(outdir) / f"{cluster}.sorted.narrowPeak"
        motifs_out = Path(outdir2) / f"motifs_in_{cluster}"
        tss_out = Path(outdir2) / f"TSS_in_{cluster}"

        if not motifs_out.exists():
            subprocess.run([
                bedtools_path, "intersect", "-a", str(peaks_file),
                "-b", motifs_file, "-wb", "-sorted"
            ], stdout=open(motifs_out, "w"), check=True)
            print(f"Generated {motifs_out}")
        else:
            print(f"Skipping {motifs_out}, already exists")

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
def run_map_motifs_to_genes(cell_types, outdir2, outdir3, motif2tf_file):
    Path(outdir3).mkdir(parents=True, exist_ok=True)
    for cluster in cell_types:
        outfile = Path(outdir3) / f"{cluster}_network.txt"
        if outfile.exists():
            print(f"Skipping {outfile}, already exists")
            continue
        subprocess.run([
            "python", "/opt/scMTNI/Scripts/genPriorNetwork/mapMot2Gene.py",
            "--mot2tf", motif2tf_file,
            "--mot2peak", str(Path(outdir2) / f"motifs_in_{cluster}"),
            "--peak2gene", str(Path(outdir2) / f"TSS_in_{cluster}"),
            "--outfile", str(outfile)
        ], check=True)
        print(f"Mapped motifs -> genes for {cluster}")


#------------------------
# 4. Filter prior network
#------------------------
def run_filter_prior_network(cell_types, datadir, outdir3, outdir4):
    Path(outdir4).mkdir(parents=True, exist_ok=True)
    for sample in cell_types:
        regfile = Path(datadir) / f"{sample}_allregulators.txt"
        genefile = Path(datadir) / f"{sample}_allGenes.txt"
        netfile0 = Path(outdir3) / f"{sample}_network.txt"
        netfile = Path(outdir3) / f"{sample}.txt"
        outfile = Path(outdir4) / f"{sample}_network.txt"

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
def run_filter_top_edges(datadir, outdir4):
    outdir4 = Path(outdir4)
    outdir4.mkdir(parents=True, exist_ok=True)
    # Here we assume the Rscript produces files in outdir4
    subprocess.run([
        "Rscript", "--vanilla", "/opt/scMTNI/Scripts/genPriorNetwork/filtertop20Pedges_mod.R",
        datadir, str(outdir4)
    ], check=True)
    print(f"Filtered top edges -> {outdir4}")


#------------------------
# 6. Percentile ranking
#------------------------
def run_percentile_ranking(cell_types, outdir4, outdir5):
    Path(outdir5).mkdir(parents=True, exist_ok=True)
    for sample in cell_types:
        networkfile = Path(outdir4) / f"{sample}_network.txt"
        outfile = Path(outdir5) / f"{sample}_network.txt"
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




# Sorted peak files 
print("Sorting peak files...")

for ct in celltypes:
    inpath = Path(path_output) / f"{ct}.narrowPeak"
    sorted_peaks = inpath.with_name(inpath.stem + ".sorted" + inpath.suffix)

    sort_bed(inpath, sorted_peaks)


print("Running bedtools intersects...")
run_bedtools_intersects(cell_types=celltypes, 
                        outdir = path_output,
                        outdir2 = path_output,
                        motifs_file = motif_file,
                        promoter_file = promoter_file)


print("Mapping motifs to genes...")
run_map_motifs_to_genes(cell_types=celltypes,
                        outdir2 = path_output,
                        outdir3 = path_output,
                        motif2tf_file= "/opt/scMTNI/ExampleData/motifs/cisbp_motif2tf.txt")

print("Filtering prior network...")
run_filter_prior_network(cell_types=celltypes, 
                        datadir= path_output, 
                        outdir3= path_output, 
                        outdir4= os.path.join(path_output, "prior_networks/"))

print("Filtering top 20% edges...")
run_filter_top_edges(datadir = path_output, 
                    outdir4 = os.path.join(path_output, "prior_networks_top0.2/"))


print("Running percentile ranking...")
run_percentile_ranking(cell_types=celltypes, 
                       outdir4 = os.path.join(path_output, "prior_networks_top0.2/"), 
                       outdir5 = os.path.join(path_output, "prior_networks_ranked/"))


#---------------------------------------------------------------------------------
# Move filtered and ranked prior networks to main output dir

def move_ranked_networks(cell_types, outdir5, path_output):
    outdir5 = Path(outdir5)
    path_output = Path(path_output)
    path_output.mkdir(parents=True, exist_ok=True)

    for sample in cell_types:
        ranked_file = outdir5 / f"{sample}_network.txt"
        dest_file = path_output / f"{sample}_network.txt"

        if not ranked_file.exists():
            print(f"Skipping {sample}, ranked file does not exist: {ranked_file}")
            continue

        if dest_file.exists():
            dest_file.unlink()  # remove existing file

        shutil.move(str(ranked_file), str(dest_file))
        print(f"Moved {ranked_file} -> {dest_file} (overwriting if existed)")


move_ranked_networks(cell_types=celltypes, 
                     outdir5 = os.path.join(path_output, "prior_networks_ranked/"),
                     path_output = path_output)

#-------------------------------------------------------------------------------------
# Run scMTNI

print("Running scMTNI...")

def run_scmtni(datadir):
    """
    Run the scMTNI command inside the container
    
    Parameters
    ----------
    datadir : str or Path
    """
    data = Path(datadir)

    cmd = [
        "/opt/scMTNI/Code/scMTNI",
        "-f", str(data / "testdata_config.txt"),
        "-x50",
        "-l", str(data / "TFs_OGs.txt"),
        "-n", str(data / "AllGenes.txt"),
        "-d", str(data / "lineage_tree.txt"),
        "-m", str(data / "testdata_ogids.txt"),
        "-s", str(data / "celltype_order.txt"),
        "-p", "0.5",
        "-c", "yes",
        "-b", "-0.9",
        "-q", "2"
    ]

    print("Running command:\n", " ".join(cmd))
    subprocess.run(cmd, check=True)


run_scmtni(
    datadir = path_output)


#-------------------------------------------------------------------------------------
# Merge networks across cell types


def merge_networks(results_dir, output_file):
    edges = {}
    files = glob.glob(results_dir + "Results/*/fold0/var_mb_pw_k50.txt")

    for file in files:
        df = pd.read_csv(file, sep="\t", header=None, names=["source", "target", "weight"])

        # strip cluster suffix to merge across cell types
        df["source"] = df["source"].str.replace(r"_cluster\d+", "", regex=True)
        df["target"] = df["target"].str.replace(r"_cluster\d+", "", regex=True)

        for _, row in df.iterrows():
            pair = (row["source"], row["target"])
            edges.setdefault(pair, []).append(row["weight"])

    consensus = []
    for (tf, target), weights in edges.items():
        arr = np.array(weights)

        # maximum absolute value
        max_abs = np.max(np.abs(arr))

        # sign product (ignoring zeros)
        signs = np.sign(arr[arr != 0])
        if len(signs) == 0:
            sign = 0
        else:
            sign = np.prod(signs)

        consensus_val = max_abs * sign
        consensus.append((tf, target, consensus_val))

    consensus_df = pd.DataFrame(consensus, columns=["source", "target", "weight"])
    consensus_df.to_csv(output_file, index=False)

    return consensus_df

# Example usage:
consensus = merge_networks(path_outfile)
