#!/usr/bin/env python3
import argparse
import sys
import subprocess
import tempfile
import anndata as ad
import pandas as pd

# Init args
parser = argparse.ArgumentParser()
parser.add_argument('-a', '--frag_in', nargs='+', required=True,
                    help='Input fragment .tsv.gz files (one per sample)')
parser.add_argument('-b', '--rna_in', nargs='+', required=True,
                    help='Input RNA .h5ad files (one per sample)')
parser.add_argument('-c', '--annot_in', required=True,
                    help='Input annotation CSV')
parser.add_argument('-d', '--frag_out', nargs='+', required=True,
                    help='Output filtered fragment .tsv.gz files (one per sample)')
parser.add_argument('-e', '--rna_out', required=True,
                    help='Output merged filtered RNA .h5ad')
parser.add_argument('-f', '--annot_out', required=True,
                    help='Output filtered annotation CSV')
args = vars(parser.parse_args())

frag_in = args['frag_in']
rna_in = args['rna_in']
annot_in = args['annot_in']
frag_out = args['frag_out']
rna_out = args['rna_out']
annot_out = args['annot_out']

if not (len(frag_in) == len(rna_in) == len(frag_out)):
    print("[filter_merge] frag_in, rna_in and frag_out must have same length",
          file=sys.stderr)
    sys.exit(1)


# Read global annotation
annot = pd.read_csv(annot_in, index_col=0)
annot_barcodes = set(annot.index)

# Read RNA and collect global RNA barcodes
rna_list = []
rna_barcodes_all = set()
for path in rna_in:
    rna = ad.read_h5ad(path)
    rna_list.append(rna)
    rna_barcodes_all.update(rna.obs_names)

# Collect global fragment barcodes (union over all fragment files)
frag_barcodes_all = set()
for path in frag_in:
    cmd = f"zcat {path} | awk '!/^#/ {{print $4}}' | sort -u"
    result = subprocess.run(
        cmd,
        shell=True,
        check=True,
        capture_output=True,
        text=True,
    )
    if result.stdout:
        frag_barcodes_all.update(result.stdout.strip().splitlines())

# Common barcodes:
#  1) in annot
#  2) in at least one RNA file
#  3) in at least one fragment file
common = annot_barcodes & rna_barcodes_all & frag_barcodes_all
print(f"[filter_merge] Common barcodes before celltype filter: {len(common)}",
      file=sys.stderr)

# Apply global celltype filter (celltypes with >= 100 cells in annot.csv)
annot_common = annot.loc[list(common)].copy()
if 'celltype' in annot_common.columns and not annot_common.empty:
    counts = annot_common['celltype'].value_counts()
    keep_types = counts[counts >= 100].index
    drop_types = counts[counts < 100].index
    if len(drop_types) > 0:
        print(
            f"[filter_merge] Removing {len(drop_types)} celltypes < 100 cells "
            "based on annot.csv",
            file=sys.stderr,
        )
    annot_filtered = annot_common[annot_common['celltype'].isin(keep_types)].copy()
else:
    annot_filtered = annot_common

common_filtered = set(annot_filtered.index)
print(f"[filter_merge] Common barcodes after celltype filter: "
      f"{len(common_filtered)}", file=sys.stderr)

# Write filtered annotation
annot_filtered.to_csv(annot_out)

# Filter each RNA to common barcodes, then merge
filtered_rna_list = []
for i, rna in enumerate(rna_list, start=1):
    keep = [bc for bc in rna.obs_names if bc in common_filtered]
    print(f"[filter_merge] RNA file {i}: keeping {len(keep)} cells",
          file=sys.stderr)
    if not keep:
        continue
    filtered_rna_list.append(rna[keep, :].copy())

if not filtered_rna_list:
    print("[filter_merge][ERROR] No RNA cells left after filtering",
          file=sys.stderr)
    sys.exit(1)

merged_rna = ad.concat(filtered_rna_list, join="outer", merge="same")
merged_rna.write(rna_out)

print(
    f"[filter_merge] Wrote merged RNA: {rna_out} | "
    f"cells={merged_rna.n_obs}, genes={merged_rna.n_vars}",
    file=sys.stderr,
)

# Write common barcodes to a temporary file for AWK-based fragment filtering
with tempfile.NamedTemporaryFile(mode='w', delete=False, suffix='.txt') as f:
    bc_file = f.name
    for bc in sorted(common_filtered):
        f.write(bc + '\n')

awk_script = (
    f'BEGIN {{while ((getline < "{bc_file}") > 0) barcodes[$0] = 1}} '
    '/^#/ {print; next} '
    '$4 in barcodes {print}'
)

# Filter each fragment file individually, preserving per-sample outputs
for i, (frag_path, frag_out_path) in enumerate(zip(frag_in, frag_out), start=1):
    print(f"[filter_merge] Filtering fragments for file {i}", file=sys.stderr)

    with open(frag_out_path + '.tmp', 'w') as fout:
        zcat = subprocess.Popen(['zcat', frag_path], stdout=subprocess.PIPE)
        awk = subprocess.Popen(['awk', awk_script], stdin=zcat.stdout, stdout=fout)
        zcat.stdout.close()
        awk.communicate()

    subprocess.run(f"bgzip -c {frag_out_path}.tmp > {frag_out_path}",
                   shell=True, check=True)
    subprocess.run(['tabix', '-p', 'bed', frag_out_path], check=True)
    subprocess.run(f"rm {frag_out_path}.tmp", shell=True, check=True)

# Clean up barcode file
subprocess.run(f"rm {bc_file}", shell=True, check=True)

print(f"[filter_merge] Finished; total cells retained: {len(common_filtered)}",
      file=sys.stderr)
