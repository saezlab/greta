#!/usr/bin/env python3
import sys, csv, re

"""
Convert a list of Ensembl mouse gene IDs (one per line) to gene symbols
using a simple two-column CSV mapping: id,symbol

Usage:
    lambert_from_ensembl.py ensembl.csv ids.txt output.csv
"""

def main():
    if len(sys.argv) != 4:
        sys.stderr.write("Usage: lambert_from_ensembl.py ensembl.csv ids.txt output.csv\n")
        sys.exit(2)

    ensembl_csv, ids_txt, out_csv = sys.argv[1:]

    # Build mapping: Ensembl ID -> symbol
    mapping = {}
    with open(ensembl_csv, newline="") as f:
        reader = csv.DictReader(f)
        if set(reader.fieldnames or []) != {"id", "symbol"}:
            sys.stderr.write(
                f"ERROR: Expected headers 'id,symbol' in {ensembl_csv}, "
                f"found {reader.fieldnames}\n"
            )
            sys.exit(1)
        for row in reader:
            gid = (row["id"] or "").strip()
            gid = re.sub(r"\.\d+$", "", gid)  # drop version if present
            sym = (row["symbol"] or "").strip()
            if gid and sym and gid not in mapping:
                mapping[gid] = sym

    # Read Ensembl IDs, map to symbols (preserve order, drop dups)
    seen_syms = set()
    symbols = []
    missing = 0
    with open(ids_txt) as f:
        for line in f:
            gid = line.strip()
            if not gid:
                continue
            gid = re.sub(r"\.\d+$", "", gid)
            sym = mapping.get(gid)
            if sym:
                if sym not in seen_syms:
                    seen_syms.add(sym)
                    symbols.append(sym)
            else:
                missing += 1

    # Write one symbol per line
    with open(out_csv, "w") as out:
        for s in symbols:
            out.write(f"{s}\n")

    if missing:
        sys.stderr.write(f"[lambert_from_ensembl] Warning: {missing} IDs had no symbol mapping.\n")

if __name__ == "__main__":
    main()