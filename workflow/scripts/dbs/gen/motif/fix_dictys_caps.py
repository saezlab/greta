#!/usr/bin/env python3
import sys

def fix_file(file_path):
    with open(file_path, "r") as f:
        lines = f.readlines()
    
    with open(file_path, "w") as f:
        for line in lines:
            if line.startswith(">"):
                parts = line[1:].rstrip("\n").split("\t")
                if len(parts) >= 2 and "_MOUSE" in parts[1]:
                    try:
                        gene, rest = parts[1].split("_MOUSE", 1)
                        parts[1] = gene.capitalize() + "_MOUSE" + rest
                        line = ">" + "\t".join(parts) + "\n"
                    except ValueError: pass
            f.write(line)

if __name__ == "__main__":
    if len(sys.argv) != 2:
        sys.exit(f"Usage: {sys.argv[0]} input.motif")
    fix_file(sys.argv[1])