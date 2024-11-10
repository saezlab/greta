import sys
import pandas as pd

b_dict = dict()
read = False

for line in sys.stdin:
    line = line.strip()
    if line.startswith('<owl:Class'):
        read = True
        continue
    elif line.startswith('<oboInOwl:id') and read:
        key = line.split('>')[1].split('<')[0]
        continue
    elif line.startswith('<rdfs:label') and read:
        val = line.split('>')[1].split('<')[0]
        continue
    elif line.startswith('</owl:Class>') and read:
        b_dict[key] = val
        read = False

b_dict = pd.DataFrame(list(b_dict.items()))
b_dict.to_csv(sys.argv[1], sep='\t', index=False, header=None)
