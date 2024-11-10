import pandas as pd
from io import BytesIO
import sys

binary_data = BytesIO(sys.stdin.buffer.read())
df = pd.read_excel(pd.ExcelFile(binary_data), sheet_name=0)
df = df[['biotype', 'identifiants/0/BTO_id']].dropna()
df = df.rename(columns={'identifiants/0/BTO_id': 'id'})
df['id'] = df['id'].str.replace('_', ':')
bto = pd.read_csv(sys.argv[1], sep='\t', header=None).set_index(0)[1].to_dict()
df['term'] = [bto[i] for i in df['id']]
df[['biotype', 'term']].to_csv(sys.argv[2], sep='\t', index=False, header=None)
