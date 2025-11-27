import pandas as pd
from io import BytesIO
import sys

binary_data = BytesIO(sys.stdin.buffer.read())
df = pd.read_excel(pd.ExcelFile(binary_data), sheet_name=0)
# Mouse remap2022 Excel doesn't have BTO_id column, use biotype directly
df = df[['biotype']].dropna()
df['term'] = df['biotype']
df[['biotype', 'term']].to_csv(sys.argv[2], sep='\t', index=False, header=None)
