#!/usr/bin/env python
# coding: utf-8

# In[22]:


import pybiomart as pbm
import argparse
import numpy as np

parser = argparse.ArgumentParser()
parser.add_argument('-j', '--path_out', required=True)
args = parser.parse_args()
out_path = args.path_out

dataset = pbm.Dataset(name='hsapiens_gene_ensembl', host='http://www.ensembl.org')

annot = dataset.query(attributes=['chromosome_name', 'start_position', 'end_position',
                                  'strand', 'external_gene_name', 'transcription_start_site', 'transcript_biotype'])
annot['chromosome_name'] = annot['chromosome_name'].to_numpy(dtype=str)
filter = annot['chromosome_name'].str.contains('CHR|GL|JH|MT', case=False)
annot = annot[~filter]
annot['chromosome_name'] = annot['chromosome_name'].str.replace(r'(\b\S)', r'chr\1')
annot.columns = ['Chromosome', 'Start', 'End', 'Strand', 'Name', 'Transcription_Start_Site', 'Transcript_type']
annot["Strand"] = annot["Strand"].replace({1: "+", -1: "-"})
annot.Start = annot.Start.astype(np.int32)
annot['Chromosome'] = 'chr' + annot['Chromosome'].astype(str)
annot.dropna(inplace=True)
annot = annot[['Chromosome', 'Start', 'End', 'Name']]

# Save the file
annot.to_csv(out_path, sep="\t", index=False)



# In[ ]:





# In[ ]:





# In[ ]:





# In[ ]:




