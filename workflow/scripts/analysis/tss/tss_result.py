#!/usr/bin/env python
# coding: utf-8

# In[42]:


import pandas as pd
from concurrent.futures import ThreadPoolExecutor
from itertools import combinations
from tqdm import tqdm
import pyranges as pr
import argparse


# Parse args

parser = argparse.ArgumentParser()
parser.add_argument('-g', '--path_granie', required=True)
parser.add_argument('-m', '--path_hummus', required=True)
parser.add_argument('-p', '--path_pando', required=True)
parser.add_argument('-s', '--path_scenicplus', required=True)
parser.add_argument('-c', '--path_celloracle', required=True)
parser.add_argument('-d', '--path_dictys', required=True)
parser.add_argument('-f', '--path_figr', required=True)
parser.add_argument('-o', '--path_out', required=True)
args = parser.parse_args()
granie_path = args.path_granie
hummus_path = args.path_hummus
scenicplus_path = args.path_scenicplus
celloracle_path = args.path_celloracle
figr_path = args.path_figr
pando_path = args.path_pando
dictys_path = args.path_dictys
out_path = args.path_out


# Read data
hummus = pd.read_csv(hummus_path, index_col=0)
scenicplus = pd.read_csv(scenicplus_path, index_col=0)
granie = pd.read_csv(granie_path, index_col=0)
celloracle = pd.read_csv(celloracle_path)
figr = pd.read_csv(figr_path, index_col=0)
pando = pd.read_csv(pando_path, index_col=0)
dictys = pd.read_csv(dictys_path, sep="\t", index_col=0)


# Define the transform function to convert and merge datasets
def transform(data):
    pyranges_obj = pr.PyRanges(data)
    merged = pyranges_obj.merge(by="Name")
    return merged


# Define the function that calculates overlap coefficient per gene
def overlap_coef_per_gene(gene, data1, data2):
    filt1 = data1[data1.Name == gene]
    filt2 = data2[data2.Name == gene]
    
    if filt1.empty or filt2.empty:
        return {'Gene': gene, 'Overlap_Coef': 0}
    overlap = filt1.intersect(filt2)
    
    # Calculate overlap coefficient
    if overlap.empty:
        overlap_coef = 0  # No overlap
    else:
        if overlap.length == 0:
            overlap_coef = 1
        else:
            overlap_coef = overlap.length / min(filt1.length, filt2.length)
    
    return {'Gene': gene, 'Overlap_Coef': overlap_coef}

# Define the function to parallelize
def calculate_all_overlap_coef(datasets, dataset_names, max_workers=32):
    transformed_datasets = [transform(data) for data in datasets]
    all_results = pd.DataFrame(columns=["Annotation1", "Annotation2", "Gene", "Overlap_Coef"])

    # Corrected indentation here
    for (data1, name1), (data2, name2) in combinations(zip(transformed_datasets, dataset_names), 2):
        intersection = set(data1.Name.unique()).intersection(set(data2.Name.unique()))
        results = []

        with ThreadPoolExecutor(max_workers=max_workers) as executor:
            futures = {executor.submit(overlap_coef_per_gene, gene, data1, data2): gene for gene in intersection}
            for future in tqdm(futures, desc=f"Calculating Overlap Coefficient for {name1} vs {name2}"):
                result = future.result()
                if result is not None:
                    results.append(result)

        if results:  
            overlap_coef_df = pd.DataFrame(results)
            overlap_coef_df["Annotation1"] = name1
            overlap_coef_df["Annotation2"] = name2

            all_results = pd.concat([all_results, overlap_coef_df], ignore_index=True)

    return all_results

datasets = [scenicplus, granie, hummus, celloracle, figr, pando, dictys]
dataset_names = ["scenicplus", "granie", "hummus", "celloracle", "figr", "pando", "dictys"]

# Calculate overlap coef
all_results_df = calculate_all_overlap_coef(datasets, dataset_names)

# Save the objects
all_results_df.to_csv(out_path, index=False)




