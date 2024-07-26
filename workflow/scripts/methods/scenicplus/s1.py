import os                                                                                                                                                                                                                
import pickle                                                                                                                                                                                                            
import numpy as np                                                                                                                                                                                                       
import pandas as pd                                                                                                                                                                                                      
from scipy.sparse import csr_matrix                                                                                                                                                                                      
import argparse                                                                                                                                                                                                          
import sys
import pandas as pd
import requests
import pickle
import dill
import warnings
import pyranges as pr
import scanpy as sc
import pybiomart as pbm
import scrublet as scr
import pycisTopic
import scenicplus
from pycisTopic.cistopic_class import *
from pycisTopic.lda_models import *
from pycisTopic.diff_features import impute_accessibility, CistopicImputedFeatures, normalize_scores, find_highly_variable_features, find_diff_features
from pycisTopic.topic_binarization import *
import pyranges as pr
from pycistarget.utils import region_names_to_coordinates
from scenicplus.wrappers.run_pycistarget import run_pycistarget

## sub motif enrichment (later)
def get_genome_annotation(db_dir, filename="gene_annotation_v1.csv"):
    return(pd.read_csv(os.path.join(db_dir, filename)))


def download_annotations(db_dir, biomart_host="http://sep2019.archive.ensembl.org/", specie='mouse'):
    if specie == 'human':
        name_biomart = 'hsapiens_gene_ensembl'
        tf_filename = "hg38_tfs.txt"
    elif specie == 'mouse':
        name_biomart = 'mmusculus_gene_ensembl'
        tf_filename = "mm10_tfs.txt"
    make_genome_annotation(db_dir, biomart_host, name_biomart = name_biomart)
    make_pr_annot(db_dir, biomart_host)
    # we need a file with tf names -> generate it from motif annotation
    make_tf_file(db_dir, tf_filename=tf_filename, specie=specie)

## sub download
def make_genome_annotation(db_dir, biomart_host, filename="gene_annotation_v1.csv", name_biomart='mmusculus_gene_ensembl'):
    full_path = os.path.join(db_dir, filename)
    if not os.path.exists(full_path):
        dataset = pbm.Dataset(name=name_biomart,  host=biomart_host)
        annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
        annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].to_numpy(dtype = str)
        filter = annot['Chromosome/scaffold name'].str.contains('CHR|GL|JH|MT')
        annot = annot[~filter]
        annot['Chromosome/scaffold name'] = annot['Chromosome/scaffold name'].str.replace(r'(\b\S)', r'chr\1')
        annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
        annot = annot[annot.Transcript_type == 'protein_coding']
        annot.to_csv(full_path)

## sub download
def make_pr_annot(db_dir, biomart_host, filename="gene_annotation_v2.csv"):
    #`pr.PyRanges` containing gene annotation, including
    # Chromosome, Start, End, Strand (as '+' and '-'), Gene name and Transcription Start Site.
    full_path = os.path.join(db_dir, filename)
    if not os.path.exists(full_path):
        dataset = pbm.Dataset(name='mmusculus_gene_ensembl',  host=biomart_host)
        if 'external_gene_name' not in dataset.attributes.keys():
            external_gene_name_query = 'hgnc_symbol'
        else:
            external_gene_name_query = 'external_gene_name'
        if 'transcription_start_site' not in dataset.attributes.keys():
            transcription_start_site_query = 'transcript_start'
        else:
            transcription_start_site_query = 'transcription_start_site'
        annot = dataset.query(attributes=[
            'chromosome_name', 'start_position', 'end_position', 'strand',
            external_gene_name_query, transcription_start_site_query, 'transcript_biotype'
        ])
        annot.columns = ['Chromosome', 'Start', 'End', 'Strand', 'Gene', 'Transcription_Start_Site', 'Transcript_type']
        annot['Chromosome'] = 'chr' + annot['Chromosome'].astype(str)
        annot = annot[annot.Transcript_type == 'protein_coding']
        annot.Strand[annot.Strand == 1] = '+'
        annot.Strand[annot.Strand == -1] = '-'
        annot.to_csv(full_path)

## sub download
def make_tf_file(db_dir, tf_filename="mm10_tfs.txt", specie='mouse'):
    full_path = os.path.join(db_dir, tf_filename)
    if not os.path.exists(full_path):
        if specie == 'mouse':
            motif_annotation = os.path.join(db_dir, 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl')
        elif specie == 'human':
            motif_annotation = os.path.join(db_dir, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')
        annot = pd.read_csv(motif_annotation, sep="\t")
        tfs = set(str(g) for g in annot["gene_name"])
        with open(full_path, 'w') as f:
            for line in tfs:
                f.write(f"{line}\n")    

## sub motif enrichment (later)
def get_genome_annotation(db_dir, filename="gene_annotation_v1.csv"):
    return(pd.read_csv(os.path.join(db_dir, filename)))


def get_pr_annot(db_dir, filename="gene_annotation_v2.csv"):
    full_path = os.path.join(db_dir, filename)
    result = pd.read_csv(full_path)
    result = pr.PyRanges(result.dropna(axis=0))
    return(result)

def get_chr_size(db_dir):
    target_file=os.path.join(db_dir, 'mm10.chrom.sizes')
    result=pd.read_csv(target_file, sep='\t', header=None)
    result.columns=['Chromosome', 'End']
    result['Start']=[0]*result.shape[0]
    result=result.loc[:,['Chromosome', 'Start', 'End']]
    # Exceptionally in this case, to agree with CellRangerARC annotations
    # adapted to mm10 (originally for hg38)
    result['Chromosome'] = [result['Chromosome'][x].split('_')[1] + ".1" if len(result['Chromosome'][x].split('_')) > 1 else result['Chromosome'][x]
                            for x in range(len(result['Chromosome']))]
    result=pr.PyRanges(result)
    return(result)

def blacklist_file(db_dir):
    return(os.path.join(db_dir, "mm10-blacklist.v2.bed.gz"))

def tf_file(db_dir, filename="mm10_tfs.txt"):
    return(os.path.join(db_dir, filename))

def find_best_ensembl_host(scplus_obj, species='mmusculus'):
    ensembl_version_dict = {
        '105': 'http://www.ensembl.org', '104': 'http://may2021.archive.ensembl.org/', '103': 'http://feb2021.archive.ensembl.org/',
        '102': 'http://nov2020.archive.ensembl.org/', '101': 'http://aug2020.archive.ensembl.org/', '100': 'http://apr2020.archive.ensembl.org/',
        '99': 'http://jan2020.archive.ensembl.org/', '98': 'http://sep2019.archive.ensembl.org/', '97': 'http://jul2019.archive.ensembl.org/',
        '96': 'http://apr2019.archive.ensembl.org/', '95': 'http://jan2019.archive.ensembl.org/', '94': 'http://oct2018.archive.ensembl.org/',
        '93': 'http://jul2018.archive.ensembl.org/', '92': 'http://apr2018.archive.ensembl.org/', '91': 'http://dec2017.archive.ensembl.org/',
        '90': 'http://aug2017.archive.ensembl.org/'
        #'89': 'http://may2017.archive.ensembl.org/', '88': 'http://mar2017.archive.ensembl.org/',
        #'87': 'http://dec2016.archive.ensembl.org/', '86': 'http://oct2016.archive.ensembl.org/', '80': 'http://may2015.archive.ensembl.org/',
        #'77': 'http://oct2014.archive.ensembl.org/', '75': 'http://feb2014.archive.ensembl.org/', '54': 'http://may2009.archive.ensembl.org/'
    }
    n_overlap = {}
    for version in ensembl_version_dict.keys():
        print(f'host: {version}')
        try:
            n_overlap[version] = test_ensembl_host(scplus_obj, ensembl_version_dict[version], species)
        except:
            print('Host not reachable')
    v = sorted(n_overlap.items(), key=lambda item: item[1], reverse=True)[0][0]
    print(f"version: {v} has the largest overlap, use {ensembl_version_dict[v]} as biomart host")
    return(ensembl_version_dict[v])

def test_ensembl_host(scplus_obj, host, species):
    dataset = pbm.Dataset(name=species+'_gene_ensembl',  host=host)
    annot = dataset.query(attributes=['chromosome_name', 'transcription_start_site', 'strand', 'external_gene_name', 'transcript_biotype'])
    annot.columns = ['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    annot['Chromosome'] = annot['Chromosome'].astype('str')
    filter = annot['Chromosome'].str.contains('CHR|GL|JH|MT')
    annot = annot[~filter]
    annot.columns=['Chromosome', 'Start', 'Strand', 'Gene', 'Transcript_type']
    gene_names_release = set(annot['Gene'].tolist())
    ov=len([x for x in scplus_obj.gene_names if x in gene_names_release])
    print('Genes recovered: ' + str(ov) + ' out of ' + str(len(scplus_obj.gene_names)))
    return ov

# Init args
parser = argparse.ArgumentParser()

parser.add_argument('-r', '--rna_path', required=True)
parser.add_argument('-a', '--atac_path', required=True)
parser.add_argument('-o', '--output_folder', required=True)
parser.add_argument('-s', '--species', required=True)

args = vars(parser.parse_args())

rna_path = args['rna_path']
atac_path = args['atac_path']
output_folder = args['output_folder']
species = args['species']

work_dir = output_folder
temp_dir = os.path.join(work_dir,'/tmp')

# Init paths to databases
if species == 'human':
    db_dir = 'human_data'
    rankings_db = os.path.join(db_dir, 'hg38_screen_v10_clust.regions_vs_motifs.rankings.feather')
    scores_db = os.path.join(db_dir, 'hg38_screen_v10_clust.regions_vs_motifs.scores.feather')
    motif_annotation = os.path.join(db_dir, 'motifs-v10nr_clust-nr.hgnc-m0.001-o0.0.tbl')

elif species == 'mouse':
    db_dir = 'mouse_data'
    rankings_db = os.path.join(db_dir, 'mm10_screen_v10_clust.regions_vs_motifs.rankings.feather')
    scores_db = os.path.join(db_dir, 'mm10_screen_v10_clust.regions_vs_motifs.scores.feather')
    motif_annotation = os.path.join(db_dir, 'motifs-v10nr_clust-nr.mgi-m0.001-o0.0.tbl')

print('downloading annotations')
download_annotations(db_dir, specie=species)

# Load atac
print('load atac data')
peak_mtx = pd.read_csv(atac_path, sep='\t')

##########################
#### Run cisTopic #######
##########################

print("init cisTopic object")
# Init cistopic object
cistopic_obj = CistopicObject(
    fragment_matrix=csr_matrix((peak_mtx > 0).astype(np.int_)),
    binary_matrix=csr_matrix((peak_mtx > 0).astype(np.int_)),
    cell_names=peak_mtx.columns,
    region_names=peak_mtx.index,
    cell_data=pd.DataFrame(index=peak_mtx.columns),
    region_data=pd.DataFrame(index=peak_mtx.index),
    path_to_fragments="",
    project = "cistopic_mouse"
    )

print('LDA model of atac data')
models = run_cgs_models(
    cistopic_obj,
    n_topics=[2,4,10,16,32,48],
    n_cpu=5,
    n_iter=500,
    random_state=555,
    alpha=50,
    alpha_by_topic=True,
    eta=0.1,
    eta_by_topic=False,
    save_path=None,
    _temp_dir = temp_dir,)


# Model accessibility through LDA ?
if not os.path.exists(os.path.join(work_dir, 'scATAC/models')):
    os.makedirs(os.path.join(work_dir, 'scATAC/models'))

pickle.dump(models,
            open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'wb'))

models = pickle.load(open(os.path.join(work_dir, 'scATAC/models/10x_pbmc_models_500_iter_LDA.pkl'), 'rb'))
#cistopic_obj = pickle.load(open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'rb'))
model = evaluate_models(models,
                       select_model=16,
                       return_model=True,
                       metrics=['Arun_2010','Cao_Juan_2009', 'Minmo_2011', 'loglikelihood'],
                       plot_metrics=False)

# Add LDA model to cistopic object
print('Add LDA model to cisTopic object')
cistopic_obj.add_LDA_model(model)

# Save cistopic object
print("Save cistopic object")
pickle.dump(cistopic_obj,
            open(os.path.join(work_dir, 'scATAC/cistopic_obj.pkl'), 'wb'))

# Impute accessibility on regions to correct sparsity
print('Impute sparse accessibility')
imputed_obj = impute_accessibility(
        cistopic_obj,
        selected_cells=None,
        selected_regions=None,
        scale_factor=10**6,
)

# Save imputed accessibility 
print('Save imputed accessibility')
pickle.dump(imputed_obj,
            open(os.path.join(work_dir, 'scATAC/imputed_accessibility.pkl'), 'wb'))

# Select highly accessible binarized regions
print('Select important regions')
region_bin_topics_otsu = binarize_topics(cistopic_obj, method='otsu')
region_bin_topics_top3k = binarize_topics(cistopic_obj, method='ntop', ntop = 3000)

if not os.path.exists(os.path.join(work_dir, 'scATAC/candidate_enhancers')):
    os.makedirs(os.path.join(work_dir, 'scATAC/candidate_enhancers'))
pickle.dump(region_bin_topics_otsu, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'wb'))
pickle.dump(region_bin_topics_top3k, open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'wb'))


# open and format as a dictionary the regions sets
region_bin_topics_otsu = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_otsu.pkl'), 'rb'))
region_bin_topics_top3k = pickle.load(open(os.path.join(work_dir, 'scATAC/candidate_enhancers/region_bin_topics_top3k.pkl'), 'rb'))
region_sets = {}
region_sets['topics_otsu'] = {}
region_sets['topics_top_3'] = {}

for topic in region_bin_topics_otsu.keys():
    regions = region_bin_topics_otsu[topic].index[region_bin_topics_otsu[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_otsu'][topic] = pr.PyRanges(region_names_to_coordinates(regions))
for topic in region_bin_topics_top3k.keys():
    regions = region_bin_topics_top3k[topic].index[region_bin_topics_top3k[topic].index.str.startswith('chr')] #only keep regions on known chromosomes
    region_sets['topics_top_3'][topic] = pr.PyRanges(region_names_to_coordinates(regions))


#######################
# Run cisTarget
#######################
# /!\ parallel computing is a disaster, keep n_cpu=1
# /!\ loading the databases on the cluster is very slow (~10 minutes) -> consider using appa partition
# /!\ we need to bypass all biomart calls -> use custom_annot instead (set species to custom, too)
enhancer_dir = os.path.join(work_dir, 'scATAC', "candidate_enhancers")
motif_dir = os.path.join(work_dir, 'scATAC', "motifs")
if not os.path.exists(motif_dir):
    os.makedirs(motif_dir)

print('run cisTarget')
run_pycistarget(
    region_sets=region_sets,
    species='custom',
    custom_annot=get_genome_annotation(db_dir,),
    save_path=motif_dir,
    ctx_db_path=rankings_db,
    dem_db_path=scores_db,
    path_to_motif_annotations=motif_annotation,
    run_without_promoters=True,
    n_cpu=1,
    _temp_dir=temp_dir,
    annotation_version='v10nr_clust',
)
