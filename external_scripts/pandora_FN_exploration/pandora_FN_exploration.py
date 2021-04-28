#!/usr/bin/env python
# coding: utf-8

# In[1]:


# TODO: RELATIVE POSITION OF PANVARS INSIDE IDENTIFIED GENES

import pandas as pd
from glob import glob
from evaluate.report import Report, PrecisionReport, RecallReport
from collections import defaultdict, namedtuple
import re
import intervaltree
import seaborn as sns
sns.set()
import json

def dump_dict(dictionary, filename):
    json_str = json.dumps(dictionary)
    with open(f"{filename}.json","w") as fout:
        fout.write(json_str)


# In[62]:


##################################################################################################################
# configs
# cluster
reports_tsv_glob_path = "/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_paper_roc/analysis_output_pandora_paper_tag1/recall/reports/*/100x/pandora_nanopore_withdenovo/coverage_filter_0/strand_bias_filter_0.0/gaps_filter_1.0/gt_conf_percentile_0/*.tsv"
samples = ['063_STEC', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_8G', 'H131800734', 'CFT073', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_9D', 'ST38', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB2_1A']
gene_localisation_dir="/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_gene_distance/gene_distance_pandora_paper_tag1/genes_from_truth_or_ref/pandora_nanopore_100x_withdenovo"
sample_data_dir="/hps/nobackup/iqbal/leandro/pdrv/data/samples/*"
probesets_dir="/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pangenome_variations/output_pangenome_variations_pandora_paper_tag1/truth_probesets"
nb_of_threads = 96
##################################################################################################################


# In[3]:


# get all refs first
from Bio import SeqIO
from collections import defaultdict
from glob import glob
from pathlib import Path
ref_to_chrom_to_seq = defaultdict(lambda: defaultdict(str))
refs = glob(f"{sample_data_dir}/*.ref.fa")
for ref in refs:
    for record in SeqIO.parse(ref, "fasta"):
        ref = Path(ref).name.replace(".ref.fa", "")
        ref_to_chrom_to_seq[ref][record.id] = record.seq
ref_to_chrom_to_seq


# In[4]:


# get the report
reports_filepaths = glob(reports_tsv_glob_path)
recall_report = RecallReport.from_files(reports_filepaths,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)
recall_report.report


# In[5]:


variation_found_nbofsamples = recall_report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(
    binary=True
)
variation_found_nbofsamples.rename(columns={"proportion_of_allele_seqs_found_binary": "FOUND"}, inplace=True)
variation_found_nbofsamples.reset_index(inplace=True)
variation_found_nbofsamples


# In[33]:


# get a df with the info needed for the pvs
chroms_array = []
samples_array = []
positions_fw_array = []
positions_rc_array = []
pv_ids_array = []
for header in recall_report.report["query_probe_header"]:
    chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", header)[0]
    chroms_array.append(chrom)
    samples_array.append(sample)
    pos_fw = int(pos)
    positions_fw_array.append(pos_fw)
    chrom_len = len(ref_to_chrom_to_seq[sample][chrom])
    pos_rc = chrom_len-pos_fw-1
    positions_rc_array.append(pos_rc)
    pv_ids_array.append(int(pv_id))
pvs_df = pd.DataFrame(data={"chrom": chroms_array, "sample": samples_array, "pos_fw": positions_fw_array, "pos_rc": positions_rc_array, "pv_id": pv_ids_array})
pvs_df.to_csv("pvs_df.csv")
pvs_df


# In[28]:


# create an index of probe seqs
pv_id_sample_chrom_to_probe = {}
for probeset_file in glob(f"{probesets_dir}/*/*.fa"):
    for record in SeqIO.parse(probeset_file, "fasta"):
        header = record.id 
        chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", header)[0]
        key = (int(pv_id), sample, chrom)
        pv_id_sample_chrom_to_probe[key] = record.seq
pv_id_sample_chrom_to_probe


# In[34]:


# improve pvs_df with probe seqs
probes=[]
for index, (chrom, sample, pos_fw, pos_rc, pv_id) in pvs_df.iterrows():
    probes.append(str(pv_id_sample_chrom_to_probe[(pv_id, sample, chrom)]))
pvs_df["probes"]=probes
pvs_df.to_csv("pvs_df.csv")
pvs_df


# In[7]:


sample_to_chrom_to_intervaltree = defaultdict(lambda: defaultdict(lambda: intervaltree.IntervalTree()))
for sample in samples:
    filename = gene_localisation_dir+f"/{sample}.csv"
    with open(filename) as gene_localisation_fh:
        for line in gene_localisation_fh:
            status,gene_name,ref_or_truth_id,contig,start,stop,sequence = line.split(",")
            if status=="Mapped":
                start, stop = int(start), int(stop)
                sample_to_chrom_to_intervaltree[ref_or_truth_id][contig][start:stop] = gene_name
sample_to_chrom_to_intervaltree


# In[8]:


# These plots show that gene annotation is fine
def get_df_of_how_many_genes_in_each_pos_of_chrom(sample, chrom):
    nb_of_pos_with_x_genes = defaultdict(int)
    for pos in range(len(ref_to_chrom_to_seq[sample][chrom])):
        nb_of_genes_in_pos = len(sample_to_chrom_to_intervaltree[sample][chrom][pos])
        nb_of_pos_with_x_genes[nb_of_genes_in_pos]+=1
    
    
    counts = []
    nbs_of_genes = []
    for nb_of_genes, count in nb_of_pos_with_x_genes.items():
        nbs_of_genes.append(nb_of_genes)
        counts.append(count)
    
    chrom_desc = "chrom" if chrom=="0" else f"pl{chrom}"
    sample = sample.replace("Escherichia_coli_", "")
    desc = f"{sample} {chrom_desc}"
    descs = [desc] * len(counts)
    genes_in_each_pos_df = pd.DataFrame(data={'desc': descs, 'nb_of_genes': nbs_of_genes, 'count': counts})
    return genes_in_each_pos_df

dfs=[]
for sample in sample_to_chrom_to_intervaltree:
    chrom_to_intervaltree = sample_to_chrom_to_intervaltree[sample]
    for chrom in chrom_to_intervaltree:
        dfs.append(get_df_of_how_many_genes_in_each_pos_of_chrom(sample, chrom))
        
genes_in_each_pos_df = pd.concat(dfs, ignore_index=True)
genes_in_each_pos_df.to_csv("genes_in_each_pos_df.csv")
genes_in_each_pos_df


# In[9]:


plot = sns.FacetGrid(genes_in_each_pos_df, col="desc", col_wrap=8)
plot.map(sns.barplot, "nb_of_genes", "count").set(yscale = 'log')
for ax in plot.axes.flatten():
    ax.tick_params(labelbottom=True)
plot.savefig("gene_frequency_distribution.png")


# In[46]:


gene_to_sequence = {}
with open("data/pandora_multisample.vcf_ref.fa") as fin:
    for line in fin:
        line = line.strip()
        if line[0] == ">":
            gene=line[1:]
        else:
            seq=line
            gene_to_sequence[gene] = seq
gene_to_sequence


# In[66]:


# choose genes that barely match the probe

import swalign
from Bio.Seq import Seq

# choose your own values here… 2 and -1 are common.
match = 2
mismatch = -1
scoring = swalign.NucleotideScoringMatrix(match, mismatch)
sw = swalign.LocalAlignment(scoring)  # you can also choose gap penalties, etc...

def probe_and_gene_has_some_similarity(probe, gene):
    alignment = sw.align(probe, gene)
    score_fw = alignment.score
    
    probe_rc = str(Seq(probe).reverse_complement())
    alignment = sw.align(probe_rc, gene)
    score_rc = alignment.score
    
    score = max(score_fw, score_rc)
    score_threshold = len(probe) * match * 0.5  # at least 50% matches  
    return score > score_threshold
        

def find_genes(chrom, sample, pos_fw, pos_rc, pv_id, probe):
    genes = set()
    for pos in [pos_fw, pos_rc]:
        for gene in sample_to_chrom_to_intervaltree[sample][chrom][pos]:
            gene = gene.data
            gene_seq = gene_to_sequence[gene]
            if probe_and_gene_has_some_similarity(probe, gene_seq):
                genes.add(gene)
    return pv_id, genes


params = []
for chrom, sample, pos_fw, pos_rc, pv_id, probe in zip(pvs_df["chrom"], pvs_df["sample"], pvs_df["pos_fw"], pvs_df["pos_rc"], pvs_df["pv_id"], pvs_df["probes"]):
    params.append((chrom, sample, pos_fw, pos_rc, pv_id, probe))

import multiprocessing
with multiprocessing.Pool(processes=nb_of_threads) as pool:
    pv_ids_and_genes = pool.starmap(find_genes, params)


pv_to_genes = {}
for pv_id, genes in pv_ids_and_genes:
    pv_to_genes[pv_id] = list(genes)  # set is not serializable
dump_dict(pv_to_genes, "pv_to_genes")
pv_to_genes



