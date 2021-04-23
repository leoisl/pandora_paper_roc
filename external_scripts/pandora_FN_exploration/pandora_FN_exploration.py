#!/usr/bin/env python
# coding: utf-8

# In[1]:


import pandas as pd
from glob import glob
from evaluate.report import Report, PrecisionReport, RecallReport
from collections import defaultdict, namedtuple
import re
import intervaltree
import seaborn as sns
sns.set()


# In[2]:


##################################################################################################################
# configs
reports_tsv_glob_path = "/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_paper_roc/analysis_output_pandora_paper_tag1/recall/reports/*/100x/pandora_nanopore_withdenovo/coverage_filter_0/strand_bias_filter_0.0/gaps_filter_1.0/gt_conf_percentile_0/*.tsv"
samples = ['063_STEC', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_8G', 'H131800734', 'CFT073', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_9D', 'ST38', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB2_1A']
gene_localisation_dir="/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_gene_distance/gene_distance_pandora_paper_tag1/genes_from_truth_or_ref/pandora_nanopore_100x_withdenovo"
##################################################################################################################


# In[3]:


reports_filepaths = glob(reports_tsv_glob_path)
recall_report = RecallReport.from_files(reports_filepaths,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)


# In[4]:


variation_found_nbofsamples = recall_report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(
    binary=True
)
variation_found_nbofsamples.rename(columns={"proportion_of_allele_seqs_found_binary": "FOUND"}, inplace=True)
variation_found_nbofsamples.reset_index(inplace=True)
variation_found_nbofsamples


# In[5]:


Coordinate = namedtuple('Coordinate', ['sample', 'chrom', 'pos'])
pv_allele_seq_id_to_coordinates = defaultdict(set)
for report_filepath in reports_filepaths:
    report_df = pd.read_csv(report_filepath, sep="\t")
    for query_probe_header in report_df.query_probe_header:
        chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", query_probe_header)[0]
        pos, pv_id, allele_seq_id = int(pos), int(pv_id), int(allele_seq_id)
        pv_allele_seq_id_to_coordinates[(pv_id, allele_seq_id)].add(Coordinate(sample, chrom, pos))
pv_allele_seq_id_to_coordinates


# In[6]:


sample_to_chrom_to_intervaltree = defaultdict(lambda: defaultdict(lambda: intervaltree.IntervalTree()))

for sample in samples:
    with open(gene_localisation_dir+f"/{sample}.csv") as gene_localisation_fh:
        for line in gene_localisation_fh:
            status,gene_name,ref_or_truth_id,contig,start,stop,sequence = line.split(",")
            if status=="Mapped":
                sample_to_chrom_to_intervaltree[ref_or_truth_id][contig][int(start) : int(stop)] = gene_name
sample_to_chrom_to_intervaltree


# In[7]:


def find_genes(coordinate):
    genes = []  # a list because a coordinate might have several genes going through it
    for gene in sample_to_chrom_to_intervaltree[coordinate.sample][coordinate.chrom][coordinate.pos]:
        genes.append(gene.data)
    return genes

pv_allele_seq_id_to_genes = defaultdict(list)
for pv_allele_seq_id, coordinates in pv_allele_seq_id_to_coordinates.items():
    for coordinate in coordinates:
        pv_allele_seq_id_to_genes[pv_allele_seq_id].extend(find_genes(coordinate))
pv_allele_seq_id_to_genes


# In[8]:


pv_allele_seq_id_to_genes_found_by_pandora = {}
for pv_allele_seq_id, gene_list in pv_allele_seq_id_to_genes.items():
    pv_allele_seq_id_to_genes_found_by_pandora[pv_allele_seq_id] = set(gene_list)
pv_allele_seq_id_to_genes_found_by_pandora


pv_ids = set()
for pv_allele_seq_id in pv_allele_seq_id_to_genes.keys():
    pv_ids.add(pv_allele_seq_id[0])

    
pv_id_to_genes_found_by_pandora = {}
for pv_id in pv_ids:
    pv_id_to_genes_found_by_pandora[pv_id] =         pv_allele_seq_id_to_genes_found_by_pandora.get((pv_id, 0), set()).union( 
        pv_allele_seq_id_to_genes_found_by_pandora.get((pv_id, 1), set()))
pv_id_to_genes_found_by_pandora

pv_id_to_nb_of_genes_found_by_pandora = {}
for pv_id, genes in pv_id_to_genes_found_by_pandora.items():
    pv_id_to_nb_of_genes_found_by_pandora[pv_id] = len(genes)
pv_id_to_nb_of_genes_found_by_pandora
    
pv_id_to_nb_of_genes_found_by_pandora_df = pd.DataFrame(data=pv_id_to_nb_of_genes_found_by_pandora.items(),
                                                 columns=["PANGENOME_VARIATION_ID", "NB_OF_GENES"])
pv_id_to_nb_of_genes_found_by_pandora_df["FOUND_IN_VCF_REF"] = (pv_id_to_nb_of_genes_found_by_pandora_df["NB_OF_GENES"] > 0).astype(int)
pv_id_to_nb_of_genes_found_by_pandora_df


# In[9]:


pv_ids_to_nb_of_genes_summary = pv_id_to_nb_of_genes_found_by_pandora_df[["PANGENOME_VARIATION_ID", "NB_OF_GENES"]].groupby("NB_OF_GENES").count()
plot = pv_ids_to_nb_of_genes_summary.plot(kind="bar")
fig = plot.get_figure()
fig.savefig("pandora_nb_of_genes_for_panvars.png")


# In[10]:


variation_found_nbofsamples = variation_found_nbofsamples.merge(pv_id_to_nb_of_genes_found_by_pandora_df)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==1)).astype(int)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==0)).astype(int)
variation_found_nbofsamples


# In[11]:


df = variation_found_nbofsamples[["NB_OF_SAMPLES", "FOUND", "NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF", "NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"]].groupby("NB_OF_SAMPLES").sum()
plot = df.plot(kind='bar', stacked=True)
fig = plot.get_figure()
fig.savefig("pandora_FN.png")

