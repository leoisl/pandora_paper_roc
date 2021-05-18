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
import json

def dump_dict(dictionary, filename):
    json_str = json.dumps(dictionary)
    with open(f"{filename}.json","w") as fout:
        fout.write(json_str)


# In[2]:


##################################################################################################################
# configs
reports_tsv_glob_path = "/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_paper_roc/analysis_output_pandora_paper_tag1/recall/reports/*/100x/pandora_nanopore_withdenovo/coverage_filter_0/strand_bias_filter_0.0/gaps_filter_1.0/gt_conf_percentile_0/*.tsv"
samples = ['063_STEC', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_8G', 'H131800734', 'CFT073', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_9D', 'ST38', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB2_1A']
gene_localisation_dir="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_gene_distance/gene_distance_pandora_paper_tag1/genes_from_truth_or_ref/pandora_nanopore_100x_withdenovo"
pandora_multisample_vcf_ref="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_workflow/pandora_output_pandora_paper_tag1/nanopore/100x/random/compare_withdenovo/pandora_multisample.vcf_ref.fa"
##################################################################################################################


# In[3]:


# get the report
def keep_report(filepath):
    for sample in samples:
        if sample in filepath:
            return True
    return False

reports_filepaths = glob(reports_tsv_glob_path)
reports_filepaths = [filepath for filepath in reports_filepaths if keep_report(filepath)]
recall_report = RecallReport.from_files(reports_filepaths,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)
recall_report.report


# In[4]:


def get_value_from_header_fast(header: str, field: str, return_type, value_to_return_if_not_found, delim: str = ";"):
    try:
        string_with_value = header[header.index(field) + len(field) + 1:]  # "+ 1" is for the "="
    except ValueError:
        assert False

    try:
        string_with_value = string_with_value[:string_with_value.index(delim)]
    except ValueError:
        assert False

    try:
        return return_type(string_with_value)
    except ValueError:
        assert False

def create_field_from_header(report, field:str, probe_header: str, field_type, default_value) -> None:
    report[field] = report[probe_header].apply(
        lambda header: get_value_from_header_fast(header, field, field_type, default_value)
    )

# improve recall report with the original pos
create_field_from_header(recall_report.report, "OR_STRAND", "query_probe_header", str, None)
create_field_from_header(recall_report.report, "OR_POS", "query_probe_header", int, None)
create_field_from_header(recall_report.report, "CHROM", "query_probe_header", str, None)
create_field_from_header(recall_report.report, "POS", "query_probe_header", int, None)

recall_report.report.to_csv("report.csv")
recall_report.report


# In[5]:


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


# In[7]:


def find_genes(chrom, sample, pos):
    genes = set()
    for gene in sample_to_chrom_to_intervaltree[sample][chrom][pos]:
        gene = gene.data
        genes.add(gene)
    return genes

pv_to_genes = defaultdict(set)
for chrom, sample, pos, pv_id in zip(recall_report.report["CHROM"], recall_report.report["sample"], recall_report.report["POS"], recall_report.report["PVID"]):
    pv_to_genes[pv_id].update(find_genes(chrom, sample, pos))

for pv in pv_to_genes:
    pv_to_genes[pv] = list(pv_to_genes[pv])

dump_dict(pv_to_genes, "pv_to_genes.json")

pv_to_genes


# In[8]:


nb_of_genes_in_pvs = list(map(len, pv_to_genes.values()))
[[elem, nb_of_genes_in_pvs.count(elem)] for elem in set(nb_of_genes_in_pvs)]


# In[9]:


plot = sns.histplot(nb_of_genes_in_pvs)
fig = plot.get_figure()
fig.savefig("panvars_and_nb_of_genes.png")

