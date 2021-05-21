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
import pickle

def dump_object_to_json(dictionary, filename):
    json_str = json.dumps(dictionary)
    with open(f"{filename}.json","w") as fout:
        fout.write(json_str)

def dump_object(obj, filename):
    with open(f"{filename}.pickle","w") as fout:
        pickle.dump(obj, fout)


# In[2]:


##################################################################################################################
# configs
# cluster
reports_tsv_glob_path = "/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_paper_roc/analysis_output_pandora_paper_tag1/recall/reports/*/100x/pandora_nanopore_withdenovo/coverage_filter_0/strand_bias_filter_0.0/gaps_filter_1.0/gt_conf_percentile_0/*.tsv"
samples = ['063_STEC', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_8G', 'H131800734', 'CFT073', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_9D', 'ST38', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB2_1A']
gene_localisation_dir="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_gene_distance/gene_distance_vcf_completed_with_prg/genes_from_truth_or_ref/pandora_vcf_ref_nanopore_100x_with_prg"
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


# In[52]:


# get a df with the info needed for the pvs
chroms_array = []
samples_array = []
positions = []
pv_ids_array = []
for header in recall_report.report["query_probe_header"]:
    chrom, sample, pv_id, pos = re.findall("CHROM=(.+?);SAMPLE=(.+?);.*;PVID=(.+?);.*;OR_POS=(.+?);.*", header)[0]
    chroms_array.append(chrom)
    samples_array.append(sample)
    pv_ids_array.append(int(pv_id))
    positions.append(int(pos))


pvs_df_leah = pd.DataFrame(data={
    "chrom": chroms_array,
    "sample": samples_array,
    "pos": positions,
    "pv_id": pv_ids_array,
    "eval": recall_report.report["good_eval"]})
pvs_df_leah.to_csv("pvs_df_leah.csv", index=False)
