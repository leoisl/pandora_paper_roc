#!/usr/bin/env python
# coding: utf-8

# In[2]:


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


# In[3]:


##################################################################################################################
# configs
reports_tsv_glob_path = "/home/leandro/git/pandora_paper_roc/external_scripts/pandora_FN_exploration/data/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_paper_roc/analysis_output_pandora_paper_tag1/recall/reports/*/100x/pandora_nanopore_withdenovo/coverage_filter_0/strand_bias_filter_0.0/gaps_filter_1.0/gt_conf_percentile_0/*.tsv"
samples = ['063_STEC', 'Escherichia_coli_MINF_1D', 'Escherichia_coli_MINF_9A', 'Escherichia_coli_MSB1_4E', 'Escherichia_coli_MSB1_7A', 'Escherichia_coli_MSB1_8G', 'H131800734', 'CFT073', 'Escherichia_coli_MINF_7C', 'Escherichia_coli_MSB1_1A', 'Escherichia_coli_MSB1_4I', 'Escherichia_coli_MSB1_7C', 'Escherichia_coli_MSB1_9D', 'ST38', 'Escherichia_coli_MINF_1A', 'Escherichia_coli_MINF_8D', 'Escherichia_coli_MSB1_3B', 'Escherichia_coli_MSB1_6C', 'Escherichia_coli_MSB1_8B', 'Escherichia_coli_MSB2_1A']
samples = ['063_STEC', 'CFT073']  # TODO: remove
gene_localisation_dir="/home/leandro/git/pandora_paper_roc/external_scripts/pandora_FN_exploration/data/hps/nobackup/iqbal/leandro/pdrv/out_20_way/pandora_gene_distance/gene_distance_pandora_paper_tag1/genes_from_truth_or_ref/pandora_nanopore_100x_withdenovo"
##################################################################################################################


# In[4]:


reports_filepaths = glob(reports_tsv_glob_path)
recall_report = RecallReport.from_files(reports_filepaths,
                                        concatenate_dfs_one_by_one_keeping_only_best_mappings=True)


# In[19]:


chroms = []
samples = []
positions = []
pv_ids = []
for header in recall_report.report["query_probe_header"]:
    chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", header)[0]
    chroms.append(chrom)
    samples.append(sample)
    positions.append(pos)
    pv_ids.append(pv_id)
pvs_df = pd.DataFrame(data={"chrom": chroms, "sample": samples, "pos": positions, "pv_id": pv_ids})
pvs_df.to_csv("pvs_df.csv")
pvs_df

