#!/usr/bin/env python
# coding: utf-8

# In[17]:


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
gene_localisation_dir="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_gene_distance/gene_distance_pandora_paper_tag1/genes_from_truth_or_ref/pandora_nanopore_100x_withdenovo"
pandora_multisample_vcf_ref="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_workflow/pandora_output_pandora_paper_tag1/nanopore/100x/random/compare_withdenovo/pandora_multisample.vcf_ref.fa"
##################################################################################################################


# In[ ]:


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


# In[ ]:


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
pvs_df = pd.DataFrame(data={"chrom": chroms_array, "sample": samples_array, "pos": positions, "pv_id": pv_ids_array})
pvs_df.to_csv("pvs_df.csv", index=False)
pvs_df


# In[ ]:


sample_to_chrom_to_intervaltree = defaultdict(lambda: defaultdict(lambda: intervaltree.IntervalTree()))
for sample in samples:
    filename = gene_localisation_dir+f"/{sample}.csv"
    with open(filename) as gene_localisation_fh:
        for line in gene_localisation_fh:
            status,gene_name,ref_or_truth_id,contig,start,stop,sequence,strand = line.strip().split(",")
            if status=="Mapped":
                start, stop = int(start), int(stop)
                sample_to_chrom_to_intervaltree[ref_or_truth_id][contig][start:stop] = gene_name

sample_to_chrom_to_intervaltree


# In[ ]:


def find_genes(chrom, sample, pos):
    genes = set()
    for gene in sample_to_chrom_to_intervaltree[sample][chrom][pos]:
        genes.add(gene.data)
    return genes

pv_to_genes = defaultdict(set)
for chrom, sample, pos, pv_id in zip(pvs_df["chrom"], pvs_df["sample"], pvs_df["pos"], pvs_df["pv_id"]):
    pv_to_genes[pv_id].update(find_genes(chrom, sample, pos))

for pv in pv_to_genes:
    pv_to_genes[pv] = list(pv_to_genes[pv])

dump_object_to_json(pv_to_genes, "pv_to_genes")

pv_to_genes


# In[ ]:


variation_found_nbofsamples = recall_report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(
    binary=True
)
variation_found_nbofsamples.rename(columns={"proportion_of_allele_seqs_found_binary": "FOUND"}, inplace=True)
variation_found_nbofsamples.reset_index(inplace=True)
variation_found_nbofsamples.to_csv("variation_found_nbofsamples.csv", index=False)
variation_found_nbofsamples


# In[ ]:


pv_to_nb_of_genes_found_by_pandora = {}
for pv, genes in pv_to_genes.items():
    pv_to_nb_of_genes_found_by_pandora[pv] = len(genes)
    
pv_id_to_nb_of_genes_found_by_pandora_df = pd.DataFrame(data=pv_to_nb_of_genes_found_by_pandora.items(),
                                                 columns=["PVID", "NB_OF_GENES"])
pv_id_to_nb_of_genes_found_by_pandora_df["FOUND_IN_VCF_REF"] = (pv_id_to_nb_of_genes_found_by_pandora_df["NB_OF_GENES"] > 0).astype(int)

pv_id_to_nb_of_genes_found_by_pandora_df.to_csv("pv_id_to_nb_of_genes_found_by_pandora_df.csv", index=False)

pv_id_to_nb_of_genes_found_by_pandora_df


# In[ ]:


pv_ids_to_nb_of_genes_summary = pv_id_to_nb_of_genes_found_by_pandora_df[["PVID", "NB_OF_GENES"]].groupby("NB_OF_GENES").count()
plot = pv_ids_to_nb_of_genes_summary.plot(kind="bar")
fig = plot.get_figure()
fig.savefig("pandora_nb_of_genes_for_panvars.png")


# In[ ]:


variation_found_nbofsamples = variation_found_nbofsamples.merge(pv_id_to_nb_of_genes_found_by_pandora_df)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==1)).astype(int)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==0)).astype(int)
variation_found_nbofsamples.to_csv("variation_found_nbofsamples.csv", index=False)
variation_found_nbofsamples


# In[ ]:


df = variation_found_nbofsamples[["NB_OF_SAMPLES", "FOUND", "NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF", "NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"]].groupby("NB_OF_SAMPLES").sum()
plot = df.plot(kind='bar', stacked=True)
fig = plot.get_figure()
fig.savefig("pandora_FN.png")

