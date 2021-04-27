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


# chroms = []
# samples = []
# positions = []
# pv_ids = []
# for header in recall_report.report["query_probe_header"]:
#     chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", header)[0]
#     chroms.append(chrom)
#     samples.append(sample)
#     positions.append(pos)
#     pv_ids.append(pv_id)
# pvs_df = pd.DataFrame(data={"chrom": chroms, "sample": samples, "pos": positions, "pv_id": pv_ids})
# pvs_df.to_csv("pvs_df.csv")
# pvs_df


# In[5]:


variation_found_nbofsamples = recall_report.get_proportion_of_allele_seqs_found_for_each_variant_with_nb_of_samples(
    binary=True
)
variation_found_nbofsamples.rename(columns={"proportion_of_allele_seqs_found_binary": "FOUND"}, inplace=True)
variation_found_nbofsamples.reset_index(inplace=True)
variation_found_nbofsamples


# In[6]:


Coordinate = namedtuple('Coordinate', ['sample', 'chrom', 'pos'])
pv_allele_seq_id_to_coordinates = defaultdict(set)

# hacky fix:
samples_to_reverse = "CFT073 Escherichia_coli_MINF_1D Escherichia_coli_MINF_8D Escherichia_coli_MSB1_1A Escherichia_coli_MSB1_3B Escherichia_coli_MSB1_4I Escherichia_coli_MSB1_6C Escherichia_coli_MSB1_7A Escherichia_coli_MSB1_7C Escherichia_coli_MSB1_8B Escherichia_coli_MSB1_8G Escherichia_coli_MSB1_9D H131800734".split()
sample_to_chrom_length = {'Escherichia_coli_MSB1_6C': 5230228,
 'Escherichia_coli_MSB1_7A': 5173732,
 'Escherichia_coli_MSB1_9D': 4963328,
 'CFT073': 5155066,
 'Escherichia_coli_MSB1_8G': 4827191,
 'Escherichia_coli_MSB1_1A': 4824417,
 'Escherichia_coli_MINF_7C': 5059386,
 'H131800734': 4725335,
 'Escherichia_coli_MINF_1D': 5092209,
 'Escherichia_coli_MINF_8D': 4843422,
 'Escherichia_coli_MSB1_4E': 5061196,
 'Escherichia_coli_MINF_1A': 5200778,
 'Escherichia_coli_MSB1_8B': 5251940,
 'Escherichia_coli_MSB2_1A': 4929229,
 'Escherichia_coli_MSB1_4I': 5278133,
 'Escherichia_coli_MSB1_7C': 5183046,
 '063_STEC': 4905296,
 'ST38': 5492930,
 'Escherichia_coli_MSB1_3B': 5173861,
 'Escherichia_coli_MINF_9A': 5121811}

for report_filepath in reports_filepaths:
    report_df = pd.read_csv(report_filepath, sep="\t")
    for query_probe_header in report_df.query_probe_header:
        chrom, sample, pos, pv_id, allele_seq_id = re.findall("CHROM=(.+?);SAMPLE=(.+?);POS=(.+?);.*;PANGENOME_VARIATION_ID=(.+?);.*;ALLELE_SEQUENCE_ID=(.+?);", query_probe_header)[0]
        if chrom=="0":  # TODO: hacky fix
            pos, pv_id, allele_seq_id = int(pos), int(pv_id), int(allele_seq_id)
            
            # TODO: hacky fix
            if sample in samples_to_reverse:
                pos = sample_to_chrom_length[sample]-pos-1
            
            pv_allele_seq_id_to_coordinates[(pv_id, allele_seq_id)].add(Coordinate(sample, chrom, pos))
pv_allele_seq_id_to_coordinates


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


# In[9]:


pv_allele_seq_id_to_genes_found_by_pandora = {}
for pv_allele_seq_id, gene_list in pv_allele_seq_id_to_genes.items():
    pv_allele_seq_id_to_genes_found_by_pandora[pv_allele_seq_id] = gene_list
pv_allele_seq_id_to_genes_found_by_pandora


pv_ids = set()
for pv_allele_seq_id in pv_allele_seq_id_to_genes.keys():
    pv_ids.add(pv_allele_seq_id[0])

    
pv_id_to_genes_found_by_pandora = {}
for pv_id in pv_ids:
    pv_id_to_genes_found_by_pandora[pv_id] =         pv_allele_seq_id_to_genes_found_by_pandora.get((pv_id, 0), list()) +         pv_allele_seq_id_to_genes_found_by_pandora.get((pv_id, 1), list())
pv_id_to_genes_found_by_pandora

pv_id_to_nb_of_genes_found_by_pandora = {}
for pv_id, genes in pv_id_to_genes_found_by_pandora.items():
    pv_id_to_nb_of_genes_found_by_pandora[pv_id] = len(genes)
pv_id_to_nb_of_genes_found_by_pandora
    
pv_id_to_nb_of_genes_found_by_pandora_df = pd.DataFrame(data=pv_id_to_nb_of_genes_found_by_pandora.items(),
                                                 columns=["PANGENOME_VARIATION_ID", "NB_OF_GENES"])
pv_id_to_nb_of_genes_found_by_pandora_df["FOUND_IN_VCF_REF"] = (pv_id_to_nb_of_genes_found_by_pandora_df["NB_OF_GENES"] > 0).astype(int)
pv_id_to_nb_of_genes_found_by_pandora_df


# In[10]:


dump_dict(pv_id_to_genes_found_by_pandora, "pv_id_to_genes_found_by_pandora")
pv_id_to_genes_found_by_pandora


# In[11]:


pv_ids_to_nb_of_genes_summary = pv_id_to_nb_of_genes_found_by_pandora_df[["PANGENOME_VARIATION_ID", "NB_OF_GENES"]].groupby("NB_OF_GENES").count()
plot = pv_ids_to_nb_of_genes_summary.plot(kind="bar")
fig = plot.get_figure()
fig.savefig("pandora_nb_of_genes_for_panvars.png")


# In[12]:


variation_found_nbofsamples = variation_found_nbofsamples.merge(pv_id_to_nb_of_genes_found_by_pandora_df)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==1)).astype(int)
variation_found_nbofsamples["NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"] = ((variation_found_nbofsamples["FOUND"]==0) & (variation_found_nbofsamples["FOUND_IN_VCF_REF"]==0)).astype(int)
variation_found_nbofsamples


# In[13]:


df = variation_found_nbofsamples[["NB_OF_SAMPLES", "FOUND", "NOT_FOUND_BY_PANDORA_BUT_IN_VCF_REF", "NOT_FOUND_BY_PANDORA_AND_NOT_IN_VCF_REF"]].groupby("NB_OF_SAMPLES").sum()
plot = df.plot(kind='bar', stacked=True)
fig = plot.get_figure()
fig.savefig("pandora_FN.png")


# with open("pv_id_to_genes_found_by_pandora.json") as fin:
#     json_lines = fin.readlines()[0]
# pv_id_to_genes_found_by_pandora_20_way = json.loads(json_lines)
# pv_id_to_genes_found_by_pandora_20_way

# gene_to_sequence = {}
# with open("data/pandora_multisample.vcf_ref.fa") as fin:
#     for line in fin:
#         line = line.strip()
#         if line[0] == ">":
#             gene=line[1:]
#         else:
#             seq=line
#             gene_to_sequence[gene] = seq
# gene_to_sequence

# import editdistance
# from collections import defaultdict
# 
# def edit_distance(gene_1, gene_2):
#     gene_1_seq = gene_to_sequence[gene_1]
#     gene_2_seq = gene_to_sequence[gene_2]
#     ed = editdistance.eval(gene_1_seq, gene_2_seq) / max(len(gene_1_seq), len(gene_2_seq))
#     if ed >= 0.75:
#         return "HIGH"
#     if ed >= 0.5:
#         return "MEDIUM"
#     return "LOW"
# 
# def ed_report(pv_id, genes):
#     ed_count = defaultdict(int)
#     for gene_1 in sorted(genes):
#         for gene_2 in sorted(genes):
#             if gene_1 < gene_2:
#                 ed_count[edit_distance(gene_1, gene_2)] += 1
#     high = ed_count["HIGH"]
#     medium = ed_count["MEDIUM"]
#     low = ed_count["LOW"]
#     ed_str = f"High: {high}, medium: {medium}, low: {low}"
#     print(f"PanVar: id = {pv_id}, nb of genes = {len(genes)}, ED: {ed_str}")
#     
#                 
# def make_ED_report_of_PVs_in_some_genes(nb_of_genes):
#     for pv_id, genes in pv_id_to_genes_found_by_pandora_20_way.items():
#         if len(genes)==nb_of_genes:
#             ed_report(pv_id, genes)
#             break

# for nb_of_genes in range(19, 1, -1):
#     make_ED_report_of_PVs_in_some_genes(nb_of_genes)
# 

# def get_gene_lengths(pv_id, genes):
#     gene_lengths = []
#     for gene in genes:
#         gene_seq = gene_to_sequence[gene]
#         gene_lengths.append(len(gene_seq))
#     return gene_lengths
# 
# def get_gene_lengths_of_PVs_in_some_genes(nb_of_genes):
#     for pv_id, genes in pv_id_to_genes_found_by_pandora_20_way.items():
#         if len(genes)==nb_of_genes:
#             gene_lengths = get_gene_lengths(pv_id, genes)
#             return gene_lengths
# 
# gene_length_array = []
# nb_of_gene_array = []
# for nb_of_genes in range(19, 1, -1):
#     gene_lengths = get_gene_lengths_of_PVs_in_some_genes(nb_of_genes)
#     gene_length_array.extend(gene_lengths)
#     nb_of_gene_array.extend([nb_of_genes] * len(gene_lengths))
# 
# import pandas as pd
# import seaborn as sns
# sns.set()
# sns.set(rc={'figure.figsize':(15,6)})
# df = pd.DataFrame(data={"gene_length": gene_length_array, "nb_of_genes": nb_of_gene_array})
# ax = sns.swarmplot(x="nb_of_genes", y="gene_length", data=df)

# def print_genes(genes):
#     for gene in genes:
#         print(f">{gene}")
#         print(gene_to_sequence[gene])
#     
# def get_seqs_of_PVs_in_some_genes(nb_of_genes):
#     for pv_id, genes in pv_id_to_genes_found_by_pandora_20_way.items():
#         if len(genes)==nb_of_genes:
#             print_genes(genes)
#             break
# 
# get_seqs_of_PVs_in_some_genes(4)

# import json
# with open("pv_id_to_genes_found_by_pandora_list.json") as fin:
#     json_lines = fin.readlines()[0]
# pv_id_to_genes_found_by_pandora_list_20_way = json.loads(json_lines)
# pv_id_to_genes_found_by_pandora_list_20_way

# from collections import defaultdict
# pv_id_to_genes_found_by_pandora_count_dict_20_way = defaultdict(dict)
# for pv_id, gene_list in pv_id_to_genes_found_by_pandora_list_20_way.items():
#     for gene in set(gene_list):
#         pv_id_to_genes_found_by_pandora_count_dict_20_way[pv_id][gene] = gene_list.count(gene)
# pv_id_to_genes_found_by_pandora_count_dict_20_way

# def get_seqs_of_PVs_in_some_genes():
#     total = 0
#     explained = 0
#     for pv_id, genes in pv_id_to_genes_found_by_pandora_count_dict_20_way.items():
#         counts_sorted = sorted(genes.values(), reverse=True)
#         if len(counts_sorted) <= 1:
#             explained += 1
#         else:
#             highest_freq, second_highest_freq = counts_sorted[:2]
#             if highest_freq >= second_highest_freq * 2:
#                 explained+=1
#         total += 1
#     print(f"Criteria explains {explained} out of {total} ({explained/total:.2f}%) PanVar gene mappings")
# get_seqs_of_PVs_in_some_genes()

# def edit_distance_core(gene_1, gene_2):
#     gene_1_seq = gene_to_sequence[gene_1]
#     gene_2_seq = gene_to_sequence[gene_2]
#     ed = editdistance.eval(gene_1_seq, gene_2_seq) / max(len(gene_1_seq), len(gene_2_seq))
#     return ed
# # edit_distance_core('GC00011350', 'Cluster_225')

# import pandas as pd
# pvs_df = pd.read_csv("pvs_df.csv")
# pvs_df = pvs_df[["chrom", "sample", "pos", "pv_id"]]
# pvs_df

# from Bio import SeqIO
# from collections import defaultdict
# from glob import glob
# ref_to_chrom_to_seq = defaultdict(lambda: defaultdict(str))
# refs = glob("data/samples/*.fa")
# for ref in refs:
#     for record in SeqIO.parse(ref, "fasta"):
#         ref = ref.replace("data/samples/", "").replace(".ref.fa", "")
#         ref_to_chrom_to_seq[ref][record.id] = record.seq
# ref_to_chrom_to_seq

# ref_to_chrom_length = {}
# for ref, chrom_to_seq in ref_to_chrom_to_seq.items():
#     chrom_seq = chrom_to_seq["0"]
#     ref_to_chrom_length[ref] = len(chrom_seq)
# ref_to_chrom_length

# def edit_distance_core(s1, s2):
#     ed = editdistance.eval(s1, s2) / max(len(s1), len(s2))
#     return ed
# the_ref = ref_to_chrom_to_seq["063_STEC"]["0"]
# for ref, chrom_to_seq in ref_to_chrom_to_seq.items():
#     chrom_seq = chrom_to_seq["0"]
#     edit_distance_fwd = edit_distance_core(chrom_seq, the_ref)
#     edit_distance_rc = edit_distance_core(chrom_seq, the_ref)
#     if edit_distance_rc < edit_distance_fwd:
#         print(f"Need to reverse: {ref}")

# import editdistance
# pv = pvs_df[pvs_df.pv_id==3]
# 
# flank_length=50
# probes = []
# ref_seq = None
# for chrom, sample, pos, pv_id in zip(pv["chrom"], pv["sample"], pv["pos"], pv["pv_id"]):
#     seq_fw = ref_to_chrom_to_seq[sample][str(chrom)]
#     seq_fw = seq_fw[pos-flank_length:pos+flank_length+1]
#     
#     if ref_seq == None:
#         ref_seq = seq_fw
#         seq = seq_fw
#     else:
#         seq_rc = ref_to_chrom_to_seq[sample][str(chrom)].reverse_complement()        
#         seq_rc = seq_rc[pos-flank_length:pos+flank_length+1]
#         if editdistance.eval(ref_seq, seq_fw) < editdistance.eval(ref_seq, seq_rc):
#             seq = seq_fw
#         else:
#             seq = seq_rc
#             print(f"Reverse {sample}")
#         
#     
#     probes.append((f">PV_ID={pv_id};Sample={sample};Chrom={chrom};Pos={pos};Strand=+",
#                   str(seq)))
# 
# 
#     
# # for header, seq in probes:
# #     print(header)
# #     print(seq)
# 
# # sample="CFT073"
# # chrom="0"
# # pos = 4413842
# # flank_length = 1000
# # chrom_seq = ref_to_chrom_to_seq[sample][chrom]
# # probe_fw = chrom_seq[pos-flank_length : pos+flank_length+1]
# # if probe_fw in chrom_seq:
# #     print("In FW")
# # else:
# #     print("In RC")
# 
# # # pos = len()-pos-1
# 
# # # print(ref_to_chrom_to_seq["CFT073"]["0"][pos-flank_length : pos+flank_length+1].reverse_complement())
