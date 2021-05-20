#!/usr/bin/env python
# coding: utf-8

# In[1]:


##################################################################################################################
# configs
# cluster
pandora_multisample_vcf_ref="/hps/nobackup/iqbal/leandro/pdrv2/paper_pandora2020_analyses/out_20_way/pandora_workflow/pandora_output_pandora_paper_tag1/nanopore/100x/random/compare_withdenovo/pandora_multisample.vcf_ref.fa"
msas_dir = "/hps/nobackup/iqbal/leandro/pdrv/data/msas"
##################################################################################################################


# In[2]:


# load pandora_multisample_vcf_ref sequences
from Bio import SeqIO

pandora_vcf_ref = {}
for record in SeqIO.parse(pandora_multisample_vcf_ref, "fasta"):
    pandora_vcf_ref[record.id] = str(record.seq)
pandora_vcf_ref


# In[3]:


# get the loci in the msas dir
from os import listdir
from os.path import isfile, join, getsize
from pathlib import Path

msa_locis = [f"{msas_dir}/{file}" for file in listdir(msas_dir) if file.endswith(".fa")]

# clean empty MSAs that sometime appear due to QC I guess
msa_locis = list(filter(lambda file: getsize(file)>0, msa_locis))

msa_locis = [Path(file).name[:-3] for file in msa_locis]
msa_locis


# In[4]:


loci_not_in_pandora_vcf_ref = set(msa_locis) - pandora_vcf_ref.keys()

assert len(loci_not_in_pandora_vcf_ref) + len(pandora_vcf_ref) == len(msa_locis)

loci_not_in_pandora_vcf_ref


# In[5]:


# get sequences of loci_not_in_pandora_vcf_ref
loci_not_in_pandora_vcf_ref_seqs = {}
for locus in loci_not_in_pandora_vcf_ref:
    fasta_file = f"{msas_dir}/{locus}.fa"
    for record in SeqIO.parse(fasta_file, "fasta"):
        loci_not_in_pandora_vcf_ref_seqs[f"{locus}_not_in_VCF_ref_only_in_PRG"] = str(record.seq.ungap())
        break
loci_not_in_pandora_vcf_ref_seqs


# In[6]:


vcf_ref_merged_with_PRG = dict(pandora_vcf_ref, **loci_not_in_pandora_vcf_ref_seqs)
assert len(vcf_ref_merged_with_PRG) == len(msa_locis)
vcf_ref_merged_with_PRG


# In[7]:


# output vcf_ref_merged_with_PRG to a fasta file
with open("vcf_ref_merged_with_PRG.fa", "w") as fout:
    for header, seq in vcf_ref_merged_with_PRG.items():
        print(f">{header}", file=fout)
        print(seq, file=fout)

