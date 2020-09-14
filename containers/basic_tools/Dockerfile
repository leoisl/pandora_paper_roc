# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1

# This container has basic software to be used by the pandora1_paper pipeline (mapping, indexing, compressing, etc)
# to build: sudo docker build . -t leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1
FROM biocontainers/biocontainers@sha256:c6f082f5658c6d31a5ccb537289d1ce777f51b7285c465fe0efa219058471c81

RUN conda install -c conda-forge -c bioconda samtools==1.10 bwa==0.7.17 bcftools==1.10.2
