# If you want to use this container, it is simpler to just pull it:
#   docker pull leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1

# This container has basic software to be used by the pandora1_paper pipeline (mapping, indexing, compressing, etc)
# to build: sudo docker build . -t leandroishilima/pandora1_paper_basic_tools:pandora_paper_tag1
FROM continuumio/miniconda3:4.8.2
ENV PATH="/opt/conda/bin/:${PATH}"
RUN apt update && apt install -y gcc
RUN conda install -c conda-forge -c bioconda samtools==1.10 bwa==0.7.17 bcftools==1.10.2
ADD ./requirements.txt .
RUN pip install -r requirements.txt
