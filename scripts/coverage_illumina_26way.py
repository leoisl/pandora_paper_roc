coverages="""Escherichia_coli_MINF_1A/Escherichia_coli_MINF_1A.illumina.fastq.gz
21389648
Escherichia_coli_MINF_1D/Escherichia_coli_MINF_1D.illumina.fastq.gz
21244952
Escherichia_coli_MINF_2E/Escherichia_coli_MINF_2E.illumina.fastq.gz
20569224
Escherichia_coli_MINF_7C/Escherichia_coli_MINF_7C.illumina.fastq.gz
18358928
Escherichia_coli_MINF_8D/Escherichia_coli_MINF_8D.illumina.fastq.gz
20060880
Escherichia_coli_MINF_9A/Escherichia_coli_MINF_9A.illumina.fastq.gz
21274512
Escherichia_coli_MSB1_1A/Escherichia_coli_MSB1_1A.illumina.fastq.gz
20106392
Escherichia_coli_MSB1_3B/Escherichia_coli_MSB1_3B.illumina.fastq.gz
22315392
Escherichia_coli_MSB1_3C/Escherichia_coli_MSB1_3C.illumina.fastq.gz
31580472
Escherichia_coli_MSB1_3I/Escherichia_coli_MSB1_3I.illumina.fastq.gz
21140856
Escherichia_coli_MSB1_4E/Escherichia_coli_MSB1_4E.illumina.fastq.gz
27039936
Escherichia_coli_MSB1_4G/Escherichia_coli_MSB1_4G.illumina.fastq.gz
29564808
Escherichia_coli_MSB1_4I/Escherichia_coli_MSB1_4I.illumina.fastq.gz
21012808
Escherichia_coli_MSB1_6C/Escherichia_coli_MSB1_6C.illumina.fastq.gz
27150248
Escherichia_coli_MSB1_6J/Escherichia_coli_MSB1_6J.illumina.fastq.gz
20115480
Escherichia_coli_MSB1_7A/Escherichia_coli_MSB1_7A.illumina.fastq.gz
19143904
Escherichia_coli_MSB1_7C/Escherichia_coli_MSB1_7C.illumina.fastq.gz
28615672
Escherichia_coli_MSB1_8B/Escherichia_coli_MSB1_8B.illumina.fastq.gz
24981680
Escherichia_coli_MSB1_8G/Escherichia_coli_MSB1_8G.illumina.fastq.gz
20936696
Escherichia_coli_MSB1_9D/Escherichia_coli_MSB1_9D.illumina.fastq.gz
24079800
Escherichia_coli_MSB1_9I/Escherichia_coli_MSB1_9I.illumina.fastq.gz
20462576
Escherichia_coli_MSB2_1A/Escherichia_coli_MSB2_1A.illumina.fastq.gz
20636896
"""

coverages_as_list = coverages.split()
print(coverages_as_list)

for index, word in enumerate(coverages_as_list):
    if index%2==0:
        print(word)
    else:
        print(int(int(word)/4))