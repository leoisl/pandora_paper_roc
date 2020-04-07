#!/usr/bin/env bash
set -eu

if [ $# -ne 4 ]
then
    echo "usage: $0 assembly.fa reads_1.fq reads_2.fq outprefix
Makes a BED file of bad regions, by mapping reads and running bam_to_low_qual_mask.pl.
BWA and samtools faidx indexes the assembly fasta if needed"
    exit
fi

assembly="$1"
reads_1="$2"
reads_2="$3"
outprefix="$4"
bam="$outprefix".mask.bam
bed="$outprefix".mask.bed

if [ ! -f "${assembly}.fai" ]; then
    samtools faidx "$assembly"
fi

if [ ! -f "${assembly}.bwt" ]; then
    bwa index "$assembly"
fi

bwa mem -x intractg "$assembly" "$reads_1" "$reads_2" | samtools sort -o "$bam"
./scripts/bam_to_low_qual_mask.pl "$bam" "$bed"

#belongs to martin
#https://github.com/martinghunt/bioinf-scripts/blob/master/bash/make_low_qual_genome_mask.sh
