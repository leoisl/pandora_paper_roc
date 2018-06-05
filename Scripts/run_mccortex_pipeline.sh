#!/usr/bin/env bash
FASTQ="$1"
OUTDIR="$2"
MEMORY="${5:-1}"
MCCORTEX_CONTAINER="$3"
CORTEXJDK="$4"

function run_mccortex() {
    memory="${4:-1}"
    fastq_sequence="$1"
    output_cortex="${2}".ctx
    output_fasta="${2}".fa
    k="$3"
    sample_name="${fastq_sequence##*.}"_K"${k}"

    # build mccortex graph with given k
    singularity exec "$MCCORTEX_CONTAINER" mccortex31 build \
        --force \
        --memory ${memory}G \
        --kmer "$k" \
        --sample "$sample_name" \
        --seq "$fastq_sequence" \
        "$output_cortex"

    # get unitigs in gfa format
    # singularity exec "$MCCORTEX_CONTAINER" mccortex31 unitigs --gfa "$output_cortex" > ${2}.gfa

    # get unitigs in fasta format
    singularity exec "$MCCORTEX_CONTAINER" mccortex31 unitigs \
        --fasta "$output_cortex" > "$output_fasta"

}

function make_bandage_gfa() {
    cortex_file=${1}.ctx
    fasta=${1}.fa
    output_gfa=${1}.gfa

    singularity exec "$MCCORTEX_CONTAINER" java -jar "$CORTEXJDK" ToGfa1 \
        --graph "$cortex_file" \
        --fasta "$fasta" \
        --output "$output_gfa"

}

function clean_graph {
    mccortex="$1"  # path to mccortex container
    graph_to_clean="$2"  # graph to clean
    outpath="$3"
    kmersize="$4"

    singularity exec "$mccortex" mccortex31 clean \
        --tips=$(($kmersize * 2)) \
        --unitigs=2 \
        --force \
        --out "$outpath" \
        "$graph_to_clean"
}

# ===========================
# MAIN CODE
# build mccortex graph and gfa for bandage for each k between 3 and 15 (uneven)
# ===========================
for K in {7..11..2}
do
    file_prefix=$(basename "${FASTQ/.fa*}")
    output_prefix=${OUTDIR}/${file_prefix}_mccortex_K${K}

    run_mccortex "$FASTQ" "$output_prefix" "$K"
    clean_graph "$MCCORTEX_CONTAINER" ${output_prefix}.ctx \
        ${outside_prefix}_clean.ctx "$K"
        
    make_bandage_gfa "$output_prefix"
done
