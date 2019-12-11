#!/usr/bin/env sh
unitigs=2


function clean_graph {
    mccortex="$1"  # path to mccortex container
    graph_to_clean="$2"  # graph to clean
    outpath="$3"
    kmersize="$4"

    singularity exec "$mccortex" mccortex31 clean \
        --tips=$(($kmersize * 2)) \
        --unitigs="$unitigs" \
        --force \
        --out "$outpath" \
        --covg-before ${outpath/.ctx}_covg_before.csv \
        --covg-after ${outpath/.ctx}_covg_after.csv \
        --len-before ${outpath/.ctx}_len_before.csv \
        --len-after ${outpath/.ctx}_len_after.csv \
        "$graph_to_clean"
}

function make_bandage_gfa {
    cortex_file="$1"
    fasta="$2"
    output_gfa="$3"

    singularity exec "$container" java -jar "$CORTEXJDK" ToGfa1 \
        --graph "$cortex_file" \
        --fasta "$fasta" \
        --output "$output_gfa"

}

container="$3"
graph_dir="$1"
out_dir="$2"
CORTEXJDK="$4"

for graph in $( realpath ${graph_dir}/*.ctx | grep -E '_K(7|9|11)\.ctx' )
do
    kmersize=$(echo $graph | grep -o -E '[0-9]+\.ctx' | grep -oE '[0-9]+')
    outpath=${out_dir}/$( basename ${graph/.ctx}_tips_cleaned_unitig${unitigs}_double_kmer.ctx )
    clean_graph "$container" "$graph" "$outpath" "$kmersize"
    fasta=${graph/.ctx}.fa
    gfa=${outpath/.ctx}.gfa
    make_bandage_gfa "$graph" "$fasta" "$gfa"
done
