#!/usr/bin/env bash
set -eux

MEMORY=20000
LOCAL_CORES=20
PROFILE="lsf"
LOG_DIR=logs/
JOB_NAME="snakemake_master_process."$(date --iso-8601='minutes')

mkdir -p $LOG_DIR

bsub -R "select[mem>$MEMORY] rusage[mem=$MEMORY] span[hosts=1]" \
    -n "$LOCAL_CORES" \
    -M "$MEMORY" \
    -o "$LOG_DIR"/"$JOB_NAME".o \
    -e "$LOG_DIR"/"$JOB_NAME".e \
    -J "$JOB_NAME" \
      snakemake --local-cores "$LOCAL_CORES" --profile "$PROFILE" \
    --groups calculate_precision_per_sample_no_gt_conf=group_1 calculate_recall_per_sample_no_gt_conf_filter=group_2 calculate_recall_per_sample_vs_nb_of_samples=group_3 create_precision_report_from_probe_mappings=group_4 create_recall_report_for_truth_variants_mappings=group_5 create_recall_report_per_sample_for_calculator=group_6 filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_medaka=group_7 filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_nanopolish=group_8 filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools=group_9 filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy=group_10 fix_medaka_vcf_for_pipeline=group_11 fix_nanopolish_vcf_for_pipeline=group_12 fix_samtools_vcf_for_pipeline=group_13 fix_snippy_vcf_for_pipeline=group_14 make_mutated_vcf_ref_for_recall=group_15 make_variant_calls_probeset_for_precision=group_16 make_vcf_for_a_single_sample=group_17 map_recall_truth_probeset_to_mutated_vcf_ref=group_18 map_variant_call_probeset_to_reference_assembly=group_19 \
    --group-components group_1=10 group_2=10 group_3=10 group_4=10 group_5=10 group_6=5 group_7=10 group_8=10 group_9=10 group_10=10 group_11=10 group_12=10 group_13=10 group_14=10 group_15=10 group_16=10 group_17=10 group_18=10 group_19=10 \
    "$@"

exit 0
