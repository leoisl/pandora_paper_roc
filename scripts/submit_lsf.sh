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
    --groups calculate_precision_per_sample_no_gt_conf=group_calculate_precision_per_sample_no_gt_conf --group-components group_calculate_precision_per_sample_no_gt_conf=10 \
    --groups calculate_recall_per_sample_no_gt_conf_filter=group_calculate_recall_per_sample_no_gt_conf_filter --group-components group_calculate_recall_per_sample_no_gt_conf_filter=10 \
    --groups calculate_recall_per_sample_vs_nb_of_samples=group_calculate_recall_per_sample_vs_nb_of_samples --group-components group_calculate_recall_per_sample_vs_nb_of_samples=10 \
    --groups create_precision_report_from_probe_mappings=group_create_precision_report_from_probe_mappings --group-components group_create_precision_report_from_probe_mappings=10 \
    --groups create_recall_report_for_truth_variants_mappings=group_create_recall_report_for_truth_variants_mappings --group-components group_create_recall_report_for_truth_variants_mappings=200 \
    --groups create_recall_report_per_sample_for_calculator=group_create_recall_report_per_sample_for_calculator --group-components group_create_recall_report_per_sample_for_calculator=10 \
    --groups filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_medaka=group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_medaka --group-components group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_medaka=10 \
    --groups filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_nanopolish=group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_nanopolish --group-components group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_nanopolish=10 \
    --groups filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools=group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools --group-components group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools=10 \
    --groups filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy=group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy --group-components group_filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy=10 \
    --groups fix_medaka_vcf_for_pipeline=group_fix_medaka_vcf_for_pipeline --group-components group_fix_medaka_vcf_for_pipeline=10 \
    --groups fix_nanopolish_vcf_for_pipeline=group_fix_nanopolish_vcf_for_pipeline --group-components group_fix_nanopolish_vcf_for_pipeline=10 \
    --groups fix_samtools_vcf_for_pipeline=group_fix_samtools_vcf_for_pipeline --group-components group_fix_samtools_vcf_for_pipeline=10 \
    --groups fix_snippy_vcf_for_pipeline=group_fix_snippy_vcf_for_pipeline --group-components group_fix_snippy_vcf_for_pipeline=10 \
    --groups gzip_vcf_file=group_gzip_vcf_file --group-components group_gzip_vcf_file=100 \
    --groups index_gzipped_vcf_file=group_index_gzipped_vcf_file --group-components group_index_gzipped_vcf_file=100 \
    --groups make_empty_depth_file=group_make_empty_depth_file --group-components group_make_empty_depth_file=2000 \
    --groups make_mutated_vcf_ref_for_recall=group_make_mutated_vcf_ref_for_recall --group-components group_make_mutated_vcf_ref_for_recall=10 \
    --groups make_variant_calls_probeset_for_precision=group_make_variant_calls_probeset_for_precision --group-components group_make_variant_calls_probeset_for_precision=10 \
    --groups make_vcf_for_a_single_sample=group_make_vcf_for_a_single_sample --group-components group_make_vcf_for_a_single_sample=10 \
    --groups map_recall_truth_probeset_to_mutated_vcf_ref=group_map_recall_truth_probeset_to_mutated_vcf_ref --group-components group_map_recall_truth_probeset_to_mutated_vcf_ref=200 \
    --groups map_variant_call_probeset_to_reference_assembly=group_map_variant_call_probeset_to_reference_assembly --group-components group_map_variant_call_probeset_to_reference_assembly=10 \
    "$@"

exit 0
