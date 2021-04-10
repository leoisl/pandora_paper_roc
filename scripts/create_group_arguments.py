groups_and_nb_of_jobs = \
"""
calculate_precision_per_sample_no_gt_conf 10
calculate_recall_per_sample_no_gt_conf_filter 10
calculate_recall_per_sample_vs_nb_of_samples 10
create_precision_report_from_probe_mappings 10
create_recall_report_for_truth_variants_mappings 200
create_recall_report_per_sample_for_calculator 10
filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_medaka 10
filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_nanopolish 10
filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_samtools 10
filter_vcf_for_a_single_sample_by_gt_conf_percentile_for_snippy 10
fix_medaka_vcf_for_pipeline 10
fix_nanopolish_vcf_for_pipeline 10
fix_samtools_vcf_for_pipeline 10
fix_snippy_vcf_for_pipeline 10
gzip_vcf_file 100
index_gzipped_vcf_file 100
make_empty_depth_file 2000
make_mutated_vcf_ref_for_recall 10
make_variant_calls_probeset_for_precision 10
make_vcf_for_a_single_sample 10
map_recall_truth_probeset_to_mutated_vcf_ref 200
map_variant_call_probeset_to_reference_assembly 10
"""

for line in groups_and_nb_of_jobs.split("\n"):
    line_split = line.split()
    if len(line_split) == 2:
        job, nb_of_jobs = line_split
        group = f"group_{job}"
        print(f"--groups {job}={group} --group-components {group}={nb_of_jobs} \\")
