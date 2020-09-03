rule merge_precision_and_recall_dfs:
    input:
         precision_file = output_folder + "/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv",
         recall_file = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    output:
         error_rate_and_recall_file = output_folder + "/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/merge_precision_and_recall_dfs/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.log"
    script:
        "../scripts/merge_precision_and_recall_dfs.py"


rule concat_all_plot_data:
    input:
         all_plot_data_intermediate_files = all_plot_data_intermediate_files
    output:
         final_plot_data_file = output_folder + "/plot_data/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_plot_data/ROC_data.log"
    script:
        "../scripts/concat_all_plot_data.py"


rule concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision:
    input:
         all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision = all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision
    output:
         nb_of_records_removed_with_mapq_sam_records_filter_for_precision_filepath = output_folder + "/plot_data/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.log"
    script:
        "../scripts/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files.py"


rule concat_all_recall_per_sample_no_gt_conf_filter:
    input:
         all_recall_per_sample_no_gt_conf_filter = all_recall_per_sample_no_gt_conf_filter
    output:
         recall_per_sample = output_folder + "/plot_data/recall_per_sample.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_recall_per_sample_no_gt_conf_filter/recall_per_sample.log"
    script:
        "../scripts/concat_all_recall_per_sample_no_gt_conf_filter.py"


rule concat_all_recall_per_sample_per_nb_of_samples:
    input:
         all_recalls_per_sample_per_nb_of_samples = cov_tool_and_filters_recall_per_sample_per_number_of_samples.values()
    output:
         aggregated_recall_per_sample_per_nb_of_samples = output_folder + "/plot_data/recall_per_sample_per_number_of_samples.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_recall_per_sample_per_nb_of_samples.log"
    run:
        import pandas as pd
        aggregated_df = pd.DataFrame(columns=["sample", "coverage", "tool", "coverage_threshold", "strand_bias_threshold",
                                   "gaps_threshold", "nb_of_samples", "recalls_wrt_truth_probes",
                                   "nbs_of_truth_probes_found", "nbs_of_truth_probes_in_total"])

        for (sample, coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold), df_filepath \
            in cov_tool_and_filters_recall_per_sample_per_number_of_samples.items():
            df = pd.read_csv(df_filepath)
            aggregated_df = pd.concat([aggregated_df, df], ignore_index=True)

        aggregated_df.to_csv(output.aggregated_recall_per_sample_per_nb_of_samples, index=False)


# rule concat_all_recall_per_sample_pair_no_gt_conf_filter:
#     input:
#          all_recall_per_sample_pair_no_gt_conf_filter = all_recall_per_sample_pair_no_gt_conf_filter
#     output:
#          recall_per_sample_pair = output_folder + "/plot_data/recall_per_sample_pair.tsv"
#     threads: 1
#     resources:
#         mem_mb = lambda wildcards, attempt: 4000 * attempt
#     log:
#         "logs/concat_all_recall_per_sample_pair_no_gt_conf_filter/recall_per_sample_pair.log"
#     script:
#         "../scripts/concat_all_recall_per_sample_pair_no_gt_conf_filter.py"


rule aggregate_recall_per_number_of_samples:
    input:
         all_recalls_per_number_of_samples = cov_tool_and_filters_recall_per_number_of_samples.values()
    output:
         aggregated_recall_per_number_of_samples = output_folder + "/plot_data/recall_per_number_of_samples.csv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/aggregate_recall_per_number_of_samples.log"
    run:
        import pandas as pd

        aggregated_df = pd.DataFrame(columns=["coverage", "tool", "coverage_threshold", "strand_bias_threshold",
                                   "gaps_threshold", "NB_OF_SAMPLES", "recall"])

        for (coverage, tool, coverage_threshold, strand_bias_threshold, gaps_threshold), df_filepath \
            in cov_tool_and_filters_recall_per_number_of_samples.items():
            df = pd.read_csv(df_filepath)
            df["coverage"] = coverage
            df["tool"] = tool
            df["coverage_threshold"] = coverage_threshold
            df["strand_bias_threshold"] = strand_bias_threshold
            df["gaps_threshold"] = gaps_threshold
            aggregated_df = pd.concat([aggregated_df, df], ignore_index=True)

        aggregated_df.to_csv(output.aggregated_recall_per_number_of_samples, index=False)


rule concat_all_precision_per_sample_no_gt_conf_filter:
    input:
         all_precision_per_sample_no_gt_conf_filter = all_precision_per_sample_no_gt_conf_filter
    output:
         precision_per_sample = output_folder + "/plot_data/precision_per_sample.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 4000 * attempt
    log:
        "logs/concat_all_precision_per_sample_no_gt_conf_filter.log"
    script:
        "../scripts/concat_all_precision_per_sample_no_gt_conf_filter.py"


rule make_enrichment_of_FPs_per_sample_plot:
    input:
        precision_per_sample = rules.concat_all_precision_per_sample_no_gt_conf_filter.output.precision_per_sample
    output:
        csv_data = output_folder + "/plot_data/enrichment_of_FPs/enrichment_of_FPs.csv",
        plot = output_folder + "/plot_data/enrichment_of_FPs/enrichment_of_FPs.png",
    threads: 1
    resources:
        mem_mb=lambda wildcards, attempt: 4000 * attempt
    notebook:
        "../../eda/enrichment_of_FPs_per_sample/enrichment_of_FPs_per_sample.ipynb"
