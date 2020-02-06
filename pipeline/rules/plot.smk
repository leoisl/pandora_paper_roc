rule merge_precision_and_recall_dfs:
    input:
         precision_file = output_folder + "/precision/precision_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/precision.tsv",
         recall_file = output_folder + "/recall/recall_files/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/recall.tsv"
    output:
         error_rate_and_recall_file = output_folder + "/plot_data/{coverage}/{tool}/coverage_filter_{coverage_threshold}/strand_bias_filter_{strand_bias_threshold}/gaps_filter_{gaps_threshold}/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
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
        mem_mb = lambda wildcards, attempt: 2000 * attempt
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
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files_for_precision/nb_of_records_removed_with_mapq_sam_records_filter_for_precision.log"
    script:
        "../scripts/concat_all_nb_of_records_removed_with_mapq_sam_records_filter_files.py"
