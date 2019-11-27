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


rule merge_all_plot_data:
    input:
         all_plot_data_intermediate_files = all_plot_data_intermediate_files
    output:
         final_plot_data_file =output_folder + "/plot_data/ROC_data.tsv"
    threads: 1
    resources:
        mem_mb = lambda wildcards, attempt: 2000 * attempt
    log:
        "logs/merge_all_plot_data/ROC_data.log"
    script:
        "../scripts/merge_all_plot_data.py"
